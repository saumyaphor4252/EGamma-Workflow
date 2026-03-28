import os
import re
import yaml
import argparse
import subprocess

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Prepare Condor jobs from a placeholder-based YAML config.")
parser.add_argument(
    "--jobFlavour",
    type=str,
    default="tomorrow",
    help="Condor jobflavor (e.g., espresso, longlunch, workday, tomorrow)",
)
parser.add_argument(
    "--n-jobs",
    type=int,
    required=True,
    help="Number of shell scripts / Condor jobs to create per project",
)
parser.add_argument(
    "--events-per-job",
    type=int,
    required=True,
    help="Number of events to run in each job (replaces __NEVENTS__ in YAML commands)",
)
parser.add_argument("--farm", type=str, default="Farm", help="Farm directory (under which project subdirs are created)")
parser.add_argument("--config", type=str, required=True, help="YAML config with cmssw-dir, eos-dir, project.*")
parser.add_argument("--proxy", type=str, default=None, help="Path to grid proxy for Condor arguments")
args = parser.parse_args()

n_events = args.events_per_job
n_jobs = args.n_jobs

# Load YAML config
with open(args.config) as f:
    cfg = yaml.safe_load(f)

cmssw_dir = os.path.normpath(cfg["cmssw-dir"].rstrip("/"))
eos_dir = cfg["eos-dir"]

shell_scripts = []


def substitute_job_tokens(cmd: str, job_idx: int) -> str:
    """Replace explicit placeholders only."""
    out = cmd
    out = out.replace("__NEVENTS__", str(n_events))
    out = out.replace("__JOBID__", str(job_idx))
    return out


def per_job_cfg_in_workdir(job_idx: int, cfg_basename: str) -> str:
    """Unique cfg path under $WORKDIR (expanded on worker) — avoids AFS races on shared names."""
    return f"$WORKDIR/j{job_idx}_{cfg_basename}"


def rewrite_cmsdriver_python_filename(cmd: str, job_idx: int) -> str:
    """Put --python_filename under $WORKDIR with a job-specific prefix (cmsCondor-style isolation)."""

    def repl(m):
        name = m.group(1)
        if name.startswith("$WORKDIR/") or name.startswith("/"):
            return m.group(0)
        return f"--python_filename {per_job_cfg_in_workdir(job_idx, name)}"

    return re.sub(r"--python_filename\s+(\S+\.py)", repl, cmd)


def wrap_cmsrun_for_workdir(cmd: str, cmssw_dir: str, job_idx: int) -> str:
    """
    cmsDriver --fileout file:$WORKDIR/... is not always embedded as an absolute path in the
    generated cfg; PoolOutput may write relative to cwd. Run cmsRun from $WORKDIR with the cfg
    in $WORKDIR (same basename as --python_filename after rewrite_cmsdriver_python_filename).
    """
    stripped = cmd.strip()
    if not stripped.startswith("cmsRun "):
        return cmd
    parts = stripped.split()
    if len(parts) < 2:
        return cmd
    cfg = parts[1]
    if cfg.startswith("/") or cfg.startswith("$WORKDIR/"):
        cfg_path = cfg
    elif cfg.endswith(".py"):
        cfg_path = per_job_cfg_in_workdir(job_idx, cfg)
    else:
        cfg_path = os.path.join(cmssw_dir, cfg)
    rest = parts[2:]
    inner = " ".join(["cmsRun", cfg_path] + rest)
    return f"cd $WORKDIR && {inner}"


def run_prepare_once(commands, cmssw_dir: str) -> None:
    """
    Run heavy/shared setup commands once (outside Condor jobs), e.g. curl + scram b.
    This avoids many jobs racing in the same CMSSW area.
    """
    if not commands:
        return
    print(f"\n🔧 Running {len(commands)} prepare-once command(s) in {cmssw_dir}")
    setup = f"cd {cmssw_dir} && eval `scramv1 runtime -sh` && "
    for idx, cmd in enumerate(commands, start=1):
        full_cmd = setup + cmd
        print(f"   [{idx}/{len(commands)}] {cmd}")
        subprocess.run(full_cmd, shell=True, check=True, executable="/bin/bash")


# -------------------------------
# MAIN LOOP
# -------------------------------

for project_name, project_cfg in cfg["project"].items():
    project_eos_dir = os.path.join(eos_dir, project_name)
    os.makedirs(project_eos_dir, exist_ok=True)

    farm_dir = os.path.join(args.farm, project_name)
    os.makedirs(farm_dir, exist_ok=True)

    process_commands = list(project_cfg["process"])
    prepare_once_cmds = []
    filtered_process_cmds = []
    for cmd in process_commands:
        stripped = cmd.strip()
        if stripped.startswith("curl ") or stripped.startswith("scram b"):
            prepare_once_cmds.append(cmd)
        else:
            filtered_process_cmds.append(cmd)

    # Execute shared build/setup commands once to avoid per-job races.
    run_prepare_once(prepare_once_cmds, cmssw_dir)

    for job_idx in range(n_jobs):
        script_name = os.path.join(farm_dir, f"run_{project_name}_{job_idx}.sh")
        workspace = f"$TMPDIR/Job_{project_name}_{job_idx}"

        print(f"Creating {script_name}")

        with open(script_name, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("set -e\n\n")

            f.write("export X509_USER_PROXY=$1\n")

            f.write(f"WORKDIR={workspace}\n")
            f.write("mkdir -p $WORKDIR\n")
            f.write("echo 'Using TMPDIR=' $WORKDIR\n\n")

            # Stay in the CMSSW area for scram/cmsDriver. cmsRun is wrapped to `cd $WORKDIR &&`
            # so PoolOutput paths match file:$WORKDIR/... (see wrap_cmsrun_for_workdir).
            f.write(f"cd {cmssw_dir}\n")
            f.write("eval `scramv1 runtime -sh`\n\n")

            for i, cmd in enumerate(filtered_process_cmds):
                cmd_tmp = substitute_job_tokens(cmd, job_idx)
                cmd_tmp = cmd_tmp.replace("file:", "file:$WORKDIR/")

                if "TnPTreeProducer" in cmd_tmp:
                    cmd_tmp = cmd_tmp.replace("outputFile=file:", "outputFile=")

                if "python3" in cmd_tmp:
                    cmd_tmp = cmd_tmp.replace(
                        "output_Phase2_HLT.root", "$WORKDIR/output_Phase2_HLT.root"
                    )
                    cmd_tmp = cmd_tmp.replace("tnpNtupler.root", "$WORKDIR/tnpNtupler.root")

                cmd_tmp = rewrite_cmsdriver_python_filename(cmd_tmp, job_idx)
                cmd_tmp = wrap_cmsrun_for_workdir(cmd_tmp, cmssw_dir, job_idx)

                f.write(f"echo 'Running step {i+1}'\n")
                f.write(f"{cmd_tmp}\n\n")

            for storefile in project_cfg.get("storefile", []):
                eos_path = os.path.join(
                    project_eos_dir,
                    storefile.replace(".root", f"_{job_idx}.root"),
                )

                f.write(f"xrdcp $WORKDIR/{storefile} root://eosuser.cern.ch/{eos_path}\n")

            f.write("rm -f $WORKDIR/*.root\n")

        os.chmod(script_name, 0o755)
        shell_scripts.append(script_name)

# -------------------------------
# CONDOR SUBMIT FILE
# -------------------------------
sub_file = f"{args.farm}/condor_jobs.sub"

condor_str = ""
condor_str += "executable = $(filename)\n"

if args.proxy is not None:
    condor_str += f"Proxy_path = {args.proxy}\n"
    condor_str += "arguments = $(Proxy_path)\n"

condor_str += "output = $Fp(filename)$(filename)_hlt.stdout\n"
condor_str += "error = $Fp(filename)$(filename)_hlt.stderr\n"
condor_str += "log = $Fp(filename)$(filename)_hlt.log\n"

condor_str += f'+JobFlavour = "{args.jobFlavour}"\n'

# Disk: bare KB is easy to under-request; Phase2 chain needs many GB on $TMPDIR.
condor_str += "request_disk = 25 GB\n"
condor_str += "request_memory = 8000\n"

condor_str += f"queue filename matching ({args.farm}/*/*.sh)\n"

with open(sub_file, "w") as condor_file:
    condor_file.write(condor_str)

print("\n✅ Condor submission file created:", sub_file)
print(f"✅ Total jobs created: {len(shell_scripts)}")
print(f"   (--n-jobs {n_jobs} per project × {len(cfg['project'])} project(s), {n_events} events per job)")
print("   (YAML placeholders supported: __NEVENTS__, __JOBID__)")
