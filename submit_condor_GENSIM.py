import os
import yaml
import argparse

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


# -------------------------------
# MAIN LOOP
# -------------------------------

for project_name, project_cfg in cfg["project"].items():
    project_eos_dir = os.path.join(eos_dir, project_name)
    os.makedirs(project_eos_dir, exist_ok=True)

    farm_dir = os.path.join(args.farm, project_name)
    os.makedirs(farm_dir, exist_ok=True)

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

            # Stay in the CMSSW area for scram/cmsDriver/cmsRun. Do NOT cd to $WORKDIR first:
            # scram/cmsDriver must see the release (SCRAM fails in /tmp).
            # Output ROOT files use file:$WORKDIR/... from the YAML substitutions below.
            f.write(f"cd {cmssw_dir}\n")
            f.write("eval `scramv1 runtime -sh`\n\n")

            for i, cmd in enumerate(project_cfg["process"]):
                cmd_tmp = substitute_job_tokens(cmd, job_idx)
                cmd_tmp = cmd_tmp.replace("file:", "file:$WORKDIR/")

                if "TnPTreeProducer" in cmd_tmp:
                    cmd_tmp = cmd_tmp.replace("outputFile=file:", "outputFile=")

                if "python3" in cmd_tmp:
                    cmd_tmp = cmd_tmp.replace(
                        "output_Phase2_HLT.root", "$WORKDIR/output_Phase2_HLT.root"
                    )
                    cmd_tmp = cmd_tmp.replace("tnpNtupler.root", "$WORKDIR/tnpNtupler.root")

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

condor_str += f"queue filename matching ({args.farm}/*/*.sh)\n"

with open(sub_file, "w") as condor_file:
    condor_file.write(condor_str)

print("\n✅ Condor submission file created:", sub_file)
print(f"✅ Total jobs created: {len(shell_scripts)}")
print(f"   (--n-jobs {n_jobs} per project × {len(cfg['project'])} project(s), {n_events} events per job)")
print("   (YAML placeholders supported: __NEVENTS__, __JOBID__)")
