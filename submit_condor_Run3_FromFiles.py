import os
import yaml
import subprocess
import argparse
import math

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Prepare Condor jobs for HLT rerun")
parser.add_argument("--jobFlavour", type=str, default="tomorrow",
                    help="Condor jobflavor (e.g., espresso, longlunch, workday, tomorrow)")
parser.add_argument("--n", type=int, default=10, help="Number of DAS files per shell script")
parser.add_argument("--farm", type=str, default="Farm",
                    help="farm directory")
parser.add_argument("--config", type=str)
parser.add_argument("--fileList", type=str, required=True, help="Text file containing one input ROOT file per line")
parser.add_argument("--proxy", type=str, default=None)
args = parser.parse_args()
if args.n < 1:
    raise ValueError("--n must be >= 1")

# ------------------------------------------------------------
# Load YAML
# ------------------------------------------------------------
with open(args.config) as f:
    cfg = yaml.safe_load(f)

cmssw_dir = cfg["cmssw-dir"]
eos_dir = cfg["eos-dir"]

shell_scripts = []

# ------------------------------------------------------------
# Loop over projects
# ------------------------------------------------------------

for project_name, project_cfg in cfg["project"].items():
    project_eos_dir = os.path.join(eos_dir, project_name)
    os.makedirs(project_eos_dir, exist_ok=True)

    farm_dir = os.path.join(args.farm, project_name)
    os.makedirs(farm_dir, exist_ok=True)
    
    # Expand blocks into files
    all_input_files = []
    with open(args.fileList) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Allow either /eos/cms/... or file:/eos/cms/...
            if not line.startswith("file:"):
                line = "file:" + line

            all_input_files.append(line)

    if len(all_input_files) == 0:
        raise RuntimeError("No input files found.")

    print(f"Loaded {len(all_input_files)} input files.")

    # Output file listing
    runlist_file = os.path.join(farm_dir, "files_to_run.txt")

    with open(runlist_file, "w") as runlist:
        # Split files into chunks of N
        file_chunks = [
            all_input_files[i:i + args.n]
            for i in range(0, len(all_input_files), args.n)
        ]
        print(f"Creating {len(file_chunks)} jobs")

        for job_idx, chunk in enumerate(file_chunks):

            script_name = os.path.join(farm_dir, f"run_{project_name}_{job_idx}.sh")
            workspace = f"$TMPDIR/Job_{project_name}_{job_idx}"

            print(f"creating {script_name}")
            for infile in chunk:
                runlist.write(f"{infile}\n")

            filelist = ",".join(chunk)

            with open(script_name, "w") as f:

                f.write("#!/bin/bash\n")
                f.write("set -e\n\n")

                f.write("export X509_USER_PROXY=$1\n")

                f.write(f"WORKDIR={workspace}\n")
                f.write("mkdir -p $WORKDIR\n")
                f.write("echo 'Using TMPDIR=' $WORKDIR\n\n")

                f.write(f"cd {cmssw_dir}\n")
                f.write("eval `scramv1 runtime -sh`\n\n")

                f.write("cd $WORKDIR\n\n")

                for i, cmd in enumerate(project_cfg["process"]):

                    cmd_tmp = cmd.replace("file:", "file:$WORKDIR/")

                    if "TnPTreeProducer" in cmd_tmp:
                        cmd_tmp = cmd_tmp.replace('outputFile=file:', 'outputFile=')

                    if "python3" in cmd_tmp:
                        cmd_tmp = cmd_tmp.replace('output_Phase2_HLT.root', '$WORKDIR/output_Phase2_HLT.root')
                        cmd_tmp = cmd_tmp.replace('tnpNtupler.root', '$WORKDIR/tnpNtupler.root')

                    if i == 0:
                        cmd_tmp += f" --filein {filelist} "

                    f.write(f"echo 'Running step {i+1}'\n")
                    f.write(f"{cmd_tmp}\n\n")

                # Copy outputs to EOS    
                for storefile in project_cfg.get('storefile', []):

                    eos_path = os.path.join(
                        project_eos_dir,
                        storefile.replace(".root", f"_{job_idx}.root")
                    )

                    f.write(
                        f"xrdcp $WORKDIR/{storefile} root://eosuser.cern.ch/{eos_path}\n"
                    )

                f.write("rm -f $WORKDIR/*.root\n")

            os.chmod(script_name, 0o755)
            shell_scripts.append(script_name)

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

# Request enough resources
#condor_str += "request_disk = 80000\n"
#condor_str += "request_memory = 4000\n"

condor_str += f"queue filename matching ({args.farm}/*/*.sh)"

with open(sub_file, "w") as condor_file:
    condor_file.write(condor_str)
