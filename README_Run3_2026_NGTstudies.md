## Relevant Links:
- CMSSW: `CMSSW_16_0_7_patch1`
- Menu used: `/dev/CMSSW_16_0_0/GRun/V63`
- Runs used: [403751, 403774, 403775, 403811, 403818, 403867, 403894, 403937]
- Dataset used: 
	- `/EGamma0/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma1/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma2/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma3/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`
    - `/EGamma4/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`
    - `/EGamma5/Run2026D-ZElectron-PromptReco-v1/RAW-RECO`

- HLT GT: 
	- Reference: `160X_dataRun3_HLT_v1`
	- Target: `160X_dataRun3_NGT_v4`

- Prompt GT for RECO step: `160X_dataRun3_Prompt_v1`

### Rucio rules for dataset if not available on disk
```
source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
voms-proxy-init -voms cms
export RUCIO_ACCOUNT=`whoami`
rucio add-rule cms:/EGamma3/Run2026D-ZElectron-PromptReco-v1/RAW-RECO 1 T2_CH_CERN --lifetime 259200 --comment "For urgent NGT studies"
```

### CMSSW set-up
```
cmsrel CMSSW_16_0_7_patch1; cd CMSSW_16_0_7_patch1/src/; cmsenv;
git cms-init 
git cms-addpkg HLTrigger/Configuration
scram b -j 10
```

### hltGetConfiguration for re-HLT step
```
### hltGetConfiguration
hltGetConfiguration /dev/CMSSW_16_0_0/GRun/V63 --output minimal --data --process MYHLT --type GRun --globaltag 160X_dataRun3_HLT_v1 --max-events 10 --unprescale --eras Run3_2026 --path HLTriggerFirstPath,HLT_Ele30_WPTight_Gsf_v*,HLT_Ele32_WPTight_Gsf_v*,HLT_Ele115_CaloIdVT_GsfTrkIdT_v*,HLT_Ele135_CaloIdVT_GsfTrkIdT_v*,HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*,HLT_DoubleEle33_CaloIdL_MW_v*,HLTriggerFinalPath --cff > "${CMSSW_BASE}"/src/HLTrigger/Configuration/python/HLT_2026D_cff.py
```

### Submit jobs on condor 
Use `submit_condor_Run3_FromFiles.py` to submit the condor jobs for `config/Run3_2026D_NGT_Studies.yaml`. 
All the commands/steps/drivers, input-output directories for both target (NGT) and reference (HLT here) are controlled using the config file `config/Run3_2026D_NGT_Studies.yaml`.
`etc/find_das_files.py` can be run to produce the text file with the list of all the input EGamma RAW-RECO files from the respective datasets and run numbers. 
The jobs are then submitted over the rootfiles in `etc/files_to_run_2026_NGT.txt`.
```
git clone git@github.com:saumyaphor4252/EGamma-Workflow.git
git clone -b 160X_2026_Studies git@github.com:saumyaphor4252/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer
scram b -j 10
voms-proxy-init --valid 100:00
cp /tmp/x509up_u<999999> /afs/cern.ch/user/s/ssaumya/private/x509up_u<999999>
python3 submit_condor_Run3_FromFiles.py --n 10 --config config/Run3_2026D_NGT_Studies.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999> --farm Jobs --fileList etc/files_to_run_2026_NGT.txt
```

### Plotting

### Run Iason's tool for plotting

#### First time setup
If you running the tool for the first time: 
```
ssh ssaumya@lxplus9.cern.ch -L8787:localhost:8787
cd /afs/cern.ch/work/s/ssaumya/private/Egamma/IasonTool/
git clone -b 2024_Studies git@github.com:saumyaphor4252/egamma-tnp.git
cd egamma-tnp/
python3 -m venv egmtnpenv
source egmtnpenv/bin/activate
pip install jupyter
pip install ipython
pip install ipykernel
ipython kernel install --user --name=egmtnpenv
python -m ipykernel install --user --name=egmtnpenv
pip install . --no-cache-dir
jupyter lab --no-browser --port 8787
```

#### After initial set-up 
Once the set-up is done, once can use the below commands from next time:
```
ssh -o ServerAliveInterval=10 ssaumya@lxplus9.cern.ch -L8787:localhost:8787
cd /afs/cern.ch/work/s/ssaumya/private/Egamma/IasonTool/egamma-tnp/
source egmtnpenv/bin/activate
jupyter lab --no-browser --port 8787

# Code template: https://github.com/saumyaphor4252/egamma-tnp/tree/2025_Studies
# Notebook template: STEAM_23ndMay2025.ipynb
```
