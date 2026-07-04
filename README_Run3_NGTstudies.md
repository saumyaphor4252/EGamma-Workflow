## Relevant Links:
- CMSSW: `CMSSW_15_0_15_patch4`
- Menu used: `/dev/CMSSW_15_0_0/GRun/V119`
- Runs used: [398675, 398680, 398681, 398682, 398683, 398787, 398797, 398801, 398802, 398803, 398827, 398828, 398858]
- Dataset used: 
	- `/EGamma0/Run2025G-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma1/Run2025G-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma2/Run2025G-ZElectron-PromptReco-v1/RAW-RECO`
	- `/EGamma3/Run2025G-ZElectron-PromptReco-v1/RAW-RECO`

- HLT GT: 
	- Reference: `150X_dataRun3_HLT_v1`
	- Target: `150X_dataRun3_NGT_v2`

### Rucio rules for dataset if not available on disk
```
source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
voms-proxy-init -voms cms
export RUCIO_ACCOUNT=`whoami`
rucio add-rule cms:/EGamma0/Run2025G-ZElectron-PromptReco-v1/RAW-RECO#7a8aa189-ef77-4039-9ba0-522527094ea6 1 T2_CH_CERN --lifetime 259200 --comment "For urgent NGT studies"
```

### CMSSW set-up
```
cmsrel CMSSW_15_0_15_patch4; cd CMSSW_15_0_15_patch4/src/; cmsenv;
git cms-init 
git cms-addpkg HLTrigger/Configuration
scram b -j 10
```


```
### hltGetConfiguration
hltGetConfiguration /dev/CMSSW_15_0_0/GRun/V119 --output minimal --data --process MYHLT --type GRun --globaltag 150X_dataRun3_HLT_v1 --max-events 100 --unprescale --eras Run3_2025 --path HLTriggerFirstPath,HLT_Ele30_WPTight_Gsf_v*,HLT_Ele32_WPTight_Gsf_v*,HLT_Ele115_CaloIdVT_GsfTrkIdT_v*,HLT_Ele135_CaloIdVT_GsfTrkIdT_v*,HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*,HLT_DoubleEle33_CaloIdL_MW_v*,HLTriggerFinalPath --cff > "${CMSSW_BASE}"/src/HLTrigger/Configuration/python/HLT_2025G_cff.py

```

### Submit jobs on condor 
Use `submit_condor_Run3_FromFiles.py` to submit the condor jobs for `config/Run3_2025G_NGT_Studies.yaml`
```
git clone git@github.com:saumyaphor4252/EGamma-Workflow.git
git clone -b 150X_2025_Studies git@github.com:saumyaphor4252/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer
scram b -j 10
voms-proxy-init --valid 100:00
cp /tmp/x509up_u<999999> /afs/cern.ch/user/s/ssaumya/private/x509up_u<999999>
python3 submit_condor_Run3_FromFiles.py --n 1 --config config/Run3_2025G_NGT_Studies.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999> --farm Jobs --fileList etc/files_to_run.txt
```

### Plotting

### Run Iason's tool for plotting
```
ssh -o ServerAliveInterval=10 ssaumya@lxplus9.cern.ch -L8777:localhost:8777
cd /afs/cern.ch/work/s/ssaumya/private/Egamma/IasonTool/egamma-tnp/
source egmtnpenv/bin/activate
jupyter lab --no-browser --port 8777

# Code template: https://github.com/saumyaphor4252/egamma-tnp/tree/2025_Studies
# Notebook template: STEAM_23ndMay2025.ipynb
```


