## Relevant Links:
- Details: https://github.com/cms-sw/cmssw/issues/42850
- Dataset used: 
	- `13_1_0`: `/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAO`
    - `17_0_0_pre3`: `/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD`

## 131X Set-up using Spring23 samples

- CMSSW: `CMSSW_13_1_0`
- Dataset used: `/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAO`
- Number of events: 98392
- Recipe details available in https://github.com/cms-sw/cmssw/issues/42850#issuecomment-1732468495

### CMSSW set-up
```
cmsrel CMSSW_13_1_0; cd CMSSW_13_1_0/src/; cmsenv;
git cms-init 
git cms-checkout-topic -u cms-l1t-offline:Phase2_prototypeSnapshot_3_CMSSW_13_1
git cms-addpkg HLTrigger/Configuration
git cms-addpkg L1Trigger/Phase2L1Taus
git cms-addpkg HLTrigger/HLTfilters
git cms-addpkg DataFormats/HLTReco
git cms-addpkg HLTrigger/HLTcore
git cms-merge-topic SohamBhattacharya:test_Phase2-HLT_with_L1T-layer2_13_1_0
git cms-merge-topic swagata87:hltSingleEleNonIso131X
git cms-merge-topic Sam-Harper:EGHLTCustomisation_1230pre6
scram b -j 10

```
YAML config:
```
cmssw-dir: /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_13_1_0/src/
eos-dir: /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/Upgrade_PR42850/CMSSW_13_1_0/

project:

  ZprimeToEE:
    dataset:
      - /ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring23DIGIRECOMiniAOD-PU200_Trk1GeV_131X_mcRun4_realistic_v5-v1/GEN-SIM-DIGI-RAW-MINIAOD
    GT: auto:phase2_realistic_T21
    process:
      - cmsDriver.py Phase2 -s RAW2DIGI,L1:RUNP2GT --conditions auto:phase2_realistic_T21 --geometry Extended2026D95 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2.addHcalTriggerPrimitives,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT --nThreads 1 --mc --inputCommands='keep *, drop l1tPFJets_*_*_*' --outputCommands='keep *P2GT*_*_*_*, drop l1tPFJets_*_*_*' --python_filename rerunL1_cfg.py --no_exec --fileout file:output_Phase2_L1T.root
      - cmsRun rerunL1_cfg.py
      - cmsDriver.py Phase2 -s HLT:75e33 --processName=HLTX --conditions auto:phase2_realistic_T21 --geometry Extended2026D95 --era Phase2C17I13M9 --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,HLTrigger/Configuration/customizeHLTforEGamma.customiseEGammaMenuDev --eventcontent FEVTDEBUGHLT --filein file:output_Phase2_L1T.root --inputCommands='keep *, drop *_hlt*_*_HLT' --nThreads 1 --no_exec --python_filename rerunHLT_cfg.py --fileout file:output_Phase2_HLT.root
      - cmsRun rerunHLT_cfg.py
      - python3 /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_13_1_0/src/EGamma-Workflow/python/lastFilterEfficiency_Upgrade.py --input-file output_Phase2_HLT.root -n 100000 -o tnpNtupler.root

    storefile:
      - tnpNtupler.root
```
Submit jobs on condor:
```
python3 submit_condor_MC.py --n 1 --config config/phase2_ZpToEE_TriggerEfficiency_GitIssue-42850_131X.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<999999> --farm Jobs --nEvent 100
```

## 170X set-up using Spring24 samples
- CMSSW: `CMSSW_17_0_0_pre3`
- Dataset used: `/ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD`
- Number of events: 100000
- Recipe details based on https://cmshltupgrade.docs.cern.ch/RunningInstructions/

### CMSSW set-up
```
cmsrel CMSSW_17_0_0_pre3;  cd CMSSW_17_0_0_pre3/src; cmsenv;
git cms-init
git cms-merge-topic Sam-Harper:EGHLTCustomisation_1230pre6
git clone git@github.com:saumyaphor4252/EGamma-Workflow.git
scram b -j 10
```
YAML config:

```
cmssw-dir: /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_17_0_0_pre3/src/
eos-dir: /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/Upgrade_PR42850/CMSSW_17_0_0_pre3/

project:

  ZprimeToEE:
    dataset:
      - /ZprimeToEE_M-6000_TuneCP5_14TeV-pythia8/Phase2Spring24DIGIRECOMiniAOD-PU200_Trk1GeV_140X_mcRun4_realistic_v4-v1/GEN-SIM-DIGI-RAW-MINIAOD
    GT: auto:phase2_realistic_T21
    process:
      - cmsDriver.py Phase2 -s L1,L1TrackTrigger --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D110 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW-MINIAOD --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,Configuration/DataProcessing/Utils.addMonitoring,L1Trigger/Configuration/customisePhase2FEVTDEBUGHLT.customisePhase2FEVTDEBUGHLT,L1Trigger/Configuration/customisePhase2TTOn110.customisePhase2TTOn110 --python_filename rerunL1_cfg.py --inputCommands="keep *, drop l1tPFJets_*_*_*, drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT" --outputCommands="drop l1tTrackerMuons_l1tTkMuonsGmt*_*_HLT" --mc --no_exec --fileout file:output_Phase2_L1T.root
      - cmsRun rerunL1_cfg.py
      - cmsDriver.py Phase2 -s L1P2GT,HLT:75e33 --processName=HLTX --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D110 --era Phase2C17I13M9 --eventcontent FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,HLTrigger/Configuration/customizeHLTforEGamma.customiseEGammaMenuDev --filein file:output_Phase2_L1T.root --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' --mc --no_exec --python_filename rerunHLT_cfg.py --fileout file:output_Phase2_HLT.root
      - cmsRun rerunHLT_cfg.py
      - python3 /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_17_0_0_pre3/src/EGamma-Workflow/python/lastFilterEfficiency_Upgrade.py --input-file output_Phase2_HLT.root -n 100000 -o tnpNtupler.root

    storefile:
      - tnpNtupler.root
```
Submit jobs on condor:
```
python3 submit_condor_MC.py --n 1 --config config/phase2_ZpToEE_TriggerEfficiency_GitIssue-42850_170X.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<999999> --farm Jobs_Spring24 --nEvent 100
```

## 20_0_X set-up using RelVals (Very low stats to provide useful info)
- CMSSW: `CMSSW_20_0_0`
- Dataset used: `/RelValZpToEE_m6000_14TeV/CMSSW_20_0_0-PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v6/GEN-SIM-DIGI-RAW`
- Number of events: 9000
- Recipe details based on: 

### CMSSW set-up
```
cmsrel CMSSW_20_0_0; cd CMSSW_20_0_0/src/; cmsenv;
git cms-init
git cms-addpkg HLT/Configuration
git cms-merge-topic Sam-Harper:EGHLTCustomisation_1230pre6
git clone git@github.com:saumyaphor4252/EGamma-Workflow.git
scram b -j 10
```
YAML config:
```
cmssw-dir: /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_20_0_0/src/
eos-dir: /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/Upgrade_PR42850/CMSSW_20_0_0/

project:

  ZprimeToEE:
    dataset:
      - /RelValZpToEE_m6000_14TeV/CMSSW_20_0_0-PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v6/GEN-SIM-DIGI-RAW
    GT: auto:phase2_realistic_T35
    process:
      - cmsDriver.py Phase2 -s L1P2GT,HLT:75e33 --processName=HLTX --conditions auto:phase2_realistic_T35 --geometry ExtendedRun4D121 --era Phase2C22I13M9 --eventcontent FEVTDEBUGHLT --customise SLHCUpgradeSimulations/Configuration/aging.customise_aging_1000,HLTrigger/Configuration/customizeHLTforEGamma.customiseEGammaMenuDev --inputCommands='keep *, drop *_hlt*_*_HLT, drop triggerTriggerFilterObjectWithRefs_l1t*_*_HLT' --mc --no_exec --python_filename rerun_L1P2GT_HLT_cfg.py --fileout file:output_Phase2_L1P2GT_HLT.root
      - cmsRun rerun_L1P2GT_HLT_cfg.py
      - python3 /afs/cern.ch/work/s/ssaumya/private/Egamma/Upgrade_PR42850/TestRecipe/CMSSW_20_0_0/src/EGamma-Workflow/python/lastFilterEfficiency_Upgrade.py --input-file output_Phase2_L1P2GT_HLT.root -n 100000 -o tnpNtupler.root

    storefile:
      - tnpNtupler.root
```
Submit jobs on condor:
```
python3 submit_condor_MC.py --n 1 --config config/phase2_ZpToEE_TriggerEfficiency_GitIssue-42850_200X.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<999999> --farm Jobs --nEvent 50
```

### Plotting Scripts in plot_func
Update the tuples location in `Inputs.py` and run plotHist.py. 
```
python3 plotHist.py
```
