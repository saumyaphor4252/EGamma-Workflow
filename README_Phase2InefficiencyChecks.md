## Relevant Links:
- See presentation: https://indico.cern.ch/event/1638160/contributions/6891732/subcontributions/595453/attachments/3243858/5786989/EGM_HLT_Upgrade_Inefficiency_Checks_24.03.2026.pdf


## For Distributions
```
git clone git@github.com:tihsu99/EGamma-Workflow.git
cd EGamma-Workflow/

voms-proxy-init --voms cms --valid 100:00

python3 submit_condor_MC.py --n 1 --config config/phase2_ZpToEE.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999> --farm Jobs --nEvent 100
python3 submit_condor_MC.py --n 1 --config config/phase2_QCD_Pt-300ToInf.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999> --farm Jobs --nEvent 100
python3 submit_condor_MC.py --n 1 --config config/phase2_QCD_Pt-15To3000.yaml --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999> --farm Jobs --nEvent 100
condor_submit Farm/condor_jobs.sub


```

### Plotting
- Plotting the sigmaVV and sigmaWW distributions in different pT ranges: 30-400, 400-1000, and >1000GeV
```
python3 PlottingScripts/plot_Phase2_EGM_Sigma_PtRanges.py /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/ZprimeToEE_ntuple.root SigmaPtBins
```

- Plotting the overlapping distributions of all the variables (QCD Vs ZpToEE) (Pg3-4)
```
python3 PlottingScripts/plot_Phase2_EGM_Variables.py \
  /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/ZprimeToEE_ntuple.root \
  /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/QCD_ntuple.root  \
  Comparison "Z'#rightarrow ee" "QCD"
```

- Plotting the 2D sigmaVV vs pT and eta distributions (Pg 5)
```
python3 PlottingScripts/plot_Phase2_EGM_Sigma_Colz_PtEta.py /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_QCD.root
python3 PlottingScripts/plot_Phase2_EGM_Sigma_Colz_PtEta.py /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_ZprimeToEE.root
```

- Plotting sigmaVV per Filter ()
```
python3 plotSigma2vvSigma2ww_PerFilter.py \
  --input-file /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/Efficinecy/sigma_FilterDistributions/sigma_per_filter.root \
  --out_dir sigma_per_filter 
```

- Cut Optimzation (See below)

### SB Optimization
```
python3 python/CutOptimization.py --signal /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_ZprimeToEE.root --background /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_QCD.root -t egHLTTree --target-eff 0.75 --pt-bins "30,5000" --var eg_sigma2vv --pt-var eg_et --var-min 1e-9

python3 python/CutOptimization.py --signal /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_ZprimeToEE.root --background /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_QCD.root -t egHLTTree --target-eff 0.75 --pt-bins "30,400 400,1200 1200,5000" --var eg_sigma2vv --pt-var eg_et --var-min 1e-9
```



## For Trigger Efficinecy Checks (Pg 12)
```
cd PlotHists/
# Update the input files in Inputs.py
python3 plotHist.py
```

### Configs Details:
- phase2_ZpToEE.yaml: For current Z'->ee signal distributions
- phase2_QCD_Pt-15To3000.yaml, phase2_QCD_Pt-300ToInf.yaml: For current Background distributions 
- phase2_ZpToEE_TriggerEfficiency.yaml: For Current/Reference Z'->ee signal trigger efficiency
- phase2_ZpToEE_TriggerEfficiency_VersionFinal.yaml: For Target Z'->ee signal trigger efficiency
- phase2_QCD_FakeRate.yaml: For Fake rate calculation using QCD samples
- phase2_ZpToEE_FilterDistributions.yaml: To get filterwise distrbutions of the Z'->ee sample (eg sigmaVV and sigmaWW here)
- phase2_ZpToEE_GENSIM.yaml: To get high stats sample for distributions of the Z'->ee samples. Currently 100k not enough to check the distrbutions for energy > 500GeV 




