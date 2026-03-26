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
```
python3 ../plot_Phase2_EGM_Sigma_PtRanges.py /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/ZprimeToEE_ntuple.root SigmaPtBins

python3 plot_Phase2_EGM_Variables.py \
  /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/ZprimeToEE_ntuple.root \
  /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/QCD_ntuple.root  \
  Comparison "Z'#rightarrow ee" "QCD"
```

### SB Optimization
```
python3 CutOptimization.py --signal /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_ZprimeToEE.root --background /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_QCD.root -t egHLTTree --target-eff 0.75 --pt-bins "30,5000" --var eg_sigma2vv --pt-var eg_et --var-min 1e-9

python3 CutOptimization.py --signal /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_ZprimeToEE.root --background /eos/cms/store/group/phys_egamma/ssaumya/Phase2_Inefficiency/tnpNtupler_QCD.root -t egHLTTree --target-eff 0.75 --pt-bins "30,400 400,1200 1200,5000" --var eg_sigma2vv --pt-var eg_et --var-min 1e-9
```
## For Trigger Efficinecy Checks
```
cd PlotHists/
# Update the input files in Inputs.py
python3 plotHist.py

```
