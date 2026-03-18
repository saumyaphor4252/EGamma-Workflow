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

```

## For Trigger Efficinecy Checks
```



```
