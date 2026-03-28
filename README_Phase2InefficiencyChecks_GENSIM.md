## Relevant Links:
- See presentation:


## For Distributions
```
git clone git@github.com:tihsu99/EGamma-Workflow.git
cd EGamma-Workflow/

voms-proxy-init --voms cms --valid 100:00

python3 submit_condor_GENSIM.py   --n-jobs 1000   --config config/phase2_ZpToEE_GENSIM.yaml   --proxy /afs/cern.ch/user/s/ssaumya/private/x509up_u<99999>   --farm Jobs_GENSIM --events-per-job 100
condor_submit Jobs_GENSIM/condor_jobs.sub
```

