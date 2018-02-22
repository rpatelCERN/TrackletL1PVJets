cmsrel CMSSW_9_2_0

cd CMSSW_9_2_0/src/

cmsenv

git cms-init

git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git

git fetch cms-l1t-offline

git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v1.14.1


