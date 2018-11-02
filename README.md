cmsrel CMSSW_9_2_0

cd CMSSW_9_2_0/src/

cmsenv

git cms-init


cp -r /nfs/hepwrk01/work/rish/L1TkNtupler/MLTraining/CMSSW_9_2_0/src/* .

git clone https://gitlab.cern.ch/L1TrackJetsPrimaryVtx/TrackletL1PVJets.git

rm -rf ./L1Trigger/TrackletL1PVJets

mv TrackletL1PVJets ./L1Trigger



