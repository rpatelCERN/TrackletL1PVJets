cmsrel CMSSW_10_1_7

cd CMSSW_10_1_7/src/

cmsenv

git cms-init


cp -r /nfs/hepwrk01/work/rish/L1TkNtupler/MLTraining/CMSSW_10_1_7/src/* .

git clone -b NtupleVar.v1 https://gitlab.cern.ch/L1TrackJetsPrimaryVtx/TrackletL1PVJets.git

rm -rf ./L1Trigger/TrackletL1PVJets

mv TrackletL1PVJets ./L1Trigger



