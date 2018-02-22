# define basic process
import FWCore.ParameterSet.Config as cms
import os
import sys
process = cms.Process("L1Tracklet")
 

# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff') ## this needs to match the geometry you are running on
process.load('Configuration.Geometry.GeometryExtended2023D13_cff')     ## this needs to match the geometry you are running on

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
Source_Files = cms.untracked.vstring(
	'file:/fdata/hepx/store/user/rish/CMSSW_9_2_0/src/L1Trigger/TrackFindingTracklet/test/TTBarPU200/%s' %sys.argv[2]
    )
process.source = cms.Source("PoolSource", fileNames = Source_Files)


# remake stubs 
# ===> IMPORTANT !!! stub window tuning as is by default in CMSSW is incorrect !!! <===
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)


# L1 tracking
process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")
process.TTTracks = cms.Path(process.L1TrackletTracks)                         #run only the tracking (no MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators) #run the tracking AND MC truth associators)


# output module
process.out = cms.OutputModule( "PoolOutputModule",
                                fileName = cms.untracked.string("Tracklets_%s.root" %sys.argv[3]),
                                fastCloning = cms.untracked.bool( False ),
                                outputCommands = cms.untracked.vstring('drop *',
                                                                       'keep *_TTTrack*_Level1TTTracks_*', 
                                                                       'keep *_*AssociatorFromPixelDigis*_*_RERUNL1',
                                                                       #'keep *_*AssociatorFromPixelDigis*_RERUNL1_*',
                                                                       'keep *_TTStubsFromPhase2TrackerDigis*_*_*',
                                                                       'keep *_gen*_*_*',
                                                                       'keep *_ak4GenJets*_*_*',
                                                                       'keep *_*_MergedTrackTruth_*'
)
)
process.FEVToutput_step = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.TTClusterStub,process.TTTracksWithTruth,process.FEVToutput_step)

