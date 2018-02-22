############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1TrackNtuple")

GEOMETRY = "D13"

 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

if GEOMETRY == "D10": 
    print "using geometry " + GEOMETRY + " (flat)"
    process.load('Configuration.Geometry.GeometryExtended2023D10Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D10_cff')
elif GEOMETRY == "D13":
    print "using geometry " + GEOMETRY + " (tilted)"
    process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')
    process.load('Configuration.Geometry.GeometryExtended2023D13_cff')
else:
    print "this is not a valid geometry!!!"

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')


############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

if GEOMETRY == "D10": 
    #D10 (flat barrel)
    Source_Files = cms.untracked.vstring(
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/044925F3-5F2E-E711-A92B-0CC47A7AB7A0.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/42543260-602E-E711-A41C-0025905A48C0.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/72D6B8B7-612E-E711-A727-0CC47A4D7692.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/786DD1A7-612E-E711-84EB-0025905B855E.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/9255C706-602E-E711-BDE0-0025905A609A.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/BE75B044-602E-E711-8DC3-0CC47A4D7692.root",
    "/store/relval/CMSSW_9_1_0_pre3/RelValSingleMuPt10Extended/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D10-v1/10000/EA054BA3-5F2E-E711-B7F4-0025905A6084.root",
)
elif GEOMETRY == "D13":
    #D13 (tilted barrel)
    Source_Files = cms.untracked.vstring(
	"file:Tracklets_0.root"
	
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/0E9A4F45-602E-E711-916C-0CC47A7C3430.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/12249F9F-5F2E-E711-849D-0CC47A4D76C6.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/4A2589F6-622E-E711-8EB2-0025905B85CC.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/4C511711-5F2E-E711-847A-0025905B857E.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/52EA294A-602E-E711-B368-0025905B85CA.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/5ECEC385-5E2E-E711-8E82-0CC47A7C3430.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/7A367D94-5D2E-E711-A3F2-0025905B8568.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/BAABF1B9-5E2E-E711-BE64-0025905A60E4.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/EC59F295-602E-E711-BAB5-0025905A48C0.root",
    #"/store/relval/CMSSW_9_1_0_pre3/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/91X_upgrade2023_realistic_v1_D13-v2/10000/ECC7C5F2-622E-E711-8704-0025905AA9F0.root",
)
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.TFileService = cms.Service("TFileService", fileName = cms.string('Muon10_'+GEOMETRY+'_PU0.root'), closeFileFast = cms.untracked.bool(True))


############################################################
# L1 tracking
############################################################

# remake stubs 

# ===> IMPORTANT !!! stub window tuning as is by default in CMSSW is incorrect !!! <===

process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *

if GEOMETRY == "D10": 
    TTStubAlgorithm_official_Phase2TrackerDigi_.zMatchingPS = cms.bool(False)
process.TTClusterStub = cms.Path(process.TrackTriggerClustersStubs)

process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")

from L1Trigger.TrackFindingTracklet.Tracklet_cfi import *
if GEOMETRY == "D10": 
    TTTracksFromTracklet.trackerGeometry = cms.untracked.string("flat")
#TTTracksFromTracklet.asciiFileName = cms.untracked.string("evlist.txt")

# run only the tracking (no MC truth associators)
process.TTTracks = cms.Path(process.L1TrackletTracks)

# run the tracking AND MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


############################################################
# Define the track ntuple process, MyProcess is the (unsigned) PDGID corresponding to the process which is run
# e.g. single electron/positron = 11
#      single pion+/pion- = 211
#      single muon+/muon- = 13 
#      pions in jets = 6
#      taus = 15
#      all TPs = 1
############################################################

process.L1TrackNtuple = cms.EDAnalyzer('L1TrackNtupleMaker',
                                       MyProcess = cms.int32(1),
                                       DebugMode = cms.bool(False),      # printout lots of debug statements
                                       SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
                                       SaveStubs = cms.bool(False),      # save some info for *all* stubs
                                       L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
                                       L1Tk_minNStub = cms.int32(4),     # L1 tracks with >= 4 stubs
                                       TP_minNStub = cms.int32(4),       # require TP to have >= X number of stubs associated with it
                                       TP_minNStubLayer = cms.int32(4),  # require TP to have stubs in >= X layers/disks
                                       TP_minPt = cms.double(2.0),       # only save TPs with pt > X GeV
                                       TP_maxEta = cms.double(2.4),      # only save TPs with |eta| < X
                                       TP_maxZ0 = cms.double(30.0),      # only save TPs with |z0| < X cm
                                       L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),               ## TTTrack input
                                       MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input 
                                       # other input collections
                                       L1StubInputTag = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                       MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                       MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                       TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                       )
process.ana = cms.Path(process.L1TrackNtuple)

process.schedule = cms.Schedule(process.TTClusterStub,process.TTTracksWithTruth,process.ana)

