############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
import sys
process = cms.Process("L1TrackJets")

GEOMETRY = "D13" 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')



from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

Source_Files = cms.untracked.vstring(             
#"file:Tracklets_0.root"
"file:/nfs/hepwrk01/work/rish/PrimaryVtx_%s.root" %sys.argv[2]  
#'file:/fdata/hepx/store/user/rish/CMSSW_9_2_0/src/L1Trigger/TrackFindingTracklet/test/TTBarPU200/%s' %sys.argv[2],
)
process.source = cms.Source("PoolSource", fileNames = Source_Files)
#process.load("L1Trigger.TrackFindingTracklet.L1TrackJetsFast_cfi")
#process.TFileService = cms.Service("TFileService", fileName = cms.string('TTBar_'+GEOMETRY+'_NoPU_%d.root' %(int(sys.argv[3]))), closeFileFast = cms.untracked.bool(True))
process.L1TrackJetsFast= cms.EDAnalyzer('L1TrackJetFastProducer',
#
# Default parameters used for the plots for the TP
#
     #HepMCInputTag = cms.InputTag("generator"),

     GenParticleInputTag = cms.InputTag("genParticles",""),
     TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
     L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks","L1Tracklet"),
     MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks","L1Tracklet"), ## MCTruth input
     MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted","L1Tracklet"),
     GenJetAK4=cms.InputTag("ak4GenJetsNoNu"),
     RecoVertexInputTag=cms.InputTag("L1TechPropPrimaryVertexProducer", "l1vertices"),
     MCVertexInputTag=cms.InputTag("L1TkPrimaryVertexMC", ""),
     TrueVertexInputTag=cms.InputTag("L1TechPropPrimaryVertexProducer","TrackingParticleVtx"),
     MyProcess = cms.int32(1),
     DebugMode = cms.bool(False),      # printout lots of debug statements
     SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
     L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
     ZMAX = cms.double ( 15. ) ,        # in cm
     CHI2MAX = cms.double( 5. ),
     PTMINTRA = cms.double( 2.),        # PTMIN of L1Tracks, in GeV
     TP_minNStub = cms.int32( 5) ,       # minimum number of stubs
     nStubsPSmin = cms.int32( -1 ),       # minimum number of stubs in PS modules
     TP_minNStubLayer=cms.int32(5),
     PTMAX = cms.double( 200. ),          # in GeV. When PTMAX > 0, tracks with PT above PTMAX are considered as
                                         # mismeasured and are treated according to HighPtTracks below.
                                         # When PTMAX < 0, no special treatment is done for high PT tracks.
                                         # If PTMAX < 0, no saturation or truncation is done.
     HighPtTracks = cms.int32( 1 ),      # when = 0 : truncation. Tracks with PT above PTMAX are ignored
                                         # when = 1 : saturation. Tracks with PT above PTMAX are set to PT=PTMAX.
     doTightChi2 = cms.bool( True ),    # chi2dof < 5 for tracks with PT > 10
     DeltaZ0Cut=cms.double(0.5),
     CONESize=cms.double(0.4)
)
process.L1TkJets = cms.Path(process.L1TrackJetsFast)

process.TFileService = cms.Service("TFileService", fileName = cms.string('TTBar_'+GEOMETRY+'_PU200_%d.root' %(int(sys.argv[2]))), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('MinBias_'+GEOMETRY+'_PU200_%d.root' %(int(sys.argv[3]))), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('QuarkGun_'+GEOMETRY+'_PU200_%d.root' %(int(sys.argv[3]))), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('BQuarkGun_'+GEOMETRY+'_PU200_%d.root' %(int(sys.argv[3]))), closeFileFast = cms.untracked.bool(True))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('QCDPtHat30_'+GEOMETRY+'_NoPU_%d.root' %(int(sys.argv[3]))), closeFileFast = cms.untracked.bool(True))
#process.schedule = cms.Schedule(process.TTTracksWithTruth,process.ana)
#process.schedule = cms.Schedule(process.TTClusterStub,process.TTTracksWithTruth,process.ana)
#process.schedule(process.L1TkJets)
