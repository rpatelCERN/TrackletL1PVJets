import FWCore.ParameterSet.Config as cms

L1TrackJetsFast= cms.EDAnalyzer('L1TrackJetFastProducer',
#
# Default parameters used for the plots for the TP
#
     #HepMCInputTag = cms.InputTag("generator"),

     GenParticleInputTag = cms.InputTag("genParticles",""),
     TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
     L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
     MCTruthTrackInputTag = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks"), ## MCTruth input
     #MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
     GenJetAK4=cms.InputTag("ak4GenJetsNoNu"), 
     RecoVertexInputTag=cms.InputTag("L1TechPropPrimaryVertexProducer", ""),
     MCVertexInputTag=cms.InputTag("L1TkPrimaryVertexMC", ""),
     TrueVertexInputTag=cms.InputTag("L1TechPropPrimaryVertexProducer","TrackingParticleVtx"),
     MyProcess = cms.int32(1),
     DebugMode = cms.bool(False),      # printout lots of debug statements
     SaveAllTracks = cms.bool(True),   # save *all* L1 tracks, not just truth matched to primary particle
     L1Tk_nPar = cms.int32(4),         # use 4 or 5-parameter L1 track fit ??
     ZMAX = cms.double ( 15. ) ,        # in cm
     CHI2MAX = cms.double( 5. ),
     PTMINTRA = cms.double( 2.),        # PTMIN of L1Tracks, in GeV
     TP_minNStub = cms.int32( 4 ) ,       # minimum number of stubs
     nStubsPSmin = cms.int32( -1 ),       # minimum number of stubs in PS modules 
     TP_minNStubLayer=cms.int32(4),
     PTMAX = cms.double( 200. ),          # in GeV. When PTMAX > 0, tracks with PT above PTMAX are considered as
					 # mismeasured and are treated according to HighPtTracks below.
					 # When PTMAX < 0, no special treatment is done for high PT tracks.
					 # If PTMAX < 0, no saturation or truncation is done.
     HighPtTracks = cms.int32( 1 ),	 # when = 0 : truncation. Tracks with PT above PTMAX are ignored 
					 # when = 1 : saturation. Tracks with PT above PTMAX are set to PT=PTMAX.
     doTightChi2 = cms.bool( True ),    # chi2dof < 5 for tracks with PT > 10
     DeltaZ0Cut=cms.double(0.5),
     CONESize=cms.double(0.4)
)
