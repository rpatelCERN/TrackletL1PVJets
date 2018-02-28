import FWCore.ParameterSet.Config as cms

L1TechPropPrimaryVertexProducer = cms.EDProducer('L1TechPropPrimaryVertexProducer',
#
# Default parameters used for the plots for the TP
#
     #HepMCInputTag = cms.InputTag("generator"),
     GenParticleInputTag = cms.InputTag("genParticles",""),
     TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
     L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks","L1Tracklet"),
     MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted","L1Tracklet"),
     ZMAX = cms.double ( 15. ) ,        # in cm
     CHI2MAX = cms.double( 5. ),
     PTMINTRA = cms.double( 2.),        # PTMIN of L1Tracks, in GeV
     nStubsmin = cms.int32( 5 ) ,       # minimum number of stubs
     nStubsPSmin = cms.int32( -1 ),       # minimum number of stubs in PS modules 
     #nBinning = cms.int32( 301 ),        # number of bins for the temp histo (from -30 cm to + 30 cm) THis is replaced by DeltaZ
     PTMAX = cms.double( 200. ),          # in GeV. When PTMAX > 0, tracks with PT above PTMAX are considered as
					 # mismeasured and are treated according to HighPtTracks below.
					 # When PTMAX < 0, no special treatment is done for high PT tracks.
					 # If PTMAX < 0, no saturation or truncation is done.
     HighPtTracks = cms.int32( 1 ),	 # when = 0 : truncation. Tracks with PT above PTMAX are ignored 
					 # when = 1 : saturation. Tracks with PT above PTMAX are set to PT=PTMAX.
     MonteCarloVertex = cms.bool( False ),    #  when True: dont run the vxt finding algo but pick up the MC generated vtx
     doPtComp = cms.bool( False ),       # track-stubs PT compatibility cut
     doTightChi2 = cms.bool( True ),    # chi2dof < 5 for tracks with PT > 10
     WEIGHT = cms.int32(1),            # WEIGHT can be set to 0, 1 or 2 for unweighted, pT weighted
                                      # or pT2 weighted tracks respectively.

     SumPtSquared=cms.bool(False),
     DeltaZ=cms.double(0.5)
#
# Other working point which works better for H -> TauTau,
# cf talk by Moshan Ather, Dec 12, 2014:

#     WEIGHT = cms.int32(2),
#     PTMAX = cms.double( 25. ),
#     nStubsmin = cms.int32( 5 ),
#     HighPtTracks = cms.int32( 1),
#     doPtComp = cms.bool( False ),     
#     CHI2MAX = cms.double( 20 )
#

)
