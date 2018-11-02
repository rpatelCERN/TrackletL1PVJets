############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import os
import sys
process = cms.Process("L1TrackVtx")

GEOMETRY = "D13"

 
############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
print "using geometry " + GEOMETRY + " (tilted)"
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

#D13 (tilted barrel)
Source_Files = cms.untracked.vstring(
    "file:/mnt/hadoop/store/user/rish/TTBarPU200/Tracklets_%s.root" %sys.argv[2]
    #"file:/fdata/hepx/store/user/rish/CMSSW_9_2_0/src/L1Trigger/TrackFindingTracklet/test/TTBarPU200/reprocess_L1T_L1TTTTracklet_step2_TT_PhaseIISpring17D-PU200_pilot_100.root"
)
process.source = cms.Source("PoolSource", fileNames = Source_Files)

process.load("L1Trigger.TrackletL1PVJets.L1TechPropPrimaryVertexProducer_cfi")
process.L1TechPropPrimaryVertexProducer.L1TrackInputTag = cms.InputTag("TTTracksFromTracklet","Level1TTTracks","L1Tracklet")
process.L1TechPropPrimaryVertexProducer.L1Tk_nPar = cms.int32(4)   ## use 4 or 5 parameter track fit
process.pL1TkPrimaryVertex = cms.Path( process.L1TechPropPrimaryVertexProducer )

process.L1TkPrimaryVertexMC = process.L1TechPropPrimaryVertexProducer.clone()
process.L1TkPrimaryVertexMC.MonteCarloVertex = cms.bool( True )
process.pL1TkPrimaryVertexMC = cms.Path( process.L1TkPrimaryVertexMC )

process.out = cms.OutputModule( "PoolOutputModule",
				fastCloning = cms.untracked.bool( False ),
                                fileName = cms.untracked.string("PrimaryVtx_%s.root" %sys.argv[2]))

process.FEVToutput_step = cms.EndPath(process.out)
process.TFileService = cms.Service("TFileService", fileName = cms.string('TTBar_'+GEOMETRY+'_PU200_%d.root' %(int(sys.argv[2]))), closeFileFast = cms.untracked.bool(True))
