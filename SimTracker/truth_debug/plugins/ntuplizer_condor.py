import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process("Ntupling", Run2_2016)

# Conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("SimTracker.VertexAssociation.VertexAssociatorByTracks_cfi")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("SimTracker.TrackHistory.VertexClassifier_cff")

options = VarParsing.VarParsing('standard')
options.register('sourcefile', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, 'Source file')
options.parseArguments()

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(options.sourcefile)
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", 
            fileName = cms.string("ntuple.root")
    )

process.trackingParticleRecoTrackAssociationByHits = process.trackingParticleRecoTrackAsssociation.clone(
    associator = 'trackAssociatorByHits'
)
process.vertexAssociatorByTracksByHits = process.VertexAssociatorByTracks.clone(
    trackAssociation = "trackingParticleRecoTrackAssociationByHits",
)

process.analyze = cms.EDAnalyzer("Ntuplizer",
    vertToTvAssociator = cms.InputTag("vertexAssociatorByTracksByHits"),
    enableDebug = cms.bool(False)
)

process.p = cms.Path(
    process.simHitTPAssocProducer+
    process.tpClusterProducer+
    process.trackAssociatorByHits+
    process.trackingParticleRecoTrackAssociationByHits+
    process.vertexAssociatorByTracksByHits+
    process.analyze
)

