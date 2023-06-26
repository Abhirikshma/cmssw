import FWCore.ParameterSet.Config as cms
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

# process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
# process.load("SimGeneral.TrackingAnalysis.Playback_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("SimTracker.VertexAssociation.VertexAssociatorByTracks_cfi")
process.load("RecoTracker.Configuration.RecoTracker_cff")
# process.load("SimGeneral.MixingModule.trackingTruthProducer_cfi")

# process.load("SimTracker.TrackHistory.SecondaryVertexTagInfoProxy_cff")
# process.load("SimTracker.TrackHistory.Playback_cff")
process.load("SimTracker.TrackHistory.VertexClassifier_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:/afs/cern.ch/user/a/abhiriks/CMSSW_12_6_4/src/1315.0_SingleElectronPt10_UP15+SingleElectronPt10_UP15+DIGIUP15+RECOUP15+HARVESTUP15/RECO_TrackingPart_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
            'file:/afs/cern.ch/user/a/abhiriks/CMSSW_12_6_4/src/1325.0_TTbar_13+TTbar_13+DIGIUP15+RECOUP15+HARVESTUP15+ALCATTUP15+NANOUP15/RECO_TrackingPart_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
# 'root://xrootd-cms.infn.it//store/relval/CMSSW_12_4_0_pre2/RelValTTbar_14TeV/FEVTDEBUGHLT/123X_upgrade2018_realistic_v2_v2-v1/2580000/20a860f4-2321-449a-bf13-20669bc27fa6.root'
            )
                            )
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
    # process.vertexClassifier,
    vertToTvAssociator = cms.InputTag("vertexAssociatorByTracksByHits"),
    # bestMatchByMaxValue = cms.untracked.bool(True),
    # trackingTruth = cms.untracked.InputTag('mix','MergedTrackTruth'),
    # vertexAssociator = cms.untracked.InputTag('vertexAssociatorByTracksByHits'),
    # vertexProducer = cms.untracked.InputTag("inclusiveSecondaryVertices"),
    # enableRecoToSim = cms.untracked.bool(True),
    # enableSimToReco = cms.untracked.bool(False),
    # hepMC = cms.untracked.InputTag("generatorSmeared"),
    # longLivedDecayLength = cms.untracked.double(1e-14),
    # vertexClusteringDistance = cms.untracked.double(0.003)
)

process.p = cms.Path(
    # process.playback+
    process.simHitTPAssocProducer+
    process.tpClusterProducer+
    # process.trackingParticles+
    # process.svTagInfoProxy+
    process.trackAssociatorByHits+
    process.trackingParticleRecoTrackAssociationByHits+
    process.vertexAssociatorByTracksByHits+
    # process.vertexAssociatorSequence+
    process.analyze
)
