import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2016_cff import Run2_2016

process = cms.Process("USER", Run2_2016)

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

process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")
# process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
# process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

# process.load("SimGeneral.MixingModule.mixNoPU_cfi")
# process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")
# process.trackingParticles.simHitCollections = cms.PSet( )
# process.mix.playback = cms.untracked.bool(True)
# process.mix.digitizers = cms.PSet(
#      mergedtruth = cms.PSet(process.trackingParticles)
# )
# for a in process.aliases: delattr(process, a)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:/afs/cern.ch/user/a/abhiriks/CMSSW_12_6_4/src/1315.0_SingleElectronPt10_UP15+SingleElectronPt10_UP15+DIGIUP15+RECOUP15+HARVESTUP15/RECO_TrackingPart_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
            'file:/afs/cern.ch/user/a/abhiriks/CMSSW_12_6_4/src/1325.0_TTbar_13+TTbar_13+DIGIUP15+RECOUP15+HARVESTUP15+ALCATTUP15+NANOUP15/RECO_TrackingPart_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
            )
                            )

process.analyze = cms.EDAnalyzer("truthDebug")

process.p = cms.Path(
    process.simHitTPAssocProducer+
    process.tpClusterProducer+
    process.trackAssociatorByHits+
    # process.trackingParticleRecoTrackAsssociation+
    process.analyze
)