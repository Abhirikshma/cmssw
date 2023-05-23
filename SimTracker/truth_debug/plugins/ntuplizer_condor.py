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

process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")

options = VarParsing.VarParsing('standard')
options.register('sourcefile', '', VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.string, 'Source file')
options.parseArguments()

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            #'file:/afs/cern.ch/user/a/abhiriks/CMSSW_12_6_4/src/1325.0_TTbar_13+TTbar_13+DIGIUP15+RECOUP15+HARVESTUP15+ALCATTUP15+NANOUP15/RECO_TrackingPart_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO.root'
            #'root://xrootd-cms.infn.it//store/relval/CMSSW_12_4_0_pre2/RelValTTbar_14TeV/FEVTDEBUGHLT/123X_upgrade2018_realistic_v2_v2-v1/2580000/20a860f4-2321-449a-bf13-20669bc27fa6.root'
            options.sourcefile
            )
                            )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.TFileService = cms.Service("TFileService", 
            fileName = cms.string("ntuple.root")
    )

process.analyze = cms.EDAnalyzer("Ntuplizer")

process.p = cms.Path(
    process.simHitTPAssocProducer+
    process.tpClusterProducer+
    process.trackAssociatorByHits+
    process.analyze
)