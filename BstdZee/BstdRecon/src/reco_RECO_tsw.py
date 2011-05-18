# Auto generated configuration file
# using: 
# Revision: 1.284.2.2 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s RECO --conditions GR_R_311_V2::All --eventcontent RECO
import FWCore.ParameterSet.Config as cms

isMC = True

process = cms.Process('BOOSTEDRECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/data/Run2011A/Photon/RECO/PromptReco-v1/000/161/312/FEEC8AA9-EC57-E011-9C88-0030487A3232.root')
    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltTest/hltTest-GRunV21_Sig-100evts_2011-05-06_secondaryInputFile.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams//Reconstruction/CMSSW_4_1_2/src/signalmcFile_2-00TeVu_GEN-SIM-RECO_test11Apr18a_100evts.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.284.2.2 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition
myCustomOutputCommands = cms.untracked.vstring(
    'keep recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*',
    'keep *_simpleSuperClusters_*_*',
    'keep *_ecalDrivenGsfElectronCores_*_*',
    'keep *_ecalDrivenGsfElectrons_*_*'
)
if(isMC):
   process.RECOSIMEventContent.outputCommands.extend(myCustomOutputCommands)
   myOutputCommands = process.RECOSIMEventContent.outputCommands
else:
   process.RECOEventContent.outputCommands.extend(myCustomOutputCommands)
   myOutputCommands = process.RECOEventContent.outputCommands

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = myOutputCommands,
#    fileName = cms.untracked.string('reco_RECO_tsw_test.root'),
#    fileName = cms.untracked.string('dataFile_bothRECO_spec412ReconV2.root'),
    fileName = cms.untracked.string('mcFile_GEN-SIM-bothRECO_spec412ReconV2.root'),
#    fileName = cms.untracked.string('localWflowTest_bothRECO_Sig-100evts_2011-05-06.root'),
#    fileName = cms.untracked.string('signalmcFile_2-00TeVu_GEN-SIM-bothRECO_test11May6a_100evts.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

# Define the global tag
if(isMC):
   process.GlobalTag.globaltag = 'MC_311_V2::All'
else:
   process.GlobalTag.globaltag = 'GR_R_311_V2::All'


# Modified electron reco
# Product flow is:
#  EB: RecHits -> Multi5x5BasicCluster
#      -> simpleSuperClusters -> Electron Reco
#  EE: RecHits -> Multi5x5BasicCluster -> simpleSuperClusters
#      -> multi5x5SuperClustersWithPreshower -> Electron Reco
process.load('JacksonJ.BoostedReco.simpleSuperClusters_cff')
process.multi5x5BasicClusters.doEndcap = True
process.multi5x5BasicClusters.doBarrel = True
process.multi5x5SuperClustersWithPreshower.endcapSClusterProducer = cms.InputTag("simpleSuperClusters", "endcap")
process.multi5x5PreshowerClusterShape.endcapSClusterProducer = cms.InputTag("simpleSuperClusters", "endcap")
process.ecalDrivenElectronSeeds.barrelSuperClusters = cms.InputTag("simpleSuperClusters", "barrel")
process.ecalDrivenElectronSeeds.endcapSuperClusters = cms.InputTag("multi5x5SuperClustersWithPreshower")

# The whole new re-reco chain
process.reconstruction_step_tmp = cms.Path(
    # New clustering
    process.multi5x5BasicClusters *
    process.simpleSuperClusters *
    process.multi5x5SuperClustersWithPreshower *
    process.multi5x5PreshowerClusterShape *
    # Required transient tracker products
    process.trackerlocalreco *
    #process.siPixelRecHits *
    process.trackingGlobalReco *
    process.electronGsfTracking *
    # Particle flow
    #process.particleFlowCluster *
    #process.particleFlowReco *
    #process.pfElectronTranslator *
    # Electrons! Phew!
    #process.gsfElectronSequence
    #process.electronSequence
    process.ecalDrivenGsfElectronCores*process.ecalDrivenGsfElectrons
)

# Remove things not required from tracking on reco
modulesToRemove = list()
modulesToRemove.append(process.siPixelClusters)
modulesToRemove.append(process.siStripZeroSuppression)
modulesToRemove.append(process.siStripClusters)
process.reconstruction_step = process.reconstruction_step_tmp.copyAndExclude(modulesToRemove)

# Validation
#process.TFileService = cms.Service("TFileService", 
#    fileName = cms.string("histo.root"),
#    closeFileFast = cms.untracked.bool(True)
#)
#process.load("JacksonJ.Analysis.analysis_cfi")
#process.analysis_step = cms.Path(process.analysis)

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.reconstruction_step,process.analysis_step,process.endjob_step)
#process.schedule = cms.Schedule(process.reconstruction_step,process.analysis_step,process.endjob_step,process.RECOoutput_step)
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.RECOoutput_step)

