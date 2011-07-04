import FWCore.ParameterSet.Config as cms

#Useful function for local running...

def DataFileLocationAdaptor(fileLocation):
   fileLocation.replace("\n","")
   filePrefix="file:"
   if(fileLocation.startswith("/pnfs/")==True):
      filePrefix="dcap://heplnx209.pp.rl.ac.uk:22125"
   elif(fileLocation.startswith("/store/")==True):
      filePrefix=""
   return filePrefix + fileLocation

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )


##### Signal MC file lists:
#fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/Boosted_Quark_To_ee_750GeV/tsw-bstd412Recon_0-75TeVu_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"
#fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/Boosted_Quark_To_ee_1TeV/tsw-bstd412Recon_1-00TeVu_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"
fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/Boosted_Quark_To_ee_2TeV/tsw-bstd412Recon_2-00TeVu_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"

##### Bkgd MC file lists:
#fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/tsw-bstd412Recon_Fall10-DYJetsToLL_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"

#Reading in a list of datafile locations from a file...
f_DatafileList = open(fileName_DatafileList)
datafilesList = f_DatafileList.readlines()
f_DatafileList.close()
datafileLocations = map(DataFileLocationAdaptor,datafilesList)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(datafileLocations)
    #fileNames = cms.untracked.vstring('/store/user/tsw/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/bstd412Recon_Fall10-DYJetsToLL_v2/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_8_1_7jk.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_2TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_1TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_0-75TeVu_2011-05-05.root')

)

#Defining the output file to store the histograms/NTuples in...
process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/sigMC_41X-Ntuple_2-00TeVu_20kEvts_2011-06-23.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/bkgdMC_41X-Ntuple_DYJetsToLL_2530516Evts_2011-06-13.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('testNTuple_2-00TeVu.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('localWflowTest_Ntuple_Sig-100evts_2011-05-06.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('test.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/scratch/williams/Data/samplemcFile_Fall10_NTuple_2011-Apr-25.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('sampleSignalmcFile_NTuple-100evts_2011-Apr-25.root'))

process.demo = cms.EDAnalyzer('BstdZeeNTupler',
                              dyJetsToLL_EventType = cms.untracked.int32(0), #==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
                              isMC = cms.untracked.bool(True),
                              printOutInfo = cms.untracked.bool(False),
                              readInNormReco = cms.untracked.bool(True),
                              readInBstdReco = cms.untracked.bool(True),
                              is2010SignalDataset = cms.untracked.bool(True)
      )


process.p = cms.Path(process.demo)
