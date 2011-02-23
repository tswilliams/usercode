import FWCore.ParameterSet.Config as cms

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

#Reading in a list of datafile locations from 
f_DatafileList = open("/opt/ppd/scratch/williams/Datafiles/Locations/BstdZee/BoostedZeeBkgdLFNs_SMZee_11-02-19.txt")
datafilesList = f_DatafileList.readlines()
f_DatafileList.close()
datafileLocations = map(DataFileLocationAdaptor,datafilesList)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(datafileLocations)
)

#Defining the output file to store the histograms in...
process.TFileService = cms.Service("TFileService", fileName=cms.string('ZandEleKinematicDistns.root'))

process.demo = cms.EDAnalyzer('EvtKinAnr')

process.p = cms.Path(process.demo)
