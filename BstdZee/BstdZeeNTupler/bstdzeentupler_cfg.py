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

## Loading Geometry modules for creation of transient geometry information used in determining eta and phi values of recHits ...
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START42_V13::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


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
      fileNames = cms.untracked.vstring('/store/mc/Summer11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S4_START42_V11-v2/0000/00FA42AB-D1A9-E011-B6A2-003048678FD6.root')
##########
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_2TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_1TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_0-75TeVu_2011-05-05.root')

)

#Defining the output file to store the histograms/NTuples in...
process.TFileService = cms.Service("TFileService", fileName=cms.string('mcNTuple-42X_v1e-Summer11.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/newscratch/williams/Datafiles/NTuples/CMSSW_41X/bkgdMC_41X-Ntuple_DYJetsToLL_2530516Evts_2011-06-13.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('testNTuple_2-00TeVu.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('localWflowTest_Ntuple_Sig-100evts_2011-05-06.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('test.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/scratch/williams/Data/samplemcFile_Fall10_NTuple_2011-Apr-25.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('sampleSignalmcFile_NTuple-100evts_2011-Apr-25.root'))


#process.refitMuons = cms.EDProducer('MuonsFromRefitTracksProducer',
#    # The input MuonCollection from which the starting Muon objects
#    # will be taken. The module will only consider globalMuons from
#    # the merged muon collection, i.e. trackerMuons, stand-alone muons
#    # will be filtered out of the merged MuonCollection.
#    src = cms.InputTag('muons'),
#
#    # The particular set of refit tracks to use. Could also be
#    # 'tevMuons:default', 'tevMuons:picky', or 'tevMuons:firstHit' to
#    # use the corresponding refits; 'none' to use this module as just
#    # a filter for globalMuons (as opposed to trackerMuons or
#    # caloMuons); to make Muons out of the cocktail tracks, 'tevMuons'
#    # by itself must be used (also specifying fromCocktail = True
#    # below).
#    tevMuonTracks = cms.string('tevMuons'),
#
#    # Exactly one of the below boolean flags may be True (determines
#    # the refit track picked for each muon).
#
#    # Whether to call muon::tevOptimized as in the above code and use
#    # the result of the cocktail choice.
#    fromCocktail = cms.bool(True),
#
#    # Whether to replace the input muons' kinematics with that of the
#    # tracker-only fit. I.e., the muon's momentum, vertex, and charge are
#    # taken from the track accessed by reco::Muon::innerTrack().
#    fromTrackerTrack = cms.bool(False),
#
#    # Whether to replace the input muons' kinematics with that of the
#    # tracker-only fit. I.e., the muon's momentum, vertex, and charge are
#    # taken from the track accessed by reco::Muon::innerTrack().
#    fromGlobalTrack = cms.bool(False),
#
#    # Whether to apply the TMR cocktail algorithm. For each muon track, we
#    # start with the first-muon-hit fit. If the difference in -ln(chi^2
#    # tail probability) between the first-muon-hit and tracker-only is
#    # greater than the below prescribed cut value, the tracker-only fit
#    # replaces the first-muon-hit. For further details see XXX.
#    fromTMR = cms.bool(False),
#
#    # The cut value used in the TMR cocktail.
#    TMRcut = cms.double(4.0),
#
#    # Whether to use Adam Everett's sigma-switch method, choosing
#    # between the global track and the tracker track.
#    fromSigmaSwitch = cms.bool(False),
#
#    # The number of sigma to switch on in the above method.
#    nSigmaSwitch = cms.double(2),
#    
#    # The pT threshold to switch at in the above method.
#    ptThreshold = cms.double(200),
#)

process.demo = cms.EDAnalyzer('BstdZeeNTupler',
                              dyJetsToLL_EventType = cms.untracked.int32(15), #==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
                              isMC = cms.untracked.bool(True),
                              printOutInfo = cms.untracked.bool(False),
                              readInNormReco = cms.untracked.bool(True),
                              readInBstdReco = cms.untracked.bool(False),
                              is2010SignalDataset = cms.untracked.bool(False),
                              useReducedRecHitsCollns = cms.untracked.bool(True),
                              hltPathA_possNames = cms.untracked.vstring("HLT_DoublePhoton33_v1",
                                                                         "HLT_DoublePhoton33_v2",
                                                                         "HLT_DoublePhoton33_v3",
                                                                         "HLT_DoubleEle33_CaloIdT_v1",
                                                                         "HLT_DoubleEle33_CaloIdT_v2",
                                                                         "HLT_DoubleEle33_CaloIdL_v1",
                                                                         "HLT_DoubleEle33_CaloIdL_v2",
                                                                         "HLT_DoubleEle33_CaloIdL_v3",
                                                                         "HLT_DoubleEle33_CaloIdL_v4",
                                                                         "HLT_DoubleEle33_CaloIdL_v5" # - run 177878 (3e33, v4.0) # (5e33, v1.4)
                                                                         ),
                              trg_emuPath_possNames = cms.untracked.vstring("HLT_Mu15_Photon20_CaloIdL_v1",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v2",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v3",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v4",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v5",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v6",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v7",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v8",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v9",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v10",# - run 177878 (3e33, v4.0)
                                                                            "HLT_Mu15_Photon20_CaloIdL_v11",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v12",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v13" # (5e33, v1.4)
                                                                            )
      )


#process.p = cms.Path(process.refitMuons * process.demo)
process.p = cms.Path(process.demo)
