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
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )


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
      fileNames = cms.untracked.vstring('/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_9_1_8Xq.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_8_1_WAc.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_7_1_Ysu.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_6_1_Cah.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_5_1_ym3.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_50_1_mVB.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_4_1_HEY.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_49_1_kmp.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_48_1_uBg.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_47_1_45m.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_46_1_gHK.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_45_1_Kg6.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_44_1_g0B.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_43_1_HsO.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_42_1_f6d.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_41_1_BsU.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_40_1_n6y.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_3_1_WSd.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_39_1_hGa.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_38_1_igl.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_37_1_ueX.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_36_1_pqo.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_35_1_Oov.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_34_1_Duw.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_33_1_SaQ.root',
       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/std42XRecon-HLTwL1rerun_2-00TeVu_v1b/8c01a0f969578bad7359cae7a66aa176/mcFile_GEN-SIM-RECO-HLT_32_1_7Tc.root')
##########
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_2TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_1TeVu_2011-05-05.root')
#    fileNames = cms.untracked.vstring('file:/opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_4/src/hltResults/2011-05-05/hltCheck-GRunV21_0-75TeVu_2011-05-05.root')

)

#Defining the output file to store the histograms/NTuples in...
process.TFileService = cms.Service("TFileService", fileName=cms.string('testNTuple_2-00TeVu.root'))
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
                              dyJetsToLL_EventType = cms.untracked.int32(0), #==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
                              isMC = cms.untracked.bool(True),
                              printOutInfo = cms.untracked.bool(False),
                              readInNormReco = cms.untracked.bool(True),
                              readInBstdReco = cms.untracked.bool(False),
                              is2010SignalDataset = cms.untracked.bool(True)
      )


#process.p = cms.Path(process.refitMuons * process.demo)
process.p = cms.Path(process.demo)
