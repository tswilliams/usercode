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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Reading in a list of datafile locations from a file...
f_DatafileList = open("/opt/ppd/scratch/williams/Datafiles/Locations/BstdZee/BstdZeeLFNs_test412stdRECO_11-04-19.txt")
datafilesList = f_DatafileList.readlines()
f_DatafileList.close()
datafileLocations = map(DataFileLocationAdaptor,datafilesList)

process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(datafileLocations)
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/store/user/tsw/Boosted_Quark_To_ee_1TeV/stdRecon_1-00TeVu_v2/98cc305c25969310e4843d40e3559cee/signalmcFile_GEN-SIM-RECO_v2_1_1_QCb.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/opt/ppd/scratch/williams/SpecRecon/CMSSW_4_1_2/src/sampleSignalmcFile_GEN-SIM-bothRECO-100evts_2011-Apr-25.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/opt/ppd/scratch/williams/SpecRecon/CMSSW_4_1_2/src/run11A-photon-DataFile_bothRECO_test11Apr26.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/opt/ppd/scratch/williams/Data/samplemcFile_Fall10_GEN-SIM-bothRECO_2011-Apr-25.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/opt/ppd/scratch/williams/SpecReconV2/CMSSW_4_1_2/src/reco_RECO_tsw_test.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/opt/ppd/scratch/williams/SpecReconV2/CMSSW_4_1_2/src/localWflowTest_bothRECO_Sig-100evts_2011-05-06.root"))
#    fileNames = cms.untracked.vstring(DataFileLocationAdaptor("/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30e/3054051540f9e6b491340ce955ddcf3c/run11A-photon-DataFile_bothRECO_11Apr30e_1_1_mFp.root"))
    fileNames = cms.untracked.vstring("/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_1_1_Kt0.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_2_1_v8q.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_3_1_C3Y.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_4_1_mP9.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_5_1_wZh.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_6_1_Ulh.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_7_1_xrB.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_8_1_5Qh.root",
                                      "/store/user/tsw/Photon/run11A-photon_bothRECO_11Apr30f/3054051540f9e6b491340ce955ddcf3c/dataFile_bothRECO_spec412ReconV2_9_1_P7u.root")
#
#    fileNames = cms.untracked.vstring('/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_9_1_7Mo.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_8_1_Q00.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_7_1_aD5.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_6_1_eXo.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_5_1_hgG.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_4_1_QeY.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_3_1_EFF.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_2_1_Um2.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_1_1_zZu.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_750GeV/bstd412Recon_0-75TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_10_1_qo4.root' )
#
#    fileNames = cms.untracked.vstring('/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_9_1_6l3.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_8_1_yDF.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_7_1_BK6.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_6_1_OyC.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_5_1_xV8.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_4_1_4cF.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_3_1_jxx.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_2_1_ZQO.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_1_1_DH6.root',
#                                      '/store/user/tsw/Boosted_Quark_To_ee_1TeV/bstd412Recon_1-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_10_1_FEw.root' )
#
#    fileNames = cms.untracked.vstring(
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_9_1_OeI.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_8_1_La9.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_7_1_Smq.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_6_1_xGM.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_5_1_uh4.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_4_1_3e4.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_3_1_Mvv.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_2_1_qGY.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_1_1_1Bo.root',
#       '/store/user/tsw/Boosted_Quark_To_ee_2TeV/bstd412Recon_2-00TeVu_v1/c17b3790b1be6b1dd249df112816fe9d/mcFile_GEN-SIM-bothRECO_spec412ReconV2_10_1_QIB.root' )
)

#Defining the output file to store the histograms/NTuples in...
#process.TFileService = cms.Service("TFileService", fileName=cms.string('signalmcNTuple_2-00TeVu_20kEvts_2011-05-13.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('signalmcNTuple_0-75TeVu_20kEvts_2011-05-13.root'))
process.TFileService = cms.Service("TFileService", fileName=cms.string('dataNTuple_170kEvts_2011-05-13.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('localWflowTest_Ntuple_Sig-100evts_2011-05-06.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('test.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('/opt/ppd/scratch/williams/Data/samplemcFile_Fall10_NTuple_2011-Apr-25.root'))
#process.TFileService = cms.Service("TFileService", fileName=cms.string('sampleSignalmcFile_NTuple-100evts_2011-Apr-25.root'))

process.demo = cms.EDAnalyzer('BstdZeeNTupler',
                              isMC = cms.untracked.bool(False),
                              printOutInfo = cms.untracked.bool(False)
      )


process.p = cms.Path(process.demo)
