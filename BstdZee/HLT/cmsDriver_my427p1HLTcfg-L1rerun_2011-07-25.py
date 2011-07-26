# Auto generated configuration file
# using: 
# Revision: 1.303.2.7 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: RelVal --step=RAW2DIGI,L1,HLT:GRun --conditions=auto:startup --filein=/store/user/tsw/Boosted_Quark_To_ee_2TeV/std424Recon_2-00TeVu_v1/30a3840207bc3be6b1c52e0833382547/sigMCFile_GEN-SIM-std424RECO_v1_2011-07-05_9_1_T2H.root --custom_conditions= --fileout=test_GEN-SIM-RECO-HLT.root --number=5 --mc --no_exec --eventcontent=RECOSIM --customise=HLTrigger/Configuration/CustomConfigs.L1THLT --python_filename=cmsDriver_my427p1HLTcfg-L1rerun_2011-07-25.py --processName=HLT
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/user/tsw/Boosted_Quark_To_ee_2TeV/std424Recon_2-00TeVu_v1/30a3840207bc3be6b1c52e0833382547/sigMCFile_GEN-SIM-std424RECO_v1_2011-07-05_9_1_T2H.root'),
    secondaryFileNames = cms.untracked.vstring('/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_9_1_IzJ.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_99_1_4c5.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_98_1_koO.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_97_1_Feh.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_96_1_5OF.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_95_1_JiK.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_94_1_yRL.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_93_1_tUU.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_92_1_W3z.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_91_1_wYG.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_90_1_KSc.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_8_1_VmT.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_89_1_jVT.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_88_1_e5Y.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_87_1_Xbv.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_86_1_eK7.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_85_1_GaU.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_84_1_iDN.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_83_1_YcQ.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_82_1_iHY.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_81_1_uib.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_80_1_WeB.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_7_1_BR5.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_79_1_Epy.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_78_1_ekj.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_77_1_GAb.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_76_1_HoN.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_75_1_v1M.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_74_1_r85.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_73_1_KQl.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_72_1_8I9.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_71_1_mdB.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_70_1_W5h.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_6_1_mND.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_69_1_Zgj.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_68_1_l0P.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_67_1_y2Y.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_66_1_TWr.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_65_1_gxZ.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_64_1_65U.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_63_1_HeD.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_62_1_uNP.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_61_1_JEM.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_60_1_gqm.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_5_1_urf.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_59_1_rTj.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_58_1_Tk5.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_57_1_P8w.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_56_1_KOD.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_55_1_74V.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_54_1_Y1K.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_53_1_vRN.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_52_1_HEp.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_51_1_c5x.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_50_1_ewd.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_4_1_9Fg.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_49_1_7ib.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_48_1_7xV.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_47_1_Qdi.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_46_1_N7N.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_45_1_xPh.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_44_1_Dav.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_43_1_agZ.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_42_1_149.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_41_1_N3l.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_40_1_HU2.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_3_1_SgU.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_39_1_7DZ.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_38_1_zzF.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_37_1_LHR.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_36_1_R7T.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_35_1_NLC.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_34_1_ZF0.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_33_1_kUn.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_32_1_ix3.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_31_1_brR.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_30_1_3Fm.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_2_1_t4k.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_29_1_mG9.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_28_1_4CP.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_27_1_fYY.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_26_1_Uuz.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_25_1_15E.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_24_1_hbw.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_23_1_CbM.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_22_1_H54.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_21_1_sJp.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_20_1_I5n.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_1_1_nIx.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_19_1_giI.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_18_1_aK2.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_17_1_oUX.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_16_1_3FH.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_15_1_Xu0.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_14_1_nPe.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_13_1_wvM.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_12_1_h0F.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_11_1_EZo.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_10_1_ExB.root',
        '/store/user/jacksonj/Boosted_Quark_To_ee_2TeV/Boosted_Quark_To_ee_2TeV/d2de6b25e2c771f4bdd5819b37d2253c/JJLHE_py_GEN_SIM_DIGI_L1_DIGI2RAW_100_1_z0L.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303.2.7 $'),
    annotation = cms.untracked.string('RelVal nevts:5'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('/opt/ppd/scratch/williams/test_GEN-SIM-RECO-HLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START42_V12::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1simulation_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.CustomConfigs
from HLTrigger.Configuration.CustomConfigs import L1THLT 

#call to customisation function L1THLT imported from HLTrigger.Configuration.CustomConfigs
process = L1THLT(process)

# End of customisation functions


###################################################################################################################################################################
# My own customisation, done in order to stop L1 emulator errors (associated with some RPC and CSC digis not being found when this config file is run with cmsRun)
# Lines taken from the 5e32 menu L1 emulator code /opt/ppd/scratch/williams/TriggerStudies/CMSSW_4_1_5/src/l1-lines-mc_2011-05-03.py on 26th July 2011
# (Details of the errors that this fixed can be found in 26th July 2011 entry of electronic logbook.)
#

### customize the L1 to run only Calo TPGs, GCT and GT
import L1Trigger.Configuration.L1Trigger_custom

process = L1Trigger.Configuration.L1Trigger_custom.customiseL1CaloAndGtEmulatorsFromRaw( process )
process = L1Trigger.Configuration.L1Trigger_custom.customiseResetPrescalesAndMasks( process )

### customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGctGtDigis( process )
