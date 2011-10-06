## /*****************************************************************************
##  * Project: CMS detector at the CERN
##  *
##  * Package: PhysicsTools/TagAndProbe
##  *
##  *
##  * Authors:
##  *
##  *   Kalanand Mishra, Fermilab - kalanand@fnal.gov
##  *
##  * Description:
##  *   - Produces tag & probe TTree for further analysis and computing efficiency
##  *
##  * History:
##  *   
##  * 
##  *****************************************************************************/


import FWCore.ParameterSet.Config as cms

##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################
## Following HLT paths are available in MC sample 
## "/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Fall10-START38_V12-v1"
## Please look into the trigger menu to find out the prescale for these paths.
## Process name used is: "HLT"
##
## HLT_Ele10_SW_L1R, HLT_Ele12_SW_TightEleId_L1R
## HLT_Ele12_SW_TightEleIdIsol_L1R
## HLT_Ele12_SW_TightEleIdIsol_NoDEtaInEE_L1R
## HLT_Ele17_SW_L1R
## HLT_Ele17_SW_CaloEleId_L1R
## HLT_Ele17_SW_LooseEleId_L1R
## HLT_Ele17_SW_EleId_L1R
## HLT_Ele22_SW_CaloEleId_L1R
## HLT_Ele40_SW_L1R
## HLT_DoubleEle4_SW_eeRes_L1R
## HLT_DoubleEle10_SW_L1R, 
## HLT_Photon10_Cleaned_L1R
## HLT_Photon15_Cleaned_L1R
## HLT_Photon20_NoHE_L1R
## HLT_Photon20_Cleaned_L1R
## HLT_Photon30_Cleaned_L1R
## HLT_Photon50_NoHE_L1R
## HLT_Photon50_NoHE_Cleaned_L1R
## HLT_DoublePhoton5_CEP_L1R
## HLT_DoublePhoton5_L1R,
## HLT_DoublePhoton10_L1R
## HLT_DoublePhoton15_L1R
## HLT_DoublePhoton17_L1R
################################################
## Following HLT paths are available in MC sample
## "/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia/Winter10-E7TeV_ProbDist_2010Data_BX156_START39_V8-v1"
## Process name used is: "REDIGI39X"
##
## HLT_Ele10_SW_L1R_v2
## HLT_Ele12_SW_TighterEleId_L1R_v2
## HLT_Ele17_SW_L1R_v2
## HLT_Ele17_SW_Isol_L1R_v2
## HLT_Ele17_SW_TighterEleIdIsol_L1R_v3
## HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2
## HLT_Ele22_SW_L1R_v2
## HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2
## HLT_Ele22_SW_TighterEleId_L1R_v3
## HLT_Ele32_SW_TighterEleId_L1R_v2
## HLT_DoubleEle4_SW_eeRes_L1R_v2
## HLT_DoubleEle5_SW_Upsilon_L1R_v2
## HLT_DoubleEle17_SW_L1R_v1
## HLT_Photon10_Cleaned_L1R
## HLT_Photon17_Isol_SC17HE_L1R_v1
## HLT_Photon20_NoHE_L1R
## HLT_Photon20_Cleaned_L1R
## HLT_Photon20_Isol_Cleaned_L1R_v1
## HLT_Photon22_SC22HE_L1R_v1
## HLT_Photon30_Cleaned_L1R
## HLT_Photon40_CaloId_Cleaned_L1R_v1
## HLT_Photon40_Isol_Cleaned_L1R_v1
## HLT_Photon50_Cleaned_L1R_v1
## HLT_Photon50_NoHE_L1R
## HLT_Photon70_Cleaned_L1R_v1
## HLT_Photon110_NoHE_Cleaned_L1R_v1
## HLT_DoublePhoton5_CEP_L1R_v3
## HLT_DoublePhoton17_SingleIsol_L1R_v1
## HLT_DoublePhoton22_L1R_v1
################################################
## Following electron/photon HLT paths are available in Run2010B data (first file)
## (replace "v1" with "v2", "v3" etc. for later runs).
## Please look into the trigger menu to find out the prescale for these paths.
## Process name used is: "HLT"
##
## HLT_Ele10_SW_L1R
## HLT_Ele12_SW_TightEleId_L1R
## HLT_Ele12_SW_TighterEleId_L1R_v1
## HLT_Ele12_SW_TighterEleIdIsol_L1R_v1
## HLT_Ele17_SW_L1R
## HLT_Ele17_SW_TightEleId_L1R
## HLT_Ele17_SW_TighterEleId_L1R_v1
## HLT_Ele17_SW_TightEleIdIsol_L1R_v1
## HLT_Ele17_SW_TighterEleIdIsol_L1R_v1
## HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1
## HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1
## HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1
## HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1  
## HLT_DoubleEle4_SW_eeRes_L1R
## HLT_DoubleEle15_SW_L1R_v1
## HLT_Photon10_Cleaned_L1R
## HLT_Photon15_Cleaned_L1R
## HLT_Photon17_SC17HE_L1R_v1
## HLT_Photon20_NoHE_L1R
## HLT_Photon20_Cleaned_L1R
## HLT_Photon30_Cleaned_L1R
## HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1
## HLT_Photon35_Isol_Cleaned_L1R_v1
## HLT_Photon50_Cleaned_L1R_v1
## HLT_Photon50_NoHE_L1R
## HLT_Photon70_NoHE_Cleaned_L1R_v1
## HLT_Photon100_NoHE_Cleaned_L1R_v1
## HLT_DoublePhoton5_CEP_L1R, HLT_DoublePhoton17_L1R
################################################

MC_flag = True
GLOBAL_TAG = 'GR_R_39X_V4::All'
if MC_flag:
    #GLOBAL_TAG = 'START38_V14::All'
    GLOBAL_TAG = 'START42_V11::All'
    
HLTPath = "HLT_Ele45_CaloIdVT_TrkIdT_v2"
HLTProcessName = "HLT"
if MC_flag:
    #HLTPath = "HLT_Ele17_SW_LooseEleId_L1R"
    #HLTPath = "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3"
    HLTPath = "HLT_Ele45_CaloIdVT_TrkIdT_v2"
    HLTProcessName = "HLT"

OUTPUT_FILE_NAME = "testTandPEffiTree.root"


ELECTRON_ET_CUT_MIN = 17.0
ELECTRON_COLL = "gsfElectrons"
ELECTRON_CUTS = "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####

PHOTON_COLL = "photons"
PHOTON_CUTS = "hadronicOverEm<0.15 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && ((isEB && sigmaIetaIeta<0.01) || (isEE && sigmaIetaIeta<0.03)) && (superCluster.energy*sin(superCluster.position.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
####

SUPERCLUSTER_COLL_EB = "hybridSuperClusters"
SUPERCLUSTER_COLL_EE = "multi5x5SuperClustersWithPreshower"
if MC_flag:
    SUPERCLUSTER_COLL_EB = "correctedHybridSuperClusters"
    SUPERCLUSTER_COLL_EE = "correctedMulti5x5SuperClustersWithPreshower"
SUPERCLUSTER_CUTS = "abs(eta)<2.5 && !(1.4442< abs(eta) <1.566) && et>" + str(ELECTRON_ET_CUT_MIN)


JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.6 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99" 
########################

##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|
##
process = cms.Process("TagProbe")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Geometry_cff")
process.GlobalTag.globaltag = GLOBAL_TAG
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

##   ____             _ ____                           
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
#           'file:/opt/ppd/scratch/williams/FEB5D990-917C-E011-A37F-003048D3C8D6.root'
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FEB5D990-917C-E011-A37F-003048D3C8D6.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FE6B9518-9D7C-E011-B606-003048D436F2.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FE0BF8DF-917C-E011-996D-002481E0DA60.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FE00E054-9A7C-E011-B016-0030487D5E81.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FCD68558-957C-E011-8C97-0030487D811F.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FCCA6184-8A7C-E011-BF91-003048CF632E.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FCBB1E8E-9B7C-E011-A1B9-003048C692FE.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FC5F2DAE-827C-E011-9685-003048D43944.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FC394557-977C-E011-A0E6-0030487D858D.root',
           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FA898555-977C-E011-893F-0030487D5DBF.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FA4B11B3-917C-E011-8256-003048D4399E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FA179724-8D7C-E011-B9E3-002481E0DDE8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FA132647-857C-E011-A29B-003048CF632E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F8660A71-857C-E011-9A8A-003048D4DCDE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F859835C-977C-E011-BA6A-0030487D70FD.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F4F4FFF7-977C-E011-8351-0030487D5E81.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F49781FE-8F7C-E011-B41D-00266CF32EAC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F45CA4DD-917C-E011-9556-003048F0E82C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F284B8A4-867C-E011-A5EA-0025901D4D6E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F2716D75-877C-E011-B6E6-003048D3C7DC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F209EEF7-8C7C-E011-8098-003048D479C0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F0ECFEB7-8F7C-E011-9402-0030487D5D91.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F07FE8C2-857C-E011-8DC8-003048C693B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/F04BF73D-807C-E011-A66B-0030487D864B.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EEB6CBBE-8A7C-E011-BEE7-003048D479C0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/ECB5137A-997C-E011-9463-003048D462BE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EC745015-947C-E011-9347-0025901D4D6E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EC5FA7AA-8A7C-E011-B646-002481E0DA96.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EC04A9F5-937C-E011-AF45-003048F0E18A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EAFD192F-7A7C-E011-AC6A-00266CF33208.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EAF84CEE-917C-E011-BB93-003048D3C886.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EAD02203-897C-E011-94AE-003048D462DC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/EA8BB1E7-917C-E011-9566-003048D436F2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E8EF5CEB-857C-E011-BA2B-003048C662B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E8D0483C-957C-E011-B920-00266CF33054.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E8A02CBD-977C-E011-9323-003048C693FA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E893E26A-917C-E011-B334-002481E0E450.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E877525D-917C-E011-8B61-0030487D5D7B.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E85C6BDA-7F7C-E011-B45B-003048C692F2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E6CD8A6C-877C-E011-B99E-003048C692F2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E66AFC4F-957C-E011-B4CC-00215AEDFC8E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E4E61F79-977C-E011-84EA-003048D4DEBC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E4AB6B2B-9D7C-E011-9DA3-002481E107A8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E49C95AF-857C-E011-9EF8-003048D43642.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E47A0080-997C-E011-BDB9-003048C693FA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E4684D61-917C-E011-AAB8-0030487D5E81.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E2A695A8-8E7C-E011-B61C-003048F0E5AA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E28E9C7F-827C-E011-9966-0030487D7B79.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E24C54C5-8F7C-E011-AF38-003048F0E18A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E2492B2B-927C-E011-B4DC-003048D4DCDE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/E20AC5C5-8D7C-E011-9D7A-003048C66BBE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DEF20FF4-627C-E011-AC48-003048D43656.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DECADA50-9A7C-E011-B1E6-003048D4610C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DE4AA5DE-8F7C-E011-AC4B-0030487D5E81.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DCCD3FF6-937C-E011-A9BD-003048D47912.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DC9D059B-957C-E011-B576-002481E0D480.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DC666B5E-947C-E011-87A9-00266CF33054.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DC61D547-827C-E011-88F2-0030487D5D95.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DC4847CB-5F7C-E011-9055-00237DE0BED6.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DC157990-917C-E011-8C13-003048C693F2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DA5DBA07-9D7C-E011-B0D5-003048D4DFAA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DA339B5D-5E7C-E011-B46F-0030487D5D53.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/DA0345E5-937C-E011-AD2D-003048F0E82C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D8651186-997C-E011-8863-002481E0DA96.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D424A048-8D7C-E011-8A0D-0030487D7B79.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D2B8BB09-8D7C-E011-BF21-003048D4393E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D26AD0CF-8C7C-E011-BC06-002481E1512E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D2110000-637C-E011-8078-00266CF33288.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D0F4480A-897C-E011-BF39-0025901D4C44.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D0F07410-927C-E011-99A0-003048D47912.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/D0ACB0C8-697C-E011-B9AC-003048C693E2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CEFBCEDC-857C-E011-86CF-003048F02CBA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CE99FA63-857C-E011-B176-003048C693B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CCABAFEE-627C-E011-AB1D-0030487D86C9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CC22E089-987C-E011-A0E4-0030487D8635.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CC1A9457-977C-E011-8266-003048C6929A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CC10FBA0-917C-E011-A044-003048D4363C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CACE62F8-937C-E011-A373-003048C693D2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CA85DBDA-8B7C-E011-8ECE-00266CF2718C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/CA285DE4-8A7C-E011-BE69-003048F0E1EA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C8C56D09-8D7C-E011-99CD-003048D437EC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C8AF2216-9B7C-E011-A928-0030487D83B9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C894AE3B-967C-E011-900C-003048C68A92.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C8346902-987C-E011-BC48-0030487D5E45.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C8343411-9B7C-E011-BB61-0030487D8635.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C818B232-637C-E011-9DBE-0030487F92E3.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C80BDEAD-917C-E011-BBF5-0030487E4EB9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C654F087-907C-E011-8FD0-003048C693D2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C6193A34-807C-E011-8497-0030487E4EB9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C4F97A0C-9D7C-E011-89CB-0030487E4EB7.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C4E885F8-937C-E011-8EF2-00237DDC5B9E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C4C07988-837C-E011-8BA5-002481E94C7A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C482A2F3-917C-E011-B279-003048D3CA06.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C47F6F3A-637C-E011-8B9A-003048D436BE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C474C0ED-8C7C-E011-A3A1-003048D462DC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C4631817-637C-E011-8FC0-003048D437F2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C447F326-807C-E011-94FF-0030487D5D95.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C2ECD149-957C-E011-84AE-0030487E4EB9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C284C86B-857C-E011-8A75-0030487D811F.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C20AFAED-917C-E011-99F7-002481E0DBE0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C0F654CF-8E7C-E011-AC02-003048F0EBB8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C09EA358-827C-E011-A16E-00266CF32E78.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/C0377127-9D7C-E011-B324-003048C693EC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BEAD1F58-7C7C-E011-BC40-003048D4DFA8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BE11872D-907C-E011-98C7-002481E0DC4C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BCF7F657-947C-E011-BCCC-0030487D811F.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BC8034A2-977C-E011-8B1E-0030487D5DB7.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BC45E207-8D7C-E011-9F35-003048D4363C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BAEA4202-677C-E011-B0D7-003048C693DA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BAD87F15-927C-E011-881A-003048D4DCDE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/BA7C5242-857C-E011-892C-0030487D811F.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B8FD8970-877C-E011-8D1F-003048F02CBA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B83EC810-907C-E011-AF64-003048D4399E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B69AEC20-9D7C-E011-9A71-002481E0D50C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B4CF9F70-997C-E011-B5B7-003048D43700.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B4A88853-907C-E011-9955-0025901D4C94.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B44B048F-5E7C-E011-9A3E-003048F0E5B4.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B437E6F6-8C7C-E011-8C4D-003048C693B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B2FCA669-917C-E011-A6FC-002481E0D974.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B2FB3D61-637C-E011-9449-003048F0E5B4.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B2DB0424-817C-E011-BA79-0030487E4EB9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B2A0A7CA-637C-E011-9404-0030487F1659.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B28F0B7D-867C-E011-BB87-003048C676E0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B09D4C96-987C-E011-8505-0030487D7103.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B06C01B7-8F7C-E011-B907-003048F0E5AA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B037A8EC-977C-E011-98B0-0030487D8541.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B0287140-857C-E011-969C-0030487D8661.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/B00A48D0-917C-E011-84A3-00215AD4D6E2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AEC6D8A6-8E7C-E011-9E28-002481E0DBE0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AEC40DC4-8E7C-E011-B258-003048C692E2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AEA27A0C-8D7C-E011-9AC8-003048CF6332.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AE129316-967C-E011-BC4B-003048F0E82C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AC89A712-897C-E011-B38C-003048D4363C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AC51B6F9-937C-E011-A477-003048D3C886.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AAF82C7F-987C-E011-AE37-0030487D5D3F.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AA55171A-9D7C-E011-9F1B-003048C69406.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/AA363683-987C-E011-9A8B-003048C693EC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A8FC7E1E-9B7C-E011-91AF-003048D436B4.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A8DB9B33-807C-E011-80EE-002481E1026A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A89FE864-8E7C-E011-90A6-003048D4DFA6.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A83CC5F1-917C-E011-BC63-003048D437EC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A6938F75-997C-E011-917F-0025901D4D6E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A60B64C4-807C-E011-AC48-003048D436EA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A4CCAB1B-9B7C-E011-933D-003048D436FE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A47EE358-887C-E011-8110-002481E0D5CE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A46ADBB1-917C-E011-8509-0030487D5D91.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A400808D-947C-E011-B583-002481E0DAB0.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A2EC361E-637C-E011-8F13-0025901D4C44.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A24EC9AC-8A7C-E011-AE06-003048F0E5A4.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A2137B0B-8D7C-E011-AA42-002481E0D5CE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/A03B5E44-877C-E011-B6A1-0030487D5D91.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9EA4471C-967C-E011-888F-0030487D5D3F.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9E241857-887C-E011-90C4-003048C693B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9CFC7668-977C-E011-8B78-003048D436FE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9C7457D6-8F7C-E011-866E-00266CF330B8.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9C556914-9B7C-E011-940A-003048D462BE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9C24BB75-647C-E011-BA2C-002481E76052.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9AF596A0-987C-E011-8C27-003048C692D6.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9AE91BD9-917C-E011-A834-0030487E4EC5.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9AE735B9-857C-E011-BFBD-0025901D4D6E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9AC5436D-637C-E011-B617-003048F0E3BA.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9A8A4A1F-907C-E011-B324-0030487E4EB9.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9A73A109-9D7C-E011-980D-003048D436FE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9A1570B4-637C-E011-930D-003048F0E51A.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/98E3100E-667C-E011-AD1D-003048D462DC.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/98B57C1E-637C-E011-99CB-0030487E510B.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/98A120D8-8C7C-E011-8E3A-00237DDC5B9E.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/986A8CB8-817C-E011-979E-002481E0D500.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/98568C8D-8E7C-E011-B1D8-003048D46122.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/985436E6-8D7C-E011-B206-003048D43980.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/981CE7C4-8B7C-E011-9F6F-002481E0D5CE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/96FBC79B-977C-E011-9D34-002481E0DE0C.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/969D1A6C-917C-E011-AC7E-003048F0E184.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9661D305-807C-E011-8228-003048C692E2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/94F08E97-9B7C-E011-B8CC-0030487D5EB3.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/94AB1618-8D7C-E011-BEB5-00237DDC5E96.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/94A02882-947C-E011-ABE8-003048D4DCDE.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/94965888-957C-E011-9B93-002481E15204.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9471E6B0-8B7C-E011-B380-0025901D4C44.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9426780D-9D7C-E011-8813-003048C69412.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/92739F25-817C-E011-B2CF-003048C692E2.root',
#           '/store/mc/Summer11/DYToEE_M-20_TuneZ2_7TeV-pythia6/AODSIM/PU_S3_START42_V11-v2/0000/9247D5F4-937C-E011-95CB-002481E0DA96.root'

        )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")

##   ____                         ____ _           _            
##  / ___| _   _ _ __   ___ _ __ / ___| |_   _ ___| |_ ___ _ __ 
##  \___ \| | | | '_ \ / _ \ '__| |   | | | | / __| __/ _ \ '__|
##   ___) | |_| | |_) |  __/ |  | |___| | |_| \__ \ ||  __/ |   
##  |____/ \__,_| .__/ \___|_|   \____|_|\__,_|___/\__\___|_|   
##  

##    ____      __ _____ _           _                   
##   / ___|___ / _| ____| | ___  ___| |_ _ __ ___  _ __  
##  | |  _/ __| |_|  _| | |/ _ \/ __| __| '__/ _ \| '_ \ 
##  | |_| \__ \  _| |___| |  __/ (__| |_| | | (_) | | | |
##   \____|___/_| |_____|_|\___|\___|\__|_|  \___/|_| |_|
##  
#  GsfElectron ################ 
process.goodElectrons = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( ELECTRON_CUTS )    
)      

##    _____ _           _                     ___    _ 
##   | ____| | ___  ___| |_ _ __ ___  _ __   |_ _|__| |
##   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \   | |/ _` |
##   | |___| |  __/ (__| |_| | | (_) | | | |  | | (_| |
##   |_____|_|\___|\___|\__|_|  \___/|_| |_| |___\__,_|
##   
# Electron ID  ######
process.PassingWP90 = process.goodElectrons.clone()
process.PassingWP90.cut = cms.string(
    process.goodElectrons.cut.value() +
    " && (gsfTrack.trackerExpectedHitsInner.numberOfHits==0 && !(-0.02<convDist<0.02 && -0.02<convDcot<0.02))"
    " && ((isEB"
    " && ( dr03TkSumPt/p4.Pt <0.12 && dr03EcalRecHitSumEt/p4.Pt < 0.09 && dr03HcalTowerSumEt/p4.Pt  < 0.1 )"
    " && (sigmaIetaIeta<0.01)"
    " && ( -0.8<deltaPhiSuperClusterTrackAtVtx<0.8 )"
    " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
    " && (hadronicOverEm<0.12)"
    ")"
    " || (isEE"
    " && ( dr03TkSumPt/p4.Pt <0.05 && dr03EcalRecHitSumEt/p4.Pt < 0.06 && dr03HcalTowerSumEt/p4.Pt  < 0.03 )"
    " && (sigmaIetaIeta<0.03)"
    " && ( -0.7<deltaPhiSuperClusterTrackAtVtx<0.7 )"
    " && ( -0.009<deltaEtaSuperClusterTrackAtVtx<0.009 )"
    " && (hadronicOverEm<0.05) "
    "))"
    ) 

process.PassingWP80 = process.goodElectrons.clone()
process.PassingWP80.cut = cms.string(
    process.goodElectrons.cut.value() +
    " && (gsfTrack.trackerExpectedHitsInner.numberOfHits==0 && !(-0.02<convDist<0.02 && -0.02<convDcot<0.02))"
    " && ((isEB"
    " && ( dr03TkSumPt/p4.Pt <0.09 && dr03EcalRecHitSumEt/p4.Pt < 0.07 && dr03HcalTowerSumEt/p4.Pt  < 0.1 )"
    " && (sigmaIetaIeta<0.01)"
    " && ( -0.06<deltaPhiSuperClusterTrackAtVtx<0.06 )"
    " && ( -0.004<deltaEtaSuperClusterTrackAtVtx<0.004 )"
    " && (hadronicOverEm<0.04)"
    ")"
    " || (isEE"
    " && ( dr03TkSumPt/p4.Pt <0.04 && dr03EcalRecHitSumEt/p4.Pt < 0.05 && dr03HcalTowerSumEt/p4.Pt  < 0.025 )"
    " && (sigmaIetaIeta<0.03)"
    " && ( -0.03<deltaPhiSuperClusterTrackAtVtx<0.03 )"
    " && ( -0.007<deltaEtaSuperClusterTrackAtVtx<0.007 )"
    " && (hadronicOverEm<0.025) "
    "))"
    ) 
                         
##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   
# Trigger  ##################
process.PassingHLT = cms.EDProducer("trgMatchedGsfElectronProducer",    
    InputProducer = cms.InputTag( ELECTRON_COLL ),                          
    hltTags = cms.VInputTag(cms.InputTag(HLTPath,"", HLTProcessName)),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",HLTProcessName),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults","",HLTProcessName)   
)

##    _____      _                        _  __     __             
##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
##   
## Here we show how to use a module to compute an external variable
## process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
## ak5PFResidual.useCondDB = False

process.superClusterDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("goodSuperClusters"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)
process.JetMultiplicityInSCEvents = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("goodSuperClusters"),
    objects = cms.InputTag(JET_COLL),
    objectSelection = cms.string(JET_CUTS + " && pt > 20.0"),
)

process.GsfDRToNearestJet = process.superClusterDRToNearestJet.clone()
process.GsfDRToNearestJet.probes = cms.InputTag( ELECTRON_COLL )
process.JetMultiplicityInGsfEvents = process.JetMultiplicityInSCEvents.clone()
process.JetMultiplicityInGsfEvents.probes = cms.InputTag( ELECTRON_COLL )

process.ext_ToNearestJet_sequence = cms.Sequence(
    #process.ak5PFResidual +   
    process.GsfDRToNearestJet +
    process.JetMultiplicityInGsfEvents
    )

##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/
## 
process.Tag = process.PassingHLT.clone()
process.Tag.InputProducer = cms.InputTag( "PassingWP80" )

process.ele_sequence = cms.Sequence(
    process.goodElectrons +
    process.PassingWP90 +
    process.PassingWP80 +
    process.PassingHLT +
    process.Tag
    )


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
process.tagGsf = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("Tag goodElectrons"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("40 < mass < 1000"),
)

process.tagWP90 = process.tagGsf.clone()
process.tagWP90.decay = cms.string("Tag PassingWP90")

process.allTagsAndProbes = cms.Sequence(
    process.tagGsf +
    process.tagWP90
)

##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("Tag"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchGsf = process.McMatchTag.clone()
process.McMatchGsf.src = cms.InputTag("goodElectrons")
process.McMatchWP90 = process.McMatchTag.clone()
process.McMatchWP90.src = cms.InputTag("PassingWP90")
    
process.mc_sequence = cms.Sequence(
   process.McMatchTag +
   process.McMatchGsf +
   process.McMatchWP90 
)

############################################################################
##    _____           _       _ ____            _            _   _  ____  ##
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| ##
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  ##
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | ##
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| ##
##              |___/                                                     ##
##                                                                        ##
############################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/
##
## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category
ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    pt  = cms.string("pt"),
    phi  = cms.string("phi"),
    et  = cms.string("et"),
    e  = cms.string("energy"),
    p  = cms.string("p"),
    px  = cms.string("px"),
    py  = cms.string("py"),
    pz  = cms.string("pz"),
    theta  = cms.string("theta"),    
    vx     = cms.string("vx"),
    vy     = cms.string("vy"),
    vz     = cms.string("vz"),
    rapidity  = cms.string("rapidity"),
    mass  = cms.string("mass"),
    mt  = cms.string("mt"),    
)   

ProbeVariablesToStore = cms.PSet(
    probe_gsfEle_eta = cms.string("eta"),
    probe_gsfEle_pt  = cms.string("pt"),
    probe_gsfEle_phi  = cms.string("phi"),
    probe_gsfEle_et  = cms.string("et"),
    probe_gsfEle_e  = cms.string("energy"),
    probe_gsfEle_p  = cms.string("p"),
    probe_gsfEle_px  = cms.string("px"),
    probe_gsfEle_py  = cms.string("py"),
    probe_gsfEle_pz  = cms.string("pz"),
    probe_gsfEle_theta  = cms.string("theta"),    
    probe_gsfEle_charge = cms.string("charge"),
    probe_gsfEle_rapidity  = cms.string("rapidity"),
    probe_gsfEle_missingHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
    probe_gsfEle_convDist = cms.string("convDist"),
    probe_gsfEle_convDcot = cms.string("convDcot"),
    probe_gsfEle_convRadius = cms.string("convRadius"),        
    probe_gsfEle_hasValidHitInFirstPixelBarrel = cms.string("gsfTrack.hitPattern.hasValidHitInFirstPixelBarrel"),
    ## super cluster quantities
    probe_sc_energy = cms.string("superCluster.energy"),
    probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
    probe_sc_x      = cms.string("superCluster.x"),
    probe_sc_y      = cms.string("superCluster.y"),
    probe_sc_z      = cms.string("superCluster.z"),
    probe_sc_eta    = cms.string("superCluster.eta"),
    probe_sc_theta  = cms.string("superClusterPosition.theta"),   
    probe_sc_phi    = cms.string("superCluster.phi"),
    probe_sc_size   = cms.string("superCluster.size"), # number of hits
    ## isolation 
    probe_gsfEle_trackiso = cms.string("dr03TkSumPt"),
    probe_gsfEle_ecaliso  = cms.string("dr03EcalRecHitSumEt"),
    probe_gsfEle_hcaliso  = cms.string("dr03HcalTowerSumEt"),
    ## classification, location, etc.    
    probe_gsfEle_classification = cms.string("classification"),
    probe_gsfEle_numberOfBrems  = cms.string("numberOfBrems"),     
    probe_gsfEle_bremFraction   = cms.string("fbrem"),
    probe_gsfEle_mva            = cms.string("mva"),        
    probe_gsfEle_deltaEta       = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    probe_gsfEle_deltaPhi       = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    probe_gsfEle_deltaPhiOut    = cms.string("deltaPhiSeedClusterTrackAtCalo"),
    probe_gsfEle_deltaEtaOut    = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    probe_gsfEle_isEB           = cms.string("isEB"),
    probe_gsfEle_isEE           = cms.string("isEE"),
    probe_gsfEle_isGap          = cms.string("isGap"),
    ## Hcal energy over Ecal Energy
    probe_gsfEle_HoverE         = cms.string("hcalOverEcal"),    
    probe_gsfEle_EoverP         = cms.string("eSuperClusterOverP"),
    probe_gsfEle_eSeedClusterOverP = cms.string("eSeedClusterOverP"),    
    ## Cluster shape information
    probe_gsfEle_sigmaEtaEta  = cms.string("sigmaEtaEta"),
    probe_gsfEle_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
    probe_gsfEle_e1x5               = cms.string("e1x5"),
    probe_gsfEle_e2x5Max            = cms.string("e2x5Max"),
    probe_gsfEle_e5x5               = cms.string("e5x5"),
    ## is ECAL driven ? is Track driven ?
    probe_gsfEle_ecalDrivenSeed     = cms.string("ecalDrivenSeed"),
    probe_gsfEle_trackerDrivenSeed  = cms.string("trackerDrivenSeed")
)


TagVariablesToStore = cms.PSet(
    gsfEle_eta = cms.string("eta"),
    gsfEle_pt  = cms.string("pt"),
    gsfEle_phi  = cms.string("phi"),
    gsfEle_et  = cms.string("et"),
    gsfEle_e  = cms.string("energy"),
    gsfEle_p  = cms.string("p"),
    gsfEle_px  = cms.string("px"),
    gsfEle_py  = cms.string("py"),
    gsfEle_pz  = cms.string("pz"),
    gsfEle_theta  = cms.string("theta"),    
    gsfEle_charge = cms.string("charge"),
    gsfEle_rapidity  = cms.string("rapidity"),
    gsfEle_missingHits = cms.string("gsfTrack.trackerExpectedHitsInner.numberOfHits"),
    gsfEle_convDist = cms.string("convDist"),
    gsfEle_convDcot = cms.string("convDcot"),
    gsfEle_convRadius = cms.string("convRadius"),     
    gsfEle_hasValidHitInFirstPixelBarrel = cms.string("gsfTrack.hitPattern.hasValidHitInFirstPixelBarrel"),
    ## isolation 
    gsfEle_trackiso = cms.string("dr03TkSumPt"),
    gsfEle_ecaliso  = cms.string("dr03EcalRecHitSumEt"),
    gsfEle_hcaliso  = cms.string("dr03HcalTowerSumEt"),
    ## classification, location, etc.    
    gsfEle_classification = cms.string("classification"),
    gsfEle_numberOfBrems  = cms.string("numberOfBrems"),     
    gsfEle_bremFraction   = cms.string("fbrem"),
    gsfEle_mva            = cms.string("mva"),        
    gsfEle_deltaEta       = cms.string("deltaEtaSuperClusterTrackAtVtx"),
    gsfEle_deltaPhi       = cms.string("deltaPhiSuperClusterTrackAtVtx"),
    gsfEle_deltaPhiOut    = cms.string("deltaPhiSeedClusterTrackAtCalo"),
    gsfEle_deltaEtaOut    = cms.string("deltaEtaSeedClusterTrackAtCalo"),
    gsfEle_isEB           = cms.string("isEB"),
    gsfEle_isEE           = cms.string("isEE"),
    gsfEle_isGap          = cms.string("isGap"),
    ## Hcal energy over Ecal Energy
    gsfEle_HoverE         = cms.string("hcalOverEcal"),    
    gsfEle_EoverP         = cms.string("eSuperClusterOverP"),
    gsfEle_eSeedClusterOverP = cms.string("eSeedClusterOverP"),  
    ## Cluster shape information
    gsfEle_sigmaEtaEta  = cms.string("sigmaEtaEta"),
    gsfEle_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
    gsfEle_e1x5               = cms.string("e1x5"),
    gsfEle_e2x5Max            = cms.string("e2x5Max"),
    gsfEle_e5x5               = cms.string("e5x5"),
    ## is ECAL driven ? is Track driven ?
    gsfEle_ecalDrivenSeed     = cms.string("ecalDrivenSeed"),
    gsfEle_trackerDrivenSeed  = cms.string("trackerDrivenSeed")
)

CommonStuffForGsfElectronProbe = cms.PSet(
    variables = cms.PSet(ProbeVariablesToStore),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    pairVariables =  cms.PSet(ZVariablesToStore),
    pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagVariablesToStore),
    tagFlags     =  cms.PSet(
          passingGsf = cms.InputTag("goodElectrons"),
          isWP90 = cms.InputTag("PassingWP90"),         
          passingHLT = cms.InputTag("PassingHLT")     
    ),    
)

if MC_flag:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_pt  = cms.string("pt"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_p  = cms.string("p"),
        probe_px  = cms.string("px"),
        probe_py  = cms.string("py"),
        probe_pz  = cms.string("pz"),
        probe_theta  = cms.string("theta"),    
        probe_vx     = cms.string("vx"),
        probe_vy     = cms.string("vy"),
        probe_vz     = cms.string("vz"),   
        probe_charge = cms.string("charge"),
        probe_rapidity  = cms.string("rapidity"),    
        probe_mass  = cms.string("mass"),
        probe_mt  = cms.string("mt"),    
        ),
        mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),      
        )
else:
     mcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )



##   ____      __       __    ___                 ___    _ 
##  / ___|___ / _|      \ \  |_ _|___  ___       |_ _|__| |
## | |  _/ __| |_   _____\ \  | |/ __|/ _ \       | |/ _` |
## | |_| \__ \  _| |_____/ /  | |\__ \ (_) |  _   | | (_| |
##  \____|___/_|        /_/  |___|___/\___/  ( ) |___\__,_|
##                                           |/            
##  gsf electron --> isolation, electron id  etc.
process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeProducer",
    mcTruthCommonStuff, CommonStuffForGsfElectronProbe,                        
    tagProbePairs = cms.InputTag("tagGsf"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet(
        probe_isWP90 = cms.InputTag("PassingWP90"),
        probe_passingHLT = cms.InputTag("PassingHLT")        
    ),
    probeMatches  = cms.InputTag("McMatchGsf"),
    allProbes     = cms.InputTag("goodElectrons")
)
process.GsfElectronToId.variables.probe_dRjet = cms.InputTag("GsfDRToNearestJet")
process.GsfElectronToId.variables.probe_nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.GsfElectronToId.tagVariables.dRjet = cms.InputTag("GsfDRToNearestJet")
process.GsfElectronToId.tagVariables.nJets = cms.InputTag("JetMultiplicityInGsfEvents")
process.GsfElectronToId.pairVariables.costheta = cms.InputTag("CSVarsTagGsf","costheta")
process.GsfElectronToId.pairVariables.sin2theta = cms.InputTag("CSVarsTagGsf","sin2theta")
process.GsfElectronToId.pairVariables.tanphi = cms.InputTag("CSVarsTagGsf","tanphi")

##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|
##
##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
if MC_flag:
    HLTmcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_flag),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_flag),
        checkMotherInUnbiasEff = cms.bool(MC_flag),
        mcVariables = cms.PSet(
          probe_eta = cms.string("eta"),
          probe_phi  = cms.string("phi"),
          probe_et  = cms.string("et"),
          probe_charge = cms.string("charge"),
        ),
        mcFlags     =  cms.PSet(
          probe_flag = cms.string("pt>0")
        ),      
        )
else:
     HLTmcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )

##  WP95 --> HLT
process.WP90ToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,                                
    variables = cms.PSet(
      probe_gsfEle_eta = cms.string("eta"),
      probe_gsfEle_phi  = cms.string("phi"),
      probe_gsfEle_et  = cms.string("et"),
      probe_gsfEle_charge = cms.string("charge"),
      probe_sc_et    = cms.string("superCluster.energy*sin(superClusterPosition.theta)"),    
      probe_sc_eta    = cms.string("superCluster.eta"), 
      probe_sc_phi    = cms.string("superCluster.phi"),
      probe_gsfEle_isEB           = cms.string("isEB"),
      probe_gsfEle_isEE           = cms.string("isEE"),
      probe_gsfEle_isGap          = cms.string("isGap"),
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (False),
    addEventVariablesInfo   =  cms.bool (False),                                                        
    tagProbePairs = cms.InputTag("tagWP90"),
    arbitration   = cms.string("Random2"),
    flags = cms.PSet( 
        probe_passingHLT = cms.InputTag("PassingHLT")        
    ),
    probeMatches  = cms.InputTag("McMatchWP90"),
    allProbes     = cms.InputTag("PassingWP90")
)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

##    _   _        ___    _  
##   | | | |      |_ _|__| | 
##   | |_| |       | |/ _` | 
##   |  _  |       | | (_| | 
##   |_| |_|EEP   |___\__,_| efficiency section
##

gsfElectronCollection = "gsfElectrons"

## Defining the strings for access to many HEEP variables ...
#  (exact methods used to calculate these variables are same methods as used in UserCode/SHarper/HEEPAnalyzer/interface/HEEPEle.h HEAD version (=1.10) on 29/07/11)
var_scEt   = "((superCluster.position.rho/superCluster.position.r)*caloEnergy)"
var_scEta  = "(superCluster.eta)"
var_scPhi  = "(superCluster.phi)"
var_heepEt = "(caloEnergy*cos(p4.theta))"
var_eta    = "(p4.eta)"
var_phi    = "(p4.phi)"

var_isEcalDriven = "ecalDrivenSeed"
var_dEtaIn       = "deltaEtaSuperClusterTrackAtVtx"
var_dPhiIn       = "deltaPhiSuperClusterTrackAtVtx"
var_hOverE       = "hadronicOverEm"
var_sigmaIEtaIEta= "sigmaIetaIeta"
var_e2x5Over5x5  = "(? e5x5!=0 ? e1x5/e5x5 : 0.0)"
var_e1x5Over5x5  = "(? e5x5!=0 ? e2x5Max/e5x5 : 0.0)"

var_isolEmHadDepth1 = "(dr03EcalRecHitSumEt+dr03HcalDepth1TowerSumEt)"
var_isolHadDepth2   = "dr03HcalDepth2TowerSumEt"
var_isolPtTrks      = "dr03TkSumPt"

var_epIn  = "eSuperClusterOverP"
var_epOut = "eSeedClusterOverPout"

##
## Defining variables to be stored for tag ...
varsToStore_heepTags = cms.PSet(
    probe_scEt   = cms.string(var_scEt),
    probe_scEta  = cms.string(var_scEta),
    probe_scPhi  =cms.string(var_scPhi),
    probe_heepEt =cms.string(var_heepEt),
    probe_eta    =cms.string(var_eta),
    probe_phi    =cms.string(var_phi),
    #
    probe_isEcalDriven  = cms.string(var_isEcalDriven),
    probe_dEtaIn        = cms.string(var_dEtaIn),
    probe_dPhiIn        = cms.string(var_dPhiIn),
    probe_hOverE        = cms.string(var_hOverE),
    probe_sigmaIEtaIEta = cms.string(var_sigmaIEtaIEta),
    probe_e2x5Over5x5   = cms.string(var_e2x5Over5x5),
    probe_e1x5Over5x5   = cms.string(var_e1x5Over5x5),
    #
    probe_isolEmHadDepth1= cms.string(var_isolEmHadDepth1),
    probe_isolHadDepth2  = cms.string(var_isolHadDepth2),
    probe_isolPtTrks     = cms.string(var_isolPtTrks),
    probe_epIn  = cms.string(var_epIn),
    probe_epOut = cms.string(var_epOut)
)
## Defining variables to be stored for probe ...
varsToStore_heepProbes = varsToStore_heepTags

## Defining variables to be stored for the T&P pair - i.e. the Z (candidate) ...
varsToStore_heepTnPPair = cms.PSet(
    gsfEta = cms.string("eta"),
    gsfPt  = cms.string("pt"),
    gsfPhi  = cms.string("phi"),
    gsfEt  = cms.string("et"),
    gsfE  = cms.string("energy"),
    gsfP  = cms.string("p"),
    gsfPx  = cms.string("px"),
    gsfPy  = cms.string("py"),
    gsfPz  = cms.string("pz"),
    gsfTheta  = cms.string("theta"),    
    gsfVx     = cms.string("vx"),
    gsfVy     = cms.string("vy"),
    gsfVz     = cms.string("vz"),
    gsfRapidity  = cms.string("rapidity"),
    gsfMass  = cms.string("mass"),
    gsfMt  = cms.string("mt"),    
)

##
## Defining PhysicsCutsParser strings for common cuts ...
cut_scIsBarrel = "( abs("+var_scEta+")<1.442 )"
cut_scIsEndcap = "( abs("+var_scEta+")>1.560 && abs("+var_scEta+")<2.5 )"

cut_heepBarrel = ("(" + cut_scIsBarrel + 
                           " && (" + var_heepEt       + ">35.0)" + 
                           " && (" + var_isEcalDriven + "==1)" + 
                           " && (abs(" + var_dEtaIn   + ")<0.005)" + 
                           " && (abs(" + var_dPhiIn   + ")<0.09)" + 
                           " && (" + var_hOverE       + "<0.05)" + 
                           " && (("+ var_e2x5Over5x5  + ">0.94) || ("+var_e1x5Over5x5+">0.83))" + 
                           " && (" + var_isolEmHadDepth1+ "<(2.0+0.03*"+var_heepEt+") )" + 
                           " && (" + var_isolPtTrks   + "<7.5)" + 
                           ")")
cut_heepEndcap = ("(" + cut_scIsEndcap + 
                           " && (" + var_heepEt       + ">40.0)" + 
                           " && (" + var_isEcalDriven + "==1)" + 
                           " && (abs(" + var_dEtaIn   + ")<0.007)" + 
                           " && (abs(" + var_dPhiIn   + ")<0.09)" + 
                           " && (" + var_hOverE       + "<0.05)" +
                           " && (" + var_sigmaIEtaIEta+ "<0.03)" + 
                           " && ( " + "(("+var_heepEt+"<50.0)&&("+var_isolEmHadDepth1+"<2.5)) || " + 
                                      "(("+var_heepEt+">=50.0)&&("+var_isolEmHadDepth1+"< (2.5+0.03*("+var_heepEt+"-50.0)) ))" + " )" +
                           #" && (? "+var_heepEt+"<50.0 ? " + var_isolEmHadDepth1 + "<2.5 : " + var_isolEmHadDepth1 + "<2.5)" +
                           #                  var_isolEmHadDepth1+ "<(2.5+0.03*("+var_heepEt+"-50.0)) )" + 
                           " && (" + var_isolHadDepth2+ "<0.5)" + 
                           " && (" + var_isolPtTrks   + "<15.0)" + 
                           ")")

## Defining the probe electron collection ...
heepProbeCuts = "(" + var_scEt + ">35.0)"
process.heepProbeEles = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( gsfElectronCollection ),
    cut = cms.string( heepProbeCuts )    
)

## Defining the AllHEEP 'passing' probe electron collection ...
process.heepProbeElesPassingAllHEEP = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( gsfElectronCollection ),
    cut = cms.string( heepProbeCuts + " && ( " + cut_heepBarrel + " || " + cut_heepEndcap + " )" )    
#    cut = cms.string( heepProbeCuts + " && ( " + cut_heepBarrel + " )" )    
)

## Defining the tag electron collection ... (Barrel HEEP ele with e/p<1.5)
heepTagCuts = ("( " + cut_heepBarrel +
               " && (" + var_epIn + "<1.5)" + 
               " )")
process.heepTagEles = cms.EDFilter("GsfElectronRefSelector",
    src = cms.InputTag( gsfElectronCollection ),
    cut = cms.string( heepTagCuts )    
)

## Defining the T & P pairs ...
process.heepTagProbePairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("heepTagEles heepProbeEles"), # charge coniugate states are implied
    checkCharge = cms.bool(False),                           
    cut   = cms.string("60 < mass < 120"), ## TODO - Should change to HEEP definition of p4 for each electron in time,
                                           ## but using non-HEEP definition at the moment is fine for the mass cut at the moment.
)

##
## (FINALLY) Actually making the tag-probe tree ...
process.heepTagProbeTreeProducer = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## TODO -- Check that (as is done in official HEEP effi measurement), events with two tags are counted twice
    # General flags
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (True),
    addEventVariablesInfo   =  cms.bool (True),
    ## *** PROBE ELECTRON ***
    # Define the probe electrons (a 'GsfElectronRefSelector' product)
    allProbes     = cms.InputTag("heepProbeEles"),
    # Define which probe variables should be stored
    variables = cms.PSet(varsToStore_heepProbes), 
    # Define probe electrons that pass cuts whose effi is being measured ('GsfElectronRefSelector' product) - sets up boolean branches in TnP tree
    flags = cms.PSet(    
        probe_passesAllHEEP = cms.InputTag("heepProbeElesPassingAllHEEP"),
        #probe_passesHEEPNoIso = cms.InputTag("heepProbeElesPassingHEEPNoIso")        
    ),
    ## *** TAG ELECTRON ***
    # Define which tag variables should be stored, and which cut-passing flags
    tagVariables   =  cms.PSet(varsToStore_heepTags),
    tagFlags     =  cms.PSet(
        probe_passesAllHEEP = cms.InputTag("heepProbeElesPassingAllHEEP"),
        #probe_passesHEEPNoIso = cms.InputTag("heepProbeElesPassingHEEPNoIso")
    ),
    ## *** TAG-PROBE PAIRS ***
    # Define the tag-probe pairs            
    tagProbePairs = cms.InputTag("heepTagProbePairs"),
    # Define which tag-probe pair variables (and cut-passing flags) should be stored
    pairVariables = cms.PSet(varsToStore_heepTnPPair),
    pairFlags     = cms.PSet( mass60to120 = cms.string("60 < mass < 120") ),
    # Define what happens when there are multiple TnP pairs for a given tag electron
    arbitration   = cms.string("Random2"),  # TODO - Check what value Laurent uses for arbitration
    #----------------------------------------------------
    ## *** MC MATCHING ***
    ## TODO - Check what Laurent does for MC matching. I will temporarily say that files aren't MC
    isMC = cms.bool(False)
    #probeMatches  = cms.InputTag("McMatchGsf"),
    # Set up use of MC truth                  
    #mcTruthCommonStuff
)

process.heepEffiPaths = cms.Sequence( process.heepProbeEles + process.heepProbeElesPassingAllHEEP + process.heepTagEles + process.heepTagProbePairs + process.heepTagProbeTreeProducer)

#################################################################################################################################################
#################################################################################################################################################

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##

process.tree_sequence = cms.Sequence(  process.GsfElectronToId + process.WP90ToHLT  )     

if MC_flag:
    process.tagAndProbe = cms.Path(
        process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence + 
        process.tree_sequence + 
        process.heepEffiPaths
        )
else:
    process.tagAndProbe = cms.Path(
        process.ele_sequence + 
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.tree_sequence
        )
    
process.TFileService = cms.Service(
    "TFileService", fileName = cms.string( OUTPUT_FILE_NAME )
    )
