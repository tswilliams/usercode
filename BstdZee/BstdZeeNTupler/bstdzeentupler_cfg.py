import FWCore.ParameterSet.Config as cms

#####################################################
## Useful function for local running...

def DataFileLocationAdaptor(fileLocation):
   fileLocation.replace("\n","")
   filePrefix="file:"
   if(fileLocation.startswith("/pnfs/")==True):
      filePrefix="dcap://heplnx209.pp.rl.ac.uk:22125"
   elif(fileLocation.startswith("/store/")==True):
      filePrefix=""
   return filePrefix + fileLocation



#####################################################
## General setup/load lines, defining nevts, input & output etc.

process = cms.Process("Demo")

## Loading Geometry modules for creation of transient geometry information used in determining eta and phi values of recHits ...
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START42_V13::All'

process.load("FWCore.MessageService.MessageLogger_cfi")
## Lines from P. Lenzi related to reading LHE info
process.MessageLogger.categories=cms.untracked.vstring('FwkJob'
                                                           ,'FwkReport'
                                                           ,'FwkSummary'
                                                           ,'Root_NoDictionary'
                                                           ,'Generator'
                                                           ,'LHEInterface'
                                                           )
process.MessageLogger.cerr.Generator = cms.untracked.PSet(limit = cms.untracked.int32(0))
process.MessageLogger.cerr.LHEInterface = cms.untracked.PSet(limit = cms.untracked.int32(10000))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Reading in a list of datafile locations from a file...
fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/Boosted_Quark_To_ee_2TeV/tsw-bstd412Recon_2-00TeVu_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"
f_DatafileList = open(fileName_DatafileList)
datafilesList = f_DatafileList.readlines()
f_DatafileList.close()
datafileLocations = map(DataFileLocationAdaptor,datafilesList)

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('/store/data/Run2011B/Photon/AOD/PromptReco-v1/000/180/252/6291CDCD-0B05-E111-8CDB-003048F117EC.root')
#                            fileNames = cms.untracked.vstring('/store/mc/Fall11/DYToLL_M-50_1jEnh2_2jEnh35_3jEnh40_4jEnh50_7TeV-sherpa/AODSIM/PU_S6_START42_V14B-v1/0000/0028E808-D407-E111-AD5B-0026189438D4.root')
                            fileNames = cms.untracked.vstring('/store/mc/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/0007A683-23F6-E011-9703-90E6BA0D09DC.root')
)

#Defining the output file to store the histograms/NTuples in...
process.TFileService = cms.Service("TFileService", fileName=cms.string('mcNTuple_42Xv1h.root'))



#####################################################
## Modified muon track reconstruction [OLD]

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



###################################################################################################
## Code for modified isolation values ...
# ModIso: 1) Creating a collection of HEEP ID (no iso) passing electrons ...
#setting up the producer to make the HEEP ID value map
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                eleLabel = cms.InputTag("gsfElectrons"),
                                barrelCuts = cms.PSet(heepBarrelCuts),
                                endcapCuts = cms.PSet(heepEndcapCuts)
                                )
process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits")
process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits")

process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer",
                                         cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")
                                         )

# ModIso: 2) Calculating the modified iso. values from IsoDeposits

from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositTk_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositEcalFromHits_cff import *
from RecoEgamma.EgammaIsolationAlgos.eleIsoDepositHcalFromTowers_cff import *
from RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi import *

process.eleIsoDepositTk = eleIsoDepositTk
process.ecalSeverityLevel = ecalSeverityLevel  # This line is required so that "EcalSeverityLevelAlgoRcd" record can be found in the EventSetup by CandIsoDepositProducer/eleIsoDepositEcalFromHits module
process.eleIsoDepositEcalFromHits = eleIsoDepositEcalFromHits
process.eleIsoDepositHcalDepth1FromTowers = eleIsoDepositHcalDepth1FromTowers

process.eleIsoDepositTk.src = "gsfElectrons"
process.eleIsoDepositEcalFromHits.src = "gsfElectrons"
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB")
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE")
process.eleIsoDepositHcalDepth1FromTowers.src = "gsfElectrons"

process.stdEleIsoFromDepsTk = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositTk"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('RectangularEtaPhiVeto(-0.015,0.015,-0.3,0.3)',
                           'Threshold(0.7)'),
       skipDefaultVeto = cms.bool(True)
   ))
)
process.modEleIsoFromDepsTk = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositTk"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('RectangularEtaPhiVeto(-0.015,0.015,-0.3,0.3)',
                           'ThresholdFromTransverse(0.7)',
                           'heepIdNoIsoEles:RectangularEtaPhiVeto(-0.015,0.015,-0.3,0.3)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

process.stdEleIsoFromDepsEcal = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositEcalFromHits"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('NumCrystalVeto(3.0)',
                           'EcalBarrel:NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThresholdFromTransverse(0.08)'),
       skipDefaultVeto = cms.bool(True)
   ))
)
                                             
process.modEleIsoFromDepsEcal = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       mode = cms.string('sum'),
       src = cms.InputTag("eleIsoDepositEcalFromHits"),
       weight = cms.string('1'),
       deltaR = cms.double(0.3),
       vetos = cms.vstring('NumCrystalVeto(3.0)',
                           'EcalBarrel:NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThresholdFromTransverse(0.08)',
                           'heepIdNoIsoEles:NumCrystalVeto(3.0)',
                           'heepIdNoIsoEles:NumCrystalEtaPhiVeto(1.5,17.25)'),
       skipDefaultVeto = cms.bool(True)
   ))
)

process.stdEleIsoFromDepsHcalD1 = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       src = cms.InputTag("eleIsoDepositHcalDepth1FromTowers"),
       deltaR = cms.double(0.3),
       weight = cms.string('1'),
       vetos = cms.vstring('0.15'),
       skipDefaultVeto = cms.bool(True),
       mode = cms.string('sum')
   ))
) 
process.modEleIsoFromDepsHcalD1 = cms.EDProducer("CandIsolatorFromDeposits",
   deposits = cms.VPSet(cms.PSet(
       src = cms.InputTag("eleIsoDepositHcalDepth1FromTowers"),
       deltaR = cms.double(0.3),
       weight = cms.string('1'),
       vetos = cms.vstring('0.15',
                           'heepIdNoIsoEles:0.15'),
       skipDefaultVeto = cms.bool(True),
       mode = cms.string('sum')
   ))
) 
process.recalculateIsoVars = cms.Sequence(process.eleIsoDepositTk * process.stdEleIsoFromDepsTk * process.modEleIsoFromDepsTk *
                                          process.eleIsoDepositEcalFromHits * process.stdEleIsoFromDepsEcal * process.modEleIsoFromDepsEcal *
                                          process.eleIsoDepositHcalDepth1FromTowers * process.stdEleIsoFromDepsHcalD1 * process.modEleIsoFromDepsHcalD1)



####################################################################################################
## Calling the NTupler , and defining sequence for running modules ...
process.demo = cms.EDAnalyzer('BstdZeeNTupler',
                              dyJetsToLL_EventType = cms.untracked.int32(0), #==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
                              isMC = cms.untracked.bool(True),
                              printOutInfo = cms.untracked.bool(False), 
                              readInNormReco = cms.untracked.bool(True),
                              readInBstdReco = cms.untracked.bool(False),
                              is2010SignalDataset = cms.untracked.bool(False),
                              useReducedRecHitsCollns = cms.untracked.bool(True),
                              vertexSrc = cms.untracked.InputTag("offlinePrimaryVertices"),
                              hltPathA_possNames = cms.untracked.vstring("HLT_DoubleEle33_CaloIdT_v2",
                                                                         "HLT_DoubleEle33_CaloIdT_v3",
                                                                         "HLT_DoubleEle33_CaloIdL_v1",
                                                                         "HLT_DoubleEle33_CaloIdL_v2",
                                                                         "HLT_DoubleEle33_CaloIdL_v3",
                                                                         "HLT_DoubleEle33_CaloIdL_v4",
                                                                         "HLT_DoubleEle33_CaloIdL_v5",# - run 177878 (3e33, v4.0) # (5e33, v1.4)
                                                                         "HLT_DoublePhoton33_v1",
                                                                         "HLT_DoublePhoton33_v2",
                                                                         "HLT_DoublePhoton33_v3"),
                              trg_emuPath_possNames = cms.untracked.vstring("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4",
                                                                            "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7",
                                                                            "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v1",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v2",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v3",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v4",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v5",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v6",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v7",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v8",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v9")
      ) 

#process.totalKinematicsFilter = cms.EDFilter('TotalKinematicsFilter',
#  src             = cms.InputTag("genParticles"),
#  tolerance       = cms.double(0.5),
#  verbose         = cms.untracked.bool(False)                                   
#)

process.p = cms.Path( process.heepIdNoIso * process.heepIdNoIsoEles * process.recalculateIsoVars * process.demo)
