import FWCore.ParameterSet.Config as cms

#####################################
# Input flags
input_isMC = False
input_is2010SignalMC = False
input_dyJetsToLLFilter = 0 #==0 =>Don't select events, ==11 =>ele, ==13 =>muon, ==15 =>tau

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

######################################################
## For reading in options from command line arguments
import FWCore.ParameterSet.VarParsing as VarParsing
import os

options = VarParsing.VarParsing ('analysis')

# setup any defaults you want
options.outputFile = 'tmp.root'
options.inputFiles = 'file1.root', 'file2.root'
options.maxEvents  = -1 # -1 means all events

# get and parse the command line arguments
options.parseArguments()

## Act appropriately if 1st specified inputFile is a directory
if os.path.isdir(options.inputFiles[0]) or os.path.isdir('/pnfs/pp.rl.ac.uk/data/cms/'+options.inputFiles[0]):
   # Grab directory location from cmd line arg, and list files within
   dirLocation = options.inputFiles.pop(0).replace('/pnfs/pp.rl.ac.uk/data/cms/store/','/store/')

   if(dirLocation).startswith('/store/'):
      options.inputFiles = os.listdir('/pnfs/pp.rl.ac.uk/data/cms'+dirLocation)
   else:
      options.inputFiles = os.listdir(dirLocation)

   # Pre-pend the directory to the file names 
   for i in range(len(options.inputFiles)):
      options.inputFiles[i] = dirLocation+"/"+options.inputFiles[i]      
else:
   for i in range(len(options.inputFiles)):
      options.inputFiles[i] = DataFileLocationAdaptor(options.inputFiles[i])

#####################################################
## General setup/load lines, defining nevts, input & output etc.

process = cms.Process("NTupler")

## Loading Geometry modules for creation of transient geometry information used in determining eta and phi values of recHits ...
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START42_V13::All'
process.GlobalTag.globaltag = 'START52_V9::All'

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

process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

#Reading in a list of datafile locations from a file...
fileName_DatafileList = "/home/ppd/nnd85574/FileLocations/BstdZee/Boosted_Quark_To_ee_2TeV/tsw-bstd412Recon_2-00TeVu_v2-c17b3790b1be6b1dd249df112816fe9d/USER/fileList.txt"
f_DatafileList = open(fileName_DatafileList)
datafilesList = f_DatafileList.readlines()
f_DatafileList.close()
datafileLocations = map(DataFileLocationAdaptor,datafilesList)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFiles) )

#Defining the output file to store the histograms/NTuples in... 
process.TFileService = cms.Service("TFileService", fileName=cms.string("BstdZeeNTuple_53X-v2pre3.root"))
#process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))

###################################################################################################
## Code for modified isolation values ...
# ModIso: 1) Creating a collection of HEEP ID (no iso) passing electrons ...
#setting up the producer to make the HEEP ID value map
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                eleLabel = cms.InputTag("gsfElectrons"),
                                barrelCuts = cms.PSet(heepBarrelCuts),
                                endcapCuts = cms.PSet(heepEndcapCuts),
                                verticesLabel = cms.InputTag("offlinePrimaryVertices"),
                                eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
                                eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                applyRhoCorrToEleIsol = cms.bool(True),
                                writeIdAsInt = cms.bool(True)
                                )
process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")  )

# ModIso: 2a) Calculating the modified iso. values using BstdZeeTools EDProducer

from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *

process.innerXSVetoModEleIso = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
process.innerXSVetoModEleIso.otherElesIntRadiusBarrelTk    = cms.double(0.015)
process.innerXSVetoModEleIso.otherElesIntRadiusEndcapTk    = cms.double(0.015)
process.innerXSVetoModEleIso.otherElesStripBarrelTk        = cms.double(0.015)
process.innerXSVetoModEleIso.otherElesStripEndcapTk        = cms.double(0.015)
process.innerXSVetoModEleIso.otherElesIntRadiusEcalBarrel  = cms.double(3.0)
process.innerXSVetoModEleIso.otherElesIntRadiusEcalEndcaps = cms.double(3.0)
process.innerXSVetoModEleIso.otherElesJurassicWidth        = cms.double(1.5)


process.innerSVetoModEleIso = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
process.innerSVetoModEleIso.otherElesIntRadiusBarrelTk    = cms.double(0.03)
process.innerSVetoModEleIso.otherElesIntRadiusEndcapTk    = cms.double(0.03)
process.innerSVetoModEleIso.otherElesStripBarrelTk        = cms.double(0.03)
process.innerSVetoModEleIso.otherElesStripEndcapTk        = cms.double(0.03)
process.innerSVetoModEleIso.otherElesIntRadiusEcalBarrel  = cms.double(3.5)
process.innerSVetoModEleIso.otherElesIntRadiusEcalEndcaps = cms.double(3.5)
process.innerSVetoModEleIso.otherElesJurassicWidth        = cms.double(1.75)


process.innerMVetoModEleIso = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
process.innerMVetoModEleIso.otherElesIntRadiusBarrelTk    = cms.double(0.045)
process.innerMVetoModEleIso.otherElesIntRadiusEndcapTk    = cms.double(0.045)
process.innerMVetoModEleIso.otherElesStripBarrelTk        = cms.double(0.045)
process.innerMVetoModEleIso.otherElesStripEndcapTk        = cms.double(0.045)
process.innerMVetoModEleIso.otherElesIntRadiusEcalBarrel  = cms.double(4.0)
process.innerMVetoModEleIso.otherElesIntRadiusEcalEndcaps = cms.double(4.0)
process.innerMVetoModEleIso.otherElesJurassicWidth        = cms.double(2.0)


process.innerLVetoModEleIso = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
process.innerLVetoModEleIso.otherElesIntRadiusBarrelTk    = cms.double(0.06)
process.innerLVetoModEleIso.otherElesIntRadiusEndcapTk    = cms.double(0.06)
process.innerLVetoModEleIso.otherElesStripBarrelTk        = cms.double(0.06)
process.innerLVetoModEleIso.otherElesStripEndcapTk        = cms.double(0.06)
process.innerLVetoModEleIso.otherElesIntRadiusEcalBarrel  = cms.double(5.0)
process.innerLVetoModEleIso.otherElesIntRadiusEcalEndcaps = cms.double(5.0)
process.innerLVetoModEleIso.otherElesJurassicWidth        = cms.double(2.5)


process.innerXLVetoModEleIso = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
process.innerXLVetoModEleIso.otherElesIntRadiusBarrelTk    = cms.double(0.09)
process.innerXLVetoModEleIso.otherElesIntRadiusEndcapTk    = cms.double(0.09)
process.innerXLVetoModEleIso.otherElesStripBarrelTk        = cms.double(0.09)
process.innerXLVetoModEleIso.otherElesStripEndcapTk        = cms.double(0.09)
process.innerXLVetoModEleIso.otherElesIntRadiusEcalBarrel  = cms.double(7.0)
process.innerXLVetoModEleIso.otherElesIntRadiusEcalEndcaps = cms.double(7.0)
process.innerXLVetoModEleIso.otherElesJurassicWidth        = cms.double(3.5)


process.vetoAreaModEleIsoPhantomDr005To010 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"),
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.05),
      phantomVetoEleDrMax = cms.double(0.10) )

process.vetoAreaModEleIsoPhantomDr010To015 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.10),
      phantomVetoEleDrMax = cms.double(0.15) )

process.vetoAreaModEleIsoPhantomDr015To020 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.15),
      phantomVetoEleDrMax = cms.double(0.20) )

process.vetoAreaModEleIsoPhantomDr020To025 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.20),
      phantomVetoEleDrMax = cms.double(0.25) )

process.vetoAreaModEleIsoPhantomDr025To030 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.25),
      phantomVetoEleDrMax = cms.double(0.30) )

process.vetoAreaModEleIsoPhantomDr030To035 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.30),
      phantomVetoEleDrMax = cms.double(0.35) )

process.vetoAreaModEleIsoPhantomDr035To040 = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles"), 
      extraPhantomVetoEle = cms.bool(True),
      phantomVetoEleDrMin = cms.double(0.35),
      phantomVetoEleDrMax = cms.double(0.40) )
      
process.modEleIsoOtherEleVetoAreaForSelf = cms.EDProducer("BstdZeeModIsolProducer",
      bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles")  )
process.modEleIsoOtherEleVetoAreaForSelf.intRadiusBarrelTk    = bstdZeeModIsolParams.otherElesIntRadiusBarrelTk
process.modEleIsoOtherEleVetoAreaForSelf.intRadiusEndcapTk    = bstdZeeModIsolParams.otherElesIntRadiusEndcapTk
process.modEleIsoOtherEleVetoAreaForSelf.stripBarrelTk        = bstdZeeModIsolParams.otherElesStripBarrelTk
process.modEleIsoOtherEleVetoAreaForSelf.stripEndcapTk        = bstdZeeModIsolParams.otherElesStripEndcapTk
process.modEleIsoOtherEleVetoAreaForSelf.intRadiusEcalBarrel  = bstdZeeModIsolParams.otherElesIntRadiusEcalBarrel
process.modEleIsoOtherEleVetoAreaForSelf.intRadiusEcalEndcaps = bstdZeeModIsolParams.otherElesIntRadiusEcalEndcaps
process.modEleIsoOtherEleVetoAreaForSelf.jurassicWidth        = bstdZeeModIsolParams.otherElesJurassicWidth


#process.innerXSVetoModEleIso.barrelRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.innerXSVetoModEleIso.endcapRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")
#process.innerSVetoModEleIso.barrelRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.innerSVetoModEleIso.endcapRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")
#process.innerMVetoModEleIso.barrelRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.innerMVetoModEleIso.endcapRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")
#process.innerLVetoModEleIso.barrelRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.innerLVetoModEleIso.endcapRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")
#process.innerXLVetoModEleIso.barrelRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.innerXLVetoModEleIso.endcapRecHitsTag = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")

process.innerVetoModEleIsos = cms.Sequence( process.innerXSVetoModEleIso * process.innerSVetoModEleIso * process.innerMVetoModEleIso 
                                            * process.innerLVetoModEleIso * process.innerXLVetoModEleIso 
                                            * process.vetoAreaModEleIsoPhantomDr005To010 * process.vetoAreaModEleIsoPhantomDr010To015 
                                            * process.vetoAreaModEleIsoPhantomDr015To020 * process.vetoAreaModEleIsoPhantomDr020To025
                                            * process.vetoAreaModEleIsoPhantomDr025To030 * process.vetoAreaModEleIsoPhantomDr030To035
                                            * process.vetoAreaModEleIsoPhantomDr035To040
                                            * process.modEleIsoOtherEleVetoAreaForSelf )


# ModIso: 2b) Calculating modified iso. values with IsoDeposits framework

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
#process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("ecalRecHit:EcalRecHitsEB:stdRECO")
#process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("ecalRecHit:EcalRecHitsEE:stdRECO")
process.eleIsoDepositEcalFromHits.ExtractorPSet.barrelEcalHits = cms.InputTag("reducedEcalRecHitsEB::RECO")
process.eleIsoDepositEcalFromHits.ExtractorPSet.endcapEcalHits = cms.InputTag("reducedEcalRecHitsEE::RECO")
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
                           'NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThreshold(0.08)',
                           'EcalEndcaps:AbsThresholdFromTransverse(0.1)'),
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
                           'NumCrystalEtaPhiVeto(1.5,9999.0)',
                           'EcalBarrel:AbsThreshold(0.08)',
                           'EcalEndcaps:AbsThresholdFromTransverse(0.1)',
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
process.recalcIsoVarsWithIsoDeps = cms.Sequence(process.eleIsoDepositTk * process.stdEleIsoFromDepsTk * process.modEleIsoFromDepsTk *
                                          process.eleIsoDepositEcalFromHits * process.stdEleIsoFromDepsEcal * process.modEleIsoFromDepsEcal *
                                          process.eleIsoDepositHcalDepth1FromTowers * process.stdEleIsoFromDepsHcalD1 * process.modEleIsoFromDepsHcalD1)

####################################################################################################
## Code for calculating the 'rho' value for PU correction ...

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5) 


####################################################################################################
## Calling the NTupler , and defining sequence for running modules ...
process.demo = cms.EDAnalyzer('BstdZeeNTupler',
                              dyJetsToLL_EventType = cms.untracked.int32(input_dyJetsToLLFilter), #==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
                              isMC = cms.untracked.bool(input_isMC),
                              printOutInfo = cms.untracked.bool(False), 
                              readInNormReco = cms.untracked.bool(True), readInBstdReco = cms.untracked.bool(False), readInTrigInfo = cms.untracked.bool(True),
                              eleCollection = cms.untracked.InputTag("gsfElectrons"),
                              useReducedRecHitsCollns = cms.untracked.bool(True) ,
                              is2010SignalDataset = cms.untracked.bool(input_is2010SignalMC),
                              vertexSrc = cms.untracked.InputTag("offlinePrimaryVertices"),
                              hltPathA_possNames = cms.untracked.vstring("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v\\\\d+",     # 2012 trigger
                                                                         "HLT_DoubleEle33_CaloIdT_v[2-3]",
                                                                         "HLT_DoubleEle33_CaloIdL_v[1-5]",
                                                                         "HLT_DoublePhoton33_v[1-3]"),
                              trg_emuPath_possNames = cms.untracked.vstring("HLT_Mu22_Photon22_CaloIdL_v\\\\d+",     # 2012 trigger
                                                                            "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v[478]",
                                                                            "HLT_Mu15_Photon20_CaloIdL_v[1-9]"),
                              puDists_mcFile       = cms.untracked.string("Summer12DR53XPileUp_true_20121112.root"),
                              puDists_dataFile     = cms.untracked.string("data12PileUp_true_20121127_r190456-206940.root"),
                              puDists_mcHistName   = cms.untracked.string("Summer12DR53XPileUpHist_true"),
                              puDists_dataHistName = cms.untracked.string("pileup")
     ) 


##################################################################################
# Useful modules for looking at details of generator event record 
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleListDrawer",
#  maxEventsToPrint = cms.untracked.int32(100),
#  printVertex = cms.untracked.bool(False),
#  src = cms.InputTag("genParticles")
#)
#
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printDecayTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#                                   printP4 = cms.untracked.bool(False),
#                                   printPtEtaPhi = cms.untracked.bool(True),
#                                   printVertex = cms.untracked.bool(False),
#                                   printStatus = cms.untracked.bool(True),
#                                   printIndex = cms.untracked.bool(True)
#                                   )

##################################################################################
# Construct the whole path ...
process.p = cms.Path( process.heepIdNoIso * process.heepIdNoIsoEles * process.innerVetoModEleIsos * process.recalcIsoVarsWithIsoDeps * process.kt6PFJets * process.demo )

## For checking transient event content ....
#process.Out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string (options.outputFile+"_AODPlus.root") )
#process.end = cms.EndPath( process.Out )
