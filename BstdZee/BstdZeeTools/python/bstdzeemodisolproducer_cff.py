import FWCore.ParameterSet.Config as cms

import RecoEgamma.EgammaElectronProducers.gsfElectrons_cfi as stdEleRecon

stdGsfEles = stdEleRecon.ecalDrivenGsfElectrons

bstdZeeModIsolParams = cms.PSet (
   inputGsfEles = cms.InputTag("gsfElectrons"),
   
   ctfTracksTag = cms.InputTag("generalTracks"),
   intRadiusBarrelTk = stdGsfEles.intRadiusBarrelTk,
   intRadiusEndcapTk = stdGsfEles.intRadiusEndcapTk,
   stripBarrelTk     = stdGsfEles.stripBarrelTk,
   stripEndcapTk     = stdGsfEles.stripEndcapTk,
   ptMinTk           = stdGsfEles.ptMinTk,
   maxVtxDistTk      = stdGsfEles.maxVtxDistTk,
   beamSpotTag       = stdGsfEles.beamSpotTag,
   maxDrbTk          = stdGsfEles.maxDrbTk,
   
   barrelRecHitsTag     = cms.InputTag("reducedEcalRecHitsEB::RECO"),
   endcapRecHitsTag     = cms.InputTag("reducedEcalRecHitsEE::RECO"),
   intRadiusEcalBarrel  = stdGsfEles.intRadiusEcalBarrel,
   intRadiusEcalEndcaps = stdGsfEles.intRadiusEcalEndcaps,
   jurassicWidth        = stdGsfEles.jurassicWidth,
   etMinBarrel          = stdGsfEles.etMinBarrel,
   eMinBarrel           = stdGsfEles.eMinBarrel,
   etMinEndcaps         = stdGsfEles.etMinEndcaps,
   eMinEndcaps          = stdGsfEles.eMinEndcaps,
   vetoClustered        = stdGsfEles.vetoClustered,
   useNumCrystals       = stdGsfEles.useNumCrystals,
   severityLevelCut     = stdGsfEles.severityLevelCut,
   
   hcalTowers    = cms.InputTag("towerMaker"),
   intRadiusHcal = stdGsfEles.intRadiusHcal,
   etMinHcal     = stdGsfEles.etMinHcal 
)
