
#ifndef TSWilliams_BstdZeeTools_BstdZeeModIsolProducer_h
#define TSWilliams_BstdZeeTools_BstdZeeModIsolProducer_h

// General C++ include files
#include <vector>

// CMSSW includes -- Basic
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"


class BstdZeeModIsolProducer : public edm::EDProducer {
	public:
		enum HcalDepths{
			HcalAllDepths=-1,
			HcalUndefined=0,
			HcalDepth1=1,
			HcalDepth2=2
		};
		struct ECALParams{
			ECALParams(const edm::InputTag& tag, double radius, double width, double otherElesRadius, double otherElesWidth, double minEt, double minE) :
				recHitsTag(tag), intRadius(radius), etaSlice(width),
				otherElesIntRadius(otherElesRadius), otherElesEtaSlice(otherElesWidth),
				etMin(minEt), eMin(minE) { }
			edm::InputTag recHitsTag;
			double intRadius;
			double etaSlice;
			double otherElesIntRadius;
			double otherElesEtaSlice;
			double etMin;
			double eMin;
		};
   public:
      explicit BstdZeeModIsolProducer(const edm::ParameterSet&);
      ~BstdZeeModIsolProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob();
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();

      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);

      // Methods for calculating isolation values
      std::vector<double> getTrackIsol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent);
      double getTrackIsol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const reco::TrackCollection* trackCollection, const reco::TrackBase::Point beamPoint);
      std::vector<double> getEcalIsol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent, const edm::EventSetup& iSetup);
      double getEcalIsol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const EcalRecHitMetaCollection* recHitsMeta, const edm::ESHandle<CaloGeometry>& theCaloGeom, const ECALParams& params, const EcalSeverityLevelAlgo* sevLevelAlgo=0);
      std::vector<double> getHcalDepth1Isol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent);
      double getHcalDepth1Isol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const CaloTowerCollection* towercollection);

      // Members grabbed from python config
      const edm::InputTag inputGsfCollnTag_;
      const edm::InputTag vetoGsfCollnTag_;

      const double extRadius_;

      const edm::InputTag ctfTracksTag_;
      const double tk_intRadiusBarrel_;
      const double tk_intRadiusEndcap_;
      const double tk_stripWidthBarrel_;
      const double tk_stripWidthEndcap_;
      const double tk_otherElesIntRadiusBarrel_;
      const double tk_otherElesIntRadiusEndcap_;
      const double tk_otherElesStripWidthBarrel_;
      const double tk_otherElesStripWidthEndcap_;
      const double tk_ptMin_;
      const double tk_maxVtxDist_;
      const edm::InputTag tk_beamSpotTag_;
      const double tk_drbMax_;

      std::vector<reco::Track> tracksInIsolSum_;

      const ECALParams ecalParams_EB_;
      const ECALParams ecalParams_EE_;
      const EcalRecHitCollection* ecalBarHits_;
      const bool ecal_vetoClustered_;
      const bool ecal_useNumCrystals_;
      const int ecal_severityLevelCut_;
      std::vector<int> recHitFlagsToBeExcluded_;

      const edm::InputTag hcalTowersTag_;
      const double hcal_intRadius_;
      const double hcal_etMin_;
};

#endif
