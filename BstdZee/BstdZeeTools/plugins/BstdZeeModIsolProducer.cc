// -*- C++ -*-
//
// Package:    BstdZeeTools
// Class:      BstdZeeModIsolProducer
// 
/**\class BstdZeeModIsolProducer BstdZeeModIsolProducer.cc TSWilliams/BstdZeeTools/plugins/BstdZeeModIsolProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams
//         Created:  Mon Mar  5 18:12:27 GMT 2012
// $Id$
//
//

#include "TSWilliams/BstdZeeTools/interface/BstdZeeModIsolProducer.h"

// C++ includes
#include <memory>

// ROOT includes
#include <Math/VectorUtil.h>

// CMSSW includes -- Basic/General
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

//
// constructors and destructor
//
BstdZeeModIsolProducer::BstdZeeModIsolProducer(const edm::ParameterSet& iConfig) :
	inputGsfCollnTag_(iConfig.getParameter<edm::InputTag>("inputGsfEles")),
	vetoGsfCollnTag_(iConfig.getParameter<edm::InputTag>("vetoGsfEles")),
	extRadius_(0.3),
	// Track isolation
	ctfTracksTag_(iConfig.getParameter<edm::InputTag>("ctfTracksTag")),
	tk_intRadiusBarrel_(iConfig.getParameter<double>("intRadiusBarrelTk")),
	tk_intRadiusEndcap_(iConfig.getParameter<double>("intRadiusEndcapTk")),
	tk_stripWidthBarrel_(iConfig.getParameter<double>("stripBarrelTk")),
	tk_stripWidthEndcap_(iConfig.getParameter<double>("stripEndcapTk")),
	tk_ptMin_(iConfig.getParameter<double>("ptMinTk")),
	tk_maxVtxDist_(iConfig.getParameter<double>("maxVtxDistTk")),
	tk_beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
	tk_drbMax_(iConfig.getParameter<double>("maxDrbTk")),
	// ECAL isolation
	ecalParams_EB_(iConfig.getParameter<edm::InputTag>("barrelRecHitsTag"),
						iConfig.getParameter<double>("intRadiusEcalBarrel"),
						iConfig.getParameter<double>("jurassicWidth"),
						iConfig.getParameter<double>("etMinBarrel"),
						iConfig.getParameter<double>("eMinBarrel") ),
	ecalParams_EE_(iConfig.getParameter<edm::InputTag>("endcapRecHitsTag"),
						iConfig.getParameter<double>("intRadiusEcalEndcaps"),
						iConfig.getParameter<double>("jurassicWidth"),
						iConfig.getParameter<double>("etMinEndcaps"),
						iConfig.getParameter<double>("eMinEndcaps") ),
	ecal_vetoClustered_(iConfig.getParameter<bool>("vetoClustered")),
	ecal_useNumCrystals_(iConfig.getParameter<bool>("useNumCrystals")),
	ecal_severityLevelCut_(iConfig.getParameter<int>("severityLevelCut")),
	// HCAL depth 1 isolation
	hcalTowersTag_(iConfig.getParameter<edm::InputTag>("hcalTowers")),
	hcal_intRadius_(iConfig.getParameter<double>("intRadiusHcal")),
	hcal_etMin_(iConfig.getParameter<double>("etMinHcal"))
{
   // Now, register the products
	produces< edm::ValueMap<double> >("track");
	produces< edm::ValueMap<double> >("ecal");
	produces< edm::ValueMap<double> >("hcalDepth1");
}


BstdZeeModIsolProducer::~BstdZeeModIsolProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void BstdZeeModIsolProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// Grab the electron collections
	edm::Handle< reco::GsfElectronCollection > inputEleHandle, vetoEleHandle;
	iEvent.getByLabel(inputGsfCollnTag_, inputEleHandle);
	iEvent.getByLabel(vetoGsfCollnTag_, vetoEleHandle);

	// Prepare the output products
	std::auto_ptr< edm::ValueMap<double> > trackMap( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > ecalMap( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > hcalD1Map( new edm::ValueMap<double> );

//	// Run over the input electrons, calculating the new isolation values for each one
//	for( reco::GsfElectronCollection::const_iterator gsfEle = inputEleHandle->begin();
//			gsfEle!=inputEleHandle->end(); gsfEle++) {
	std::vector<double> trackIsolVec  = getTrackIsol(*inputEleHandle, *vetoEleHandle, iEvent);
	std::vector<double> ecalIsolVec   = getEcalIsol(*inputEleHandle, *vetoEleHandle, iEvent, iSetup);
	std::vector<double> hcalD1IsolVec = getHcalDepth1Isol(*inputEleHandle, *vetoEleHandle, iEvent);


	edm::ValueMap<double>::Filler trackMapFiller(*trackMap);
	trackMapFiller.insert(inputEleHandle, trackIsolVec.begin(), trackIsolVec.end());
	trackMapFiller.fill();
	edm::ValueMap<double>::Filler ecalMapFiller(*ecalMap);
	ecalMapFiller.insert(inputEleHandle, ecalIsolVec.begin(), ecalIsolVec.end());
	ecalMapFiller.fill();
	edm::ValueMap<double>::Filler hcalD1MapFiller(*hcalD1Map);
	hcalD1MapFiller.insert(inputEleHandle, hcalD1IsolVec.begin(), hcalD1IsolVec.end());
	hcalD1MapFiller.fill();

	// Store the products
	iEvent.put(trackMap, "track");
	iEvent.put(ecalMap, "ecal");
	iEvent.put(hcalD1Map, "hcalDepth1");
 
}

// ------------ method called once each job just before starting event loop  ------------
void BstdZeeModIsolProducer::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void BstdZeeModIsolProducer::endJob() { }

// ------------ method called when starting to processes a run  ------------
void BstdZeeModIsolProducer::beginRun(edm::Run&, edm::EventSetup const&)
{ }

// ------------ method called when ending the processing of a run  ------------
void BstdZeeModIsolProducer::endRun(edm::Run&, edm::EventSetup const&)
{ }

std::vector<double> BstdZeeModIsolProducer::getTrackIsol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent)
{
	// Grab the CTF track collection
	edm::Handle<reco::TrackCollection> ctfTracksHandle;
	iEvent.getByLabel(ctfTracksTag_,ctfTracksHandle);

	edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
	iEvent.getByLabel(tk_beamSpotTag_,recoBeamSpotHandle);

	std::vector<double> isolValues;
	// Run over the input electrons, calculating the new isolation values for each one
	for( reco::GsfElectronCollection::const_iterator gsfIt = inputEles.begin(); gsfIt!=inputEles.end(); gsfIt++)
		isolValues.push_back( getTrackIsol(*gsfIt, vetoEles, ctfTracksHandle.product(), recoBeamSpotHandle->position()) );

	return isolValues;
}

double BstdZeeModIsolProducer::getTrackIsol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const reco::TrackCollection* trackCollection, const reco::TrackBase::Point beamPoint)
{
	const int dzOption = egammaisolation::EgammaTrackSelector::vz;
	const double lip = tk_maxVtxDist_;

	double ptSum =0.;
	//Take the electron track
	reco::GsfTrackRef tmpTrack = theEle.gsfTrack() ;
	math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).momentum () ;
	double tmpElectronEtaAtVertex = (*tmpTrack).eta();

	for ( reco::TrackCollection::const_iterator itrTr  = (*trackCollection).begin() ;
			itrTr != (*trackCollection).end(); ++itrTr ) {
		math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).momentum () ;

		double this_pt  = (*itrTr).pt();
		if ( this_pt < tk_ptMin_ ) continue;

		double dzCut = 0;
		switch( dzOption ) {
			case egammaisolation::EgammaTrackSelector::dz : dzCut = fabs( (*itrTr).dz() - (*tmpTrack).dz() ); break;
			case egammaisolation::EgammaTrackSelector::vz : dzCut = fabs( (*itrTr).vz() - (*tmpTrack).vz() ); break;
			case egammaisolation::EgammaTrackSelector::bs : dzCut = fabs( (*itrTr).dz(beamPoint) - (*tmpTrack).dz(beamPoint) ); break;
			case egammaisolation::EgammaTrackSelector::vtx: dzCut = fabs( (*itrTr).dz(tmpTrack->vertex()) ); break;
			default : dzCut = fabs( (*itrTr).vz() - (*tmpTrack).vz() ); break;
		}
		if (dzCut > lip ) continue;
		if ( fabs( (*itrTr).dxy(beamPoint) ) > tk_drbMax_ ) continue;
		double dr = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(),tmpElectronMomentumAtVtx) ;
		double deta = (*itrTr).eta() - tmpElectronEtaAtVertex;

		bool includeInIsolSum = false;

		if (fabs(tmpElectronEtaAtVertex) < 1.479) {
			if ( fabs(dr) < extRadius_ && fabs(dr) >= tk_intRadiusBarrel_ && fabs(deta) >= tk_stripWidthBarrel_)
				includeInIsolSum=true;
		}
		else {
			if ( fabs(dr) < extRadius_ && fabs(dr) >= tk_intRadiusEndcap_ && fabs(deta) >= tk_stripWidthEndcap_)
				includeInIsolSum = true;
		}

		//TODO -- Check if the track is in the inner veto area of any of the vetoEles

		if(includeInIsolSum)
			ptSum += this_pt;

	  } //end loop over tracks

	  return ptSum;
}

std::vector<double> BstdZeeModIsolProducer::getEcalIsol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// Grab the recHit collections
	edm::Handle<EcalRecHitCollection> recHitsHandleEB, recHitsHandleEE;
	iEvent.getByLabel(ecalParams_EB_.recHitsTag, recHitsHandleEB);
	iEvent.getByLabel(ecalParams_EE_.recHitsTag, recHitsHandleEE);

	EcalRecHitMetaCollection* recHitsMetaEB = new EcalRecHitMetaCollection(*recHitsHandleEB);
	EcalRecHitMetaCollection* recHitsMetaEE = new EcalRecHitMetaCollection(*recHitsHandleEE);

	// Grab the ECAL geometry
	edm::ESHandle<CaloGeometry> theCaloGeom;
	iSetup.get<CaloGeometryRecord>().get(theCaloGeom);

	std::vector<double> isolValues;
	for( reco::GsfElectronCollection::const_iterator eleIt = inputEles.begin(); eleIt!=inputEles.end(); eleIt++)
		isolValues.push_back( getEcalIsol(*eleIt, vetoEles, recHitsMetaEB, theCaloGeom, ecalParams_EB_) + getEcalIsol(*eleIt, vetoEles, recHitsMetaEE, theCaloGeom, ecalParams_EE_) );

	delete recHitsMetaEB;
	delete recHitsMetaEE;

	return isolValues;
}


double BstdZeeModIsolProducer::getEcalIsol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const EcalRecHitMetaCollection* recHitsMeta, const edm::ESHandle<CaloGeometry>& theCaloGeom, const BstdZeeModIsolProducer::ECALParams& params)
{
	const CaloGeometry* caloGeom = theCaloGeom.product();
	const CaloSubdetectorGeometry* subdetGeoms[2];
	subdetGeoms[0] = theCaloGeom->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
	subdetGeoms[1] = theCaloGeom->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);

	double energySum = 0.;
	if (recHitsMeta){
		//Take the SC position
		reco::SuperClusterRef sc = theEle.get<reco::SuperClusterRef>();
		math::XYZPoint theCaloPosition = sc.get()->position();
		GlobalPoint pclu (theCaloPosition.x () ,
				theCaloPosition.y () ,
				theCaloPosition.z () );
		double etaclus = pclu.eta();
		double phiclus = pclu.phi();
		double r2 = params.intRadius*params.intRadius;

		std::vector< std::pair<DetId, float> >::const_iterator rhIt;

		for(int subdetnr=0; subdetnr<=1 ; subdetnr++){  // look in barrel and endcap
			CaloSubdetectorGeometry::DetIdSet chosen = subdetGeoms[subdetnr]->getCells(pclu,extRadius_);// select cells around cluster
			CaloRecHitMetaCollectionV::const_iterator j=recHitsMeta->end();
			for (CaloSubdetectorGeometry::DetIdSet::const_iterator  i = chosen.begin ();i!= chosen.end ();++i){//loop selected cells

				j=recHitsMeta->find(*i); // find selected cell among rechits
				if( j!=recHitsMeta->end()){ // add rechit only if available
					const  GlobalPoint & position = theCaloGeom.product()->getPosition(*i);
					double eta = position.eta();
					double phi = position.phi();
					double etaDiff = eta - etaclus;
					double phiDiff= reco::deltaPhi(phi,phiclus);
					double energy = j->energy();

					if(ecal_useNumCrystals_) {
						if( fabs(etaclus) < 1.479 ) {  // Barrel num crystals, crystal width = 0.0174
							if ( fabs(etaDiff) < 0.0174*params.etaSlice) continue;
							if ( sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.0174*params.intRadius) continue;
						} else {                       // Endcap num crystals, crystal width = 0.00864*fabs(sinh(eta))
							if ( fabs(etaDiff) < 0.00864*fabs(sinh(eta))*params.etaSlice) continue;
							if	 ( sqrt(etaDiff*etaDiff + phiDiff*phiDiff) < 0.00864*fabs(sinh(eta))*params.intRadius) continue;
						}
					} else {
						if ( fabs(etaDiff) < params.etaSlice) continue;  // jurassic strip cut
						if ( etaDiff*etaDiff + phiDiff*phiDiff < r2) continue; // jurassic exclusion cone cut
					}

					//Check if RecHit is in SC
					if(ecal_vetoClustered_) {
						//Loop over basic clusters:
						bool isClustered = false;
						for(reco::CaloCluster_iterator bcIt = sc->clustersBegin();bcIt != sc->clustersEnd(); ++bcIt) {
							for(rhIt = (*bcIt)->hitsAndFractions().begin();rhIt != (*bcIt)->hitsAndFractions().end(); ++rhIt) {
								if( rhIt->first == *i ) isClustered = true;
								if( isClustered ) break;
							}
							if( isClustered ) break;
						} //end loop over basic clusters
						if(isClustered) continue;
					}  //end if vetoClustered_

					//TODO -- Check if recHit is in inner veto area around any of the 'veto' eles

					//TODO -- Uncomment spike-removal & add relevant members/vars
//					//Severity level check
//					//make sure we have a barrel rechit
//					//call the severity level method
//					//passing the EBDetId
//					//the rechit collection in order to calculate the swiss crss
//					//and the EcalChannelRecHitRcd
//					//only consider rechits with ET >
//					//the SpikeId method (currently kE1OverE9 or kSwissCross)
//					//cut value for above
//					//then if the severity level is too high, we continue to the next rechit
//
//					if( ecal_severityLevelCut_!=-1 && ecalBarHits_ &&
//							sevLevel_->severityLevel(EBDetId(j->detid()), *ecalBarHits_) >= ecal_severityLevelCut_)
//						continue;
//					//                            *chStatus_,
//					//        severityRecHitThreshold_,
//					//        spId_,
//					//        spIdThreshold_
//					//    ) >= severityLevelCut_) continue;
//
//					//Check based on flags to protect from recovered channels from non-read towers
//					//Assumption is that v_chstatus_ is empty unless doFlagChecks() has been called
//					std::vector<int>::const_iterator vit = std::find( v_chstatus_.begin(), v_chstatus_.end(),  ((EcalRecHit*)(&*j))->recoFlag() );
//					if ( vit != v_chstatus_.end() ) continue; // the recHit has to be excluded from the iso sum

					double et = energy*position.perp()/position.mag();
					if ( fabs(et) > params.etMin && fabs(energy) > params.eMin ){ //Changed energy --> fabs(energy)
						energySum+=et;
					}

				} //End if not end of list
			} //End loop over rechits
		} //End loop over barrel/endcap
	} //End if recHitsMeta
	return energySum;
}

std::vector<double> BstdZeeModIsolProducer::getHcalDepth1Isol(const reco::GsfElectronCollection& inputEles, const reco::GsfElectronCollection& vetoEles, const edm::Event& iEvent)
{
	// Grab the CaloTowers collection
	edm::Handle<CaloTowerCollection> caloTowers;
	iEvent.getByLabel(hcalTowersTag_, caloTowers);

	// Run over the input electrons, calculating the new isolation values for each one
	std::vector<double> isolValues;
	for( reco::GsfElectronCollection::const_iterator gsfEle = inputEles.begin(); gsfEle!=inputEles.end(); gsfEle++)
		isolValues.push_back( getHcalDepth1Isol(*gsfEle, vetoEles, caloTowers.product()) );

	return isolValues;
}

double BstdZeeModIsolProducer::getHcalDepth1Isol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const CaloTowerCollection* towercollection)
{
	signed int depth = HcalDepth1;

	const reco::SuperCluster* sc = theEle.get<reco::SuperClusterRef>().get();
	math::XYZPoint theCaloPosition = sc->position();
	double candEta=sc->eta();
	double candPhi=sc->phi();
	double ptSum =0.;

	//loop over caloTowers
	for(CaloTowerCollection::const_iterator trItr = towercollection->begin(); trItr != towercollection->end(); ++trItr){

		double this_pt=0;
		switch(depth){
		case HcalAllDepths: this_pt = trItr->hadEt();break;
		case HcalDepth1: this_pt = trItr->ietaAbs()<18 || trItr->ietaAbs()>29 ? trItr->hadEt() : trItr->hadEnergyHeInnerLayer()*sin(trItr->p4().theta());break;
		case HcalDepth2: this_pt = trItr->hadEnergyHeOuterLayer()*sin(trItr->p4().theta());break;
		default: throw cms::Exception("Configuration Error") << "EgammaTowerIsolation: Depth " << depth << " not known. "; break;
		}

		if ( this_pt < hcal_etMin_ )
			continue ;

		double towerEta=trItr->eta();
		double towerPhi=trItr->phi();
		double twoPi= 2*M_PI;
		if(towerPhi<0) towerPhi+=twoPi;
		if(candPhi<0) candPhi+=twoPi;
		double deltaPhi=fabs(towerPhi-candPhi);
		if(deltaPhi>twoPi) deltaPhi-=twoPi;
		if(deltaPhi>M_PI) deltaPhi=twoPi-deltaPhi;
		double deltaEta = towerEta - candEta;

		// TODO -- Check if caloTower is within inner veto area of any of the other eles

		double dr = deltaEta*deltaEta + deltaPhi*deltaPhi;
		if( dr < extRadius_*extRadius_ &&
				dr >= hcal_intRadius_*hcal_intRadius_ )
		{
			ptSum += this_pt;
		}
	}//end loop over caloTowers

	return ptSum;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BstdZeeModIsolProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BstdZeeModIsolProducer);
