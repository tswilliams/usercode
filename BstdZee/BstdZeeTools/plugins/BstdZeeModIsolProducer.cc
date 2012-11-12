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
// $Id: BstdZeeModIsolProducer.cc,v 1.5 2012/10/07 23:55:33 tsw Exp $
//
//

#include "TSWilliams/BstdZeeTools/interface/BstdZeeModIsolProducer.h"

// C++ includes
#include <memory>

// ROOT includes
#include "TMath.h"
#include <Math/VectorUtil.h>
#include "TRandom3.h"

// CMSSW includes -- Basic/General
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTrackSelector.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
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
   tk_otherElesIntRadiusBarrel_(iConfig.getParameter<double>("otherElesIntRadiusBarrelTk")),
   tk_otherElesIntRadiusEndcap_(iConfig.getParameter<double>("otherElesIntRadiusEndcapTk")),
   tk_otherElesStripWidthBarrel_(iConfig.getParameter<double>("otherElesStripBarrelTk")),
   tk_otherElesStripWidthEndcap_(iConfig.getParameter<double>("otherElesStripEndcapTk")),
	tk_ptMin_(iConfig.getParameter<double>("ptMinTk")),
	tk_maxVtxDist_(iConfig.getParameter<double>("maxVtxDistTk")),
	tk_beamSpotTag_(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
	tk_drbMax_(iConfig.getParameter<double>("maxDrbTk")),
	// ECAL isolation
	ecalParams_EB_(iConfig.getParameter<edm::InputTag>("barrelRecHitsTag"),
						iConfig.getParameter<double>("intRadiusEcalBarrel"),
						iConfig.getParameter<double>("jurassicWidth"),
						iConfig.getParameter<double>("otherElesIntRadiusEcalBarrel"),
						iConfig.getParameter<double>("otherElesJurassicWidth"),
						iConfig.getParameter<double>("etMinBarrel"),
						iConfig.getParameter<double>("eMinBarrel"),
						iConfig.getParameter< std::vector<std::string> >("recHitFlagsToBeExcludedBarrel"),
						iConfig.getParameter< std::vector<std::string> >("recHitSeverityToBeExcludedBarrel") ),
	ecalParams_EE_(iConfig.getParameter<edm::InputTag>("endcapRecHitsTag"),
						iConfig.getParameter<double>("intRadiusEcalEndcaps"),
						iConfig.getParameter<double>("jurassicWidth"),
						iConfig.getParameter<double>("otherElesIntRadiusEcalEndcaps"),
						iConfig.getParameter<double>("otherElesJurassicWidth"),
						iConfig.getParameter<double>("etMinEndcaps"),
						iConfig.getParameter<double>("eMinEndcaps"),
						iConfig.getParameter< std::vector<std::string> >("recHitFlagsToBeExcludedEndcaps"),
						iConfig.getParameter< std::vector<std::string> >("recHitSeverityToBeExcludedEndcaps") ),
	ecal_vetoClustered_(iConfig.getParameter<bool>("vetoClustered")),
	ecal_useNumCrystals_(iConfig.getParameter<bool>("useNumCrystals")),
	//ecal_severityLevelCut_(iConfig.getParameter<int>("severityLevelCut")),
//   recHitFlagsToBeExcluded_(iConfig.getParameter< std::vector<int> >("recHitFlagsToBeExcluded")),
	// HCAL depth 1 isolation
	hcalTowersTag_(iConfig.getParameter<edm::InputTag>("hcalTowers")),
	hcal_intRadius_(iConfig.getParameter<double>("intRadiusHcal")),
	hcal_etMin_(iConfig.getParameter<double>("etMinHcal")),
	// Phantom ele params
   extraPhantomVetoEle_( iConfig.exists("extraPhantomVetoEle") ? iConfig.getParameter<bool>("extraPhantomVetoEle") : false ),
   phantomVetoEleDrMin_( iConfig.exists("phantomVetoEleDrMin") ? iConfig.getParameter<double>("phantomVetoEleDrMin") : 0.0 ),
   phantomVetoEleDrMax_( iConfig.exists("phantomVetoEleDrMax") ? iConfig.getParameter<double>("phantomVetoEleDrMax") : 0.0 ),
   phantomVetoEleDEta_(0.0), phantomVetoEleDPhi_(0.0)

{
   // Now, register the products
	produces< std::vector<reco::Track> >("tracksInIsolSum");
	produces< edm::ValueMap<double> >("track");
	produces< edm::ValueMap<double> >("ecal");
	produces< edm::ValueMap<double> >("hcalDepth1");
	if( extraPhantomVetoEle_ ){
		produces< edm::ValueMap<double> >("dEtaPhantomEle");
		produces< edm::ValueMap<double> >("dPhiPhantomEle");
	}
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
	tracksInIsolSum_.clear();
	// Grab the electron collections
	edm::Handle< reco::GsfElectronCollection > inputEleHandle, vetoEleHandle;
	iEvent.getByLabel(inputGsfCollnTag_, inputEleHandle);
	iEvent.getByLabel(vetoGsfCollnTag_, vetoEleHandle);

	// Prepare the output products
	std::auto_ptr< edm::ValueMap<double> > trackMap( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > ecalMap( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > hcalD1Map( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > dEtaPhantomEleMap( new edm::ValueMap<double> );
	std::auto_ptr< edm::ValueMap<double> > dPhiPhantomEleMap( new edm::ValueMap<double> );

	// If phantom veto ele requested, then set up it's dEta & dPhi w.r.t. ele whose isolation value is being calculated.
	TRandom3 randNumGen( iEvent.id().event() );
	const double phantomVetoEleDR = phantomVetoEleDrMin_ + (phantomVetoEleDrMax_-phantomVetoEleDrMin_)*randNumGen.Rndm();
	const double phantomVetoEleEtaPhiAngle = 2.0*TMath::Pi()*randNumGen.Rndm();
	phantomVetoEleDEta_ = phantomVetoEleDR*cos(phantomVetoEleEtaPhiAngle);
	phantomVetoEleDPhi_ = phantomVetoEleDR*sin(phantomVetoEleEtaPhiAngle);

	std::vector<double> dEtaPhantomEleVec(inputEleHandle->size(), phantomVetoEleDEta_);
	std::vector<double> dPhiPhantomEleVec(inputEleHandle->size(), phantomVetoEleDPhi_);

//	// Run over the input electrons, calculating the new isolation values for each one
//	for( reco::GsfElectronCollection::const_iterator gsfEle = inputEleHandle->begin();
//			gsfEle!=inputEleHandle->end(); gsfEle++) {
	edm::Handle<EcalRecHitCollection> recHitsHandleEB;
	iEvent.getByLabel(ecalParams_EB_.recHitsTag, recHitsHandleEB);
	ecalBarHits_ = recHitsHandleEB.product();

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

	edm::ValueMap<double>::Filler dEtaPhantomEleMapFiller(*dEtaPhantomEleMap);
	dEtaPhantomEleMapFiller.insert(inputEleHandle, dEtaPhantomEleVec.begin(), dEtaPhantomEleVec.end());
	dEtaPhantomEleMapFiller.fill();
	edm::ValueMap<double>::Filler dPhiPhantomEleMapFiller(*dPhiPhantomEleMap);
	dPhiPhantomEleMapFiller.insert(inputEleHandle, dPhiPhantomEleVec.begin(), dPhiPhantomEleVec.end());
	dPhiPhantomEleMapFiller.fill();

	// Store the products
	std::auto_ptr< std::vector<reco::Track> > tracksToStore( new std::vector<reco::Track>(tracksInIsolSum_) );
	iEvent.put(tracksToStore, "tracksInIsolSum");
	iEvent.put(trackMap, "track");
	iEvent.put(ecalMap, "ecal");
	iEvent.put(hcalD1Map, "hcalDepth1");
	if( extraPhantomVetoEle_ ){
		iEvent.put(dEtaPhantomEleMap, "dEtaPhantomEle");
		iEvent.put(dPhiPhantomEleMap, "dPhiPhantomEle");
	}
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
		//math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).momentum () ;

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
		if(!includeInIsolSum)
			continue;

		// *** Check if the track is in the inner veto area of any of the vetoEles
		//   If extraPhantomVetoEle_ true, then instead just check if track is w/i inner veto area of this phantom ele
		for( reco::GsfElectronCollection::const_iterator otherEleIt = vetoEles.begin(); otherEleIt!=vetoEles.end(); otherEleIt++){
			reco::GsfTrackRef otherEleTrack = otherEleIt->gsfTrack();

			// Skip if this is the electron who's isol value is being calculated
			if( !extraPhantomVetoEle_ && fabs(otherEleTrack->eta()-tmpTrack->eta())<0.001 && fabs(otherEleTrack->phi()-tmpTrack->phi())<0.001 )
				continue;

			math::XYZVector otherEleTrackP3 = otherEleTrack->momentum();
			if(extraPhantomVetoEle_){
				math::RhoEtaPhiVector rhoEtaPhiTrackP3;
				rhoEtaPhiTrackP3.SetEta(tmpTrack->eta()+phantomVetoEleDEta_).SetPhi(tmpTrack->phi()+phantomVetoEleDPhi_).SetRho(1.0);
				otherEleTrackP3.SetXYZ(rhoEtaPhiTrackP3.x(), rhoEtaPhiTrackP3.y(), rhoEtaPhiTrackP3.z());
			}
			double vsOtherEle_dR = ROOT::Math::VectorUtil::DeltaR(itrTr->momentum(), otherEleTrackP3);
			double vsOtherEle_dEta = (*itrTr).eta() - (extraPhantomVetoEle_ ? (tmpTrack->eta()+phantomVetoEleDEta_) : otherEleTrack->eta() );
			if (fabs(tmpElectronEtaAtVertex) < 1.479) {
				if ( fabs(vsOtherEle_dR) < tk_otherElesIntRadiusBarrel_ || ( fabs(vsOtherEle_dR) < extRadius_ && fabs(vsOtherEle_dEta) < tk_otherElesStripWidthBarrel_ ) )
					includeInIsolSum = false;
			}
			else {
				if ( fabs(vsOtherEle_dR) < tk_otherElesIntRadiusEndcap_ || ( fabs(vsOtherEle_dR) < extRadius_ && fabs(vsOtherEle_dEta) < tk_otherElesStripWidthEndcap_ ) )
					includeInIsolSum = false;
			}
			if(!includeInIsolSum)
				break;
		} // End loop over vetoEles

		if(includeInIsolSum){
			ptSum += this_pt;
			tracksInIsolSum_.push_back(*itrTr);
		}

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

	// Grab the ECAL severity level algo
	edm::ESHandle<EcalSeverityLevelAlgo> sevLevelAlgo;
   iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevLevelAlgo);

	std::vector<double> isolValues;
	for( reco::GsfElectronCollection::const_iterator eleIt = inputEles.begin(); eleIt!=inputEles.end(); eleIt++)
		isolValues.push_back( getEcalIsol(*eleIt, vetoEles, recHitsMetaEB, theCaloGeom, ecalParams_EB_, sevLevelAlgo.product()) + getEcalIsol(*eleIt, vetoEles, recHitsMetaEE, theCaloGeom, ecalParams_EE_, sevLevelAlgo.product()) );

	delete recHitsMetaEB;
	delete recHitsMetaEE;

	return isolValues;
}


double BstdZeeModIsolProducer::getEcalIsol(const reco::GsfElectron& theEle, const reco::GsfElectronCollection& vetoEles, const EcalRecHitMetaCollection* recHitsMeta, const edm::ESHandle<CaloGeometry>& theCaloGeom, const BstdZeeModIsolProducer::ECALParams& params, const EcalSeverityLevelAlgo* sevLevelAlgo)
{
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

					// *** Check if recHit is in inner veto area around any of the 'veto' eles
					//   If extraPhantomVetoEle_ true, then instead just check if recHit is w/i inner veto area of this phantom ele
					bool vetoedByVetoEles = false;
					for( reco::GsfElectronCollection::const_iterator otherEleIt = vetoEles.begin(); otherEleIt!=vetoEles.end(); otherEleIt++){
						reco::SuperClusterRef otherEleSC = otherEleIt->get<reco::SuperClusterRef>();
						math::XYZPoint otherEle_caloPosition = otherEleSC.get()->position();
						GlobalPoint otherEle_pclu (otherEle_caloPosition.x () ,
								otherEle_caloPosition.y () ,
								otherEle_caloPosition.z () );
						double otherEle_etaclus = extraPhantomVetoEle_ ? (etaclus+phantomVetoEleDEta_) : otherEle_pclu.eta();
						double otherEle_phiclus = extraPhantomVetoEle_ ? (phiclus+phantomVetoEleDPhi_) : otherEle_pclu.phi();

						// Skip if this is the electron who's isol value is being calculated
						if( fabs(otherEle_etaclus-etaclus)<0.001 && fabs(otherEle_phiclus-phiclus)<0.001 )
							continue;

						double otherEle_etaDiff = eta - otherEle_etaclus;
						double otherEle_phiDiff = reco::deltaPhi(phi, otherEle_phiclus);
						double otherEle_dR = sqrt(otherEle_etaDiff*otherEle_etaDiff + otherEle_phiDiff*otherEle_phiDiff );

						if(otherEle_dR>extRadius_)
							continue;
						if(ecal_useNumCrystals_) {
							if( fabs(etaclus) < 1.479 ) {  // Barrel num crystals, crystal width = 0.0174
								if ( fabs(otherEle_etaDiff) < 0.0174*params.otherElesEtaSlice ) vetoedByVetoEles=true;
								if ( otherEle_dR < 0.0174*params.otherElesIntRadius) vetoedByVetoEles=true;
							} else {                       // Endcap num crystals, crystal width = 0.00864*fabs(sinh(eta))
								if ( fabs(otherEle_etaDiff) < 0.00864*fabs(sinh(eta))*params.otherElesEtaSlice) vetoedByVetoEles=true;
								if	 ( otherEle_dR < 0.00864*fabs(sinh(eta))*params.otherElesIntRadius) vetoedByVetoEles=true;
							}
						} else {
							if ( fabs(otherEle_etaDiff) < params.otherElesEtaSlice) vetoedByVetoEles=true;  // jurassic strip cut
							if ( otherEle_dR < params.otherElesIntRadius) vetoedByVetoEles=true; // jurassic exclusion cone cut
						}
						if(vetoedByVetoEles)
							break;
					}
					if(vetoedByVetoEles)
						continue;

					// RecHit severity level / 'reco' flag checks
					int severityFlag = ecalBarHits_ == 0 ? -1 : sevLevelAlgo->severityLevel(((EcalRecHit*)(&*j))->detid(), *ecalBarHits_);
					std::vector<int>::const_iterator sit = std::find(params.recHitSeveritiesExcl_.begin(), 
											 params.recHitSeveritiesExcl_.end(), severityFlag);
					if( sit != params.recHitSeveritiesExcl_.end() )
						continue;

					std::vector<int>::const_iterator vit = std::find( params.recHitFlagsExcl_.begin(), params.recHitFlagsExcl_.end(),  ((EcalRecHit*)(&*j))->recoFlag() );
					if ( vit != params.recHitFlagsExcl_.end() )
						continue; // the recHit has to be excluded from the iso sum


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
	//math::XYZPoint theCaloPosition = sc->position();
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

		double dr2 = deltaEta*deltaEta + deltaPhi*deltaPhi;

		bool includeInIsolSum=false;
		if( dr2 < extRadius_*extRadius_ && dr2 >= hcal_intRadius_*hcal_intRadius_ )
			includeInIsolSum = true;

		// *** Check if caloTower is within inner veto area of any of the other eles
		//   If extraPhantomVetoEle_ true, then instead just check if caloTower is w/i inner veto area of this phantom ele
		for( reco::GsfElectronCollection::const_iterator otherEleIt = vetoEles.begin(); otherEleIt!=vetoEles.end(); otherEleIt++){
			const reco::SuperCluster* otherEleSC = otherEleIt->get<reco::SuperClusterRef>().get();
			//math::XYZPoint otherEleCaloPosition = otherEleSC->position();
			double otherEle_eta = extraPhantomVetoEle_ ? (candEta+phantomVetoEleDEta_) : otherEleSC->eta();
			double otherEle_phi = extraPhantomVetoEle_ ? (candPhi+phantomVetoEleDPhi_) : otherEleSC->phi();

			double otherEle_dEta = towerEta - otherEle_eta;
			double otherEle_dPhi = reco::deltaPhi(towerPhi, otherEle_phi);
			double otherEle_dR   = sqrt(otherEle_dEta*otherEle_dEta + otherEle_dPhi*otherEle_dPhi);

			if(otherEle_dR < hcal_intRadius_)
				includeInIsolSum=false;

			if(!includeInIsolSum)
				break;
		}

		if(includeInIsolSum)
			ptSum += this_pt;

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
