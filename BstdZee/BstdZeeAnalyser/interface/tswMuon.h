#ifndef tswMuon_h
#define tswMuon_h

#include "TSWilliams/BstdZeeNTupler/interface/tswMuStruct.h"
#include "TSWilliams/BstdZeeNTupler/interface/tswUsefulFunctions.h"

namespace tsw{
	class Muon{
		private:
			MuStruct muStr_;
			TLorentzVector stdP4_;  ///< 4-momentum returned by reco::Muon::p4
			TLorentzVector corrP4_; ///< 4-momentum of reco muon, corrected to pT of TuneP (i.e. 'best') track if exists
		public:
			//CTOR and DTOR ...
			Muon(): muStr_()
				{}
			Muon(tsw::MuStruct struct_theMuon) :
				muStr_(struct_theMuon),
				stdP4_( ConvertToTLorentzVector( &(muStr_.p4) ) ),
				corrP4_(bestTrk_exists() ? bestTrk_pT()/stdP4_.Pt()*stdP4_ : stdP4_)
			{}
			~Muon(){}

			// Methods for access to the cut variables ...
			const TLorentzVector& p4() const { return corrP4_; }
			double p()   const {   return p4().P(); }
			double pT()  const {  return p4().Pt(); }
			double eta() const { return p4().Eta(); }
			double phi() const { return p4().Phi(); }
			const TLorentzVector& stdRecoP4() const { return stdP4_; }
			int charge() const { return muStr_.charge; }
			bool isGlobalMuon()      const { return muStr_.isGlobalMuon; }
			bool isTrackerMuon()     const { return muStr_.isTrackerMuon; }
			bool isStandAloneMuon()  const { return muStr_.isStandAloneMuon; }
			int numMatchedMuonStns() const { return muStr_.numMatchedMuonStns; }
			double isolR03_sumPt()   const { return muStr_.isolR03_sumPt; }

			bool glob_exists() const { return muStr_.globTrk_exists; }
			double glob_pT()   const { return muStr_.globTrk_pT; }
	   	double glob_eta()  const { return muStr_.globTrk_eta; }
	   	double glob_phi()  const { return muStr_.globTrk_phi; }
	   	int glob_charge()  const { return muStr_.globTrk_charge; }
	   	int glob_numValidMuonHits() const { return muStr_.globTrk_numberOfValidMuonHits; }
	   	double glob_normalisedChi2() const { return muStr_.globTrk_normalisedChi2; }

	   	bool inner_exists() const { return muStr_.inTrk_exists; }
	   	double inner_pT()   const { return muStr_.inTrk_pT; }
	   	double inner_eta()  const { return muStr_.inTrk_eta; }
	   	double inner_phi()  const { return muStr_.inTrk_phi; }
	   	int inner_charge()  const { return muStr_.inTrk_charge;}
	   	int inner_numValidPixHits()  const { return muStr_.inTrk_numValidPixHits; }
	   	int inner_numValidTrkrHits() const { return muStr_.inTrk_numValidTrkrHits; }
	   	double inner_dxyVsOrigin()   const { return muStr_.inTrk_dxyVsOrigin; }

	   	int trk_trkrLayerWHits() const { return muStr_.trk_trkrLayersWHits; }

	   	bool outer_exists()  const { return muStr_.outTrk_exists; }
	   	double outer_pT()    const { return muStr_.outTrk_pT; }
	   	double outer_eta()   const { return muStr_.outTrk_eta; }
	   	double outer_phi()   const { return muStr_.outTrk_phi; }
	   	int outer_charge()   const { return muStr_.outTrk_charge; }

	   	bool   bestTrk_exists()    const { return muStr_.bestTrk_exists; }
	   	double bestTrk_pT()        const { return muStr_.bestTrk_pT; }
	   	double bestTrk_ptError()   const { return muStr_.bestTrk_ptError; }
	   	double bestTrk_dxy_bspot() const { return muStr_.bestTrk_dxy_bspot; }
			double bestTrk_dxy_vtx()   const { return muStr_.bestTrk_dxy_vtx; }
			double bestTrk_dz_vtx()    const { return muStr_.bestTrk_dz_vtx; }

	   	// Methods for which region the muon is in - ALL regions are MUTUALLY exclusive ...
	   	bool isInBarrel() const { return (fabs(eta())<=0.9); }
	   	bool isInOverlap() const { return (fabs(eta())>0.9 && fabs(eta())<=1.2); }
	   	bool isInEndcap()  const { return (fabs(eta())>1.2 && fabs(eta())<=2.4); }
	   	bool isWithinCMS() const { return (fabs(eta())<=2.4); }
	   	bool isWithinAcc() const { return isWithinCMS(); }

			// Methods for applying cuts to the muon ...
			int highPtIdCutCode(float minPt, float maxEta) const;
			int highPtIdCut(float minPt, float maxEta) const { return (highPtIdCutCode(minPt, maxEta)==0); }

			//  STATIC CONSTS - CUT CODES  //

			static const int mCutCode_fiducial    = 0x0001;
			static const int mCutCode_isGlobal    = 0x0002;
			static const int mCutCode_isPFMuon    = 0x0004;
			static const int mCutCode_normChi2    = 0x0008;
			static const int mCutCode_nrValidHits = 0x0010;
			static const int mCutCode_nrStations  = 0x0020;
			static const int mCutCode_bestTrkDxy  = 0x0040;
			static const int mCutCode_bestTrkDz   = 0x0080;
			static const int mCutCode_nrPixHits   = 0x0100;
			static const int mCutCode_nrTkLayers  = 0x0200;
			static const int mCutCode_dPtOverPt   = 0x0400;
	};
}//end namespace tsw


int tsw::Muon::highPtIdCutCode(float minPt, float maxEta) const
{
	int cutCode = 0;

	if( fabs(eta())>maxEta || pT()<minPt  )
		cutCode |= mCutCode_fiducial;

	// is global muon
	if( !isGlobalMuon() )
		cutCode |= mCutCode_isGlobal;

	// is PF muon -- NOT IN HIGH PT
	// norm'd chi^2 -- NOT IN HIGH PT

	// number muon chamber hits in global track-fit
	if( !glob_exists() || glob_numValidMuonHits()<=0 )
		cutCode |= mCutCode_nrValidHits;

	// Num muon segements in muon stations
	if( numMatchedMuonStns()<=1 )
		cutCode |= mCutCode_nrStations;

	// Best track impact parameters
	if( !bestTrk_exists() || fabs(bestTrk_dxy_vtx())>=0.2 )
		cutCode |= mCutCode_bestTrkDxy;
	if( !bestTrk_exists() || fabs(bestTrk_dz_vtx())>=0.5 )
		cutCode |= mCutCode_bestTrkDz;

	// Number pixel hits
	if( !inner_exists() || inner_numValidPixHits()<=0 )
		cutCode |= mCutCode_nrPixHits;

	// Number tracker layers with hits
	if( trk_trkrLayerWHits()<=5 )
		cutCode |= mCutCode_nrTkLayers;

	// dPt/Pt
	if( !bestTrk_exists() || (bestTrk_ptError()/bestTrk_pT())>=0.3 )
		cutCode |= mCutCode_dPtOverPt;

	return cutCode;
}

#endif
