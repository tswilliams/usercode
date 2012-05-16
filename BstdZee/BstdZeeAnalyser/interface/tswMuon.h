#ifndef tswMuon_h
#define tswMuon_h

#include "BstdZeeFirst/Analyser/interface/tswMuStruct.h"
#include "BstdZeeFirst/Analyser/interface/tswUsefulFunctions.h"

namespace tsw{
	class Muon{
		private:
			MuStruct muStr_;
		public:
			//CTOR and DTOR ...
			Muon(): muStr_()
				{}
			Muon(tsw::MuStruct struct_theMuon): muStr_(struct_theMuon)
				{}
			~Muon(){}

			// Methods for access to the cut variables ...
			TLorentzVector p4(){
				return ConvertToTLorentzVector( &(muStr_.p4) );}
			double p(){   return p4().P(); }
			double pT(){  return p4().Pt(); }
			double eta(){ return p4().Eta(); }
			double phi(){ return p4().Phi(); }
			int charge(){
				return muStr_.charge;}
			bool isGlobalMuon(){
				return muStr_.isGlobalMuon;}
			bool isTrackerMuon(){
				return muStr_.isTrackerMuon;}
			bool isStandAloneMuon(){
				return muStr_.isStandAloneMuon;}
			int numMatchedMuonStns(){
				return muStr_.numMatchedMuonStns;}
			double isolR03_sumPt(){
				return muStr_.isolR03_sumPt;}

			bool glob_exists(){
				return muStr_.globTrk_exists;}
			double glob_pT(){
				return muStr_.globTrk_pT;}
	   	double glob_eta(){
	   		return muStr_.globTrk_eta;}
	   	double glob_phi(){
	   		return muStr_.globTrk_phi;}
	   	int glob_charge(){
	   		return muStr_.globTrk_charge;}
	   	int glob_numValidMuonHits(){
	   		return muStr_.globTrk_numberOfValidMuonHits;}
	   	double glob_normalisedChi2(){
	   		return muStr_.globTrk_normalisedChi2;}

	   	bool inner_exists(){
	   		return muStr_.inTrk_exists;}
	   	double inner_pT(){
	   		return muStr_.inTrk_pT;}
	   	double inner_eta(){
	   		return muStr_.inTrk_eta;}
	   	double inner_phi(){
	   		return muStr_.inTrk_phi;}
	   	int inner_charge(){
	   		return muStr_.inTrk_charge;}
	   	int inner_numValidPixHits(){
	   		return muStr_.inTrk_numValidPixHits;}
	   	int inner_numValidTrkrHits(){
	   		return muStr_.inTrk_numValidTrkrHits;}
	   	double inner_dxyVsOrigin(){
	   		return muStr_.inTrk_dxyVsOrigin;}

	   	bool outer_exists(){
	   		return muStr_.outTrk_exists;}
	   	double outer_pT(){
	   		return muStr_.outTrk_pT;}
	   	double outer_eta(){
	   		return muStr_.outTrk_eta;}
	   	double outer_phi(){
	   		return muStr_.outTrk_phi;}
	   	int outer_charge(){
	   		return muStr_.outTrk_charge;}

	   	// Methods for which region the muon is in - ALL regions are MUTUALLY exclusive ...
	   	bool isInBarrel(){
	   		if( fabs(eta())<=0.9)
	   			return true;
	   		else
	   			return false;
	   	}
	   	bool isInOverlap(){
	   		if( fabs(eta())>0.9 && fabs(eta())<=1.2 )
	   			return true;
	   		else
	   			return false;
	   	}
	   	bool isInEndcap(){
	   		if( fabs(eta())>1.2 && fabs(eta())<=2.4 )
	   			return true;
	   		else
	   			return false;
	   	}
	   	bool isWithinCMS(){
	   		if( fabs(eta())<=2.4 )
	   			return true;
	   		else
	   			return false;
	   	}

			// Methods for applying cuts to the muon ...
			bool isWithinAcc();
			bool isTightMuon();
	};

	bool Muon::isWithinAcc(){
		bool passFailFlag = true;
		passFailFlag = passFailFlag && isWithinCMS();
		return passFailFlag;
	}

	bool Muon::isTightMuon(){
		bool passFailFlag = true;

		passFailFlag = passFailFlag && ( fabs(eta())<1.442 );
		passFailFlag = passFailFlag && ( pT()>35.0 );
		passFailFlag = passFailFlag && isGlobalMuon();
		passFailFlag = passFailFlag && ( glob_normalisedChi2()<10 );
		passFailFlag = passFailFlag && ( glob_numValidMuonHits()>0 );
		passFailFlag = passFailFlag && ( numMatchedMuonStns()>1 );
		// Isolation ...
		passFailFlag = passFailFlag && ( isolR03_sumPt()<10.0 );
		// Tracker stuff ...
		passFailFlag = passFailFlag && ( inner_exists() );
		passFailFlag = passFailFlag && ( inner_numValidPixHits()>0 );
		passFailFlag = passFailFlag && ( inner_numValidTrkrHits()>10 );
		// --- N.B: dxy cut has temporarily been removed as it is highly inefficient under the current way that I'm calculating the value of dxy
		// -------- [i.e. calculating dxy rel. to (0,0) RATHER THAN prim vertex/beamspot ]
		// -------- This cut also currently results in weird phi distributions (i.e. tight muons only having phi from 0.4->1.4 and -1.7 to -2.7)
		//passFailFlag = passFailFlag && ( fabs(inner_dxyVsOrigin())<0.2 ); // N.B: It may be this cut that's responsible for weird behaviour

		return passFailFlag;
	}
}

#endif
