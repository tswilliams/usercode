#ifndef tswEventHelper_h
#define tswEventHelper_h

#include "tswMuonCollection.h"
#include "NTupler/BstdZeeNTupler/interface/tswEvent.h"

namespace tsw{
	class EventHelper{
		private:
			tsw::Event* theEvent_;
		public:
			//CTORs and DTOR ...
			EventHelper() : theEvent_(0) {}
			EventHelper(tsw::Event* ptr_event) : theEvent_(ptr_event){}
			// Methods for getting information out of tsw::Event ...
			tsw::MuonCollection GetNormMuons();
	};

	tsw::MuonCollection EventHelper::GetNormMuons(){
		std::vector<tsw::Muon> vecOfMuons; vecOfMuons.clear();
		unsigned int numOfMuons = theEvent_->normMuons_charge_.size();

		//Loop over the normMuon vectors in the event, storing the information for each muon in a tsw::MuStruct ...
		//... and using the muStruct to initialise instances of tsw::Muon ...
		for(unsigned int muIdx=0; muIdx<numOfMuons; muIdx++){
			MuStruct ithMuStruct;
			// General variables ...
			ithMuStruct.p4                 = theEvent_->normMuons_p4_.at(muIdx);
			ithMuStruct.charge             = theEvent_->normMuons_charge_.at(muIdx);
			ithMuStruct.isGlobalMuon       = theEvent_->normMuons_isGlobalMuon_.at(muIdx);
			ithMuStruct.isTrackerMuon      = theEvent_->normMuons_isTrackerMuon_.at(muIdx);
			ithMuStruct.isStandAloneMuon   = theEvent_->normMuons_isStandAloneMuon_.at(muIdx);
			ithMuStruct.numMatchedMuonStns = theEvent_->normMuons_numMatchedMuonStns_.at(muIdx);
			ithMuStruct.isolR03_sumPt      = theEvent_->normMuons_isolR03_sumPt_.at(muIdx);
			// Kinematic variables from the global track ...
			ithMuStruct.globTrk_exists = theEvent_->normMuons_globTrk_exists_.at(muIdx);
			ithMuStruct.globTrk_pT     = theEvent_->normMuons_globTrk_pT_.at(muIdx);
			ithMuStruct.globTrk_eta    = theEvent_->normMuons_globTrk_eta_.at(muIdx);
			ithMuStruct.globTrk_phi    = theEvent_->normMuons_globTrk_phi_.at(muIdx);
			ithMuStruct.globTrk_charge = theEvent_->normMuons_globTrk_charge_.at(muIdx);
			ithMuStruct.globTrk_numberOfValidMuonHits = theEvent_->normMuons_globTrk_numberOfValidMuonHits_.at(muIdx);
			ithMuStruct.globTrk_normalisedChi2 = theEvent_->normMuons_globTrk_normalisedChi2_.at(muIdx);
			// ... and from the inner track ...
			ithMuStruct.inTrk_exists           = theEvent_->normMuons_inTrk_exists_.at(muIdx);
			ithMuStruct.inTrk_pT               = theEvent_->normMuons_inTrk_pT_.at(muIdx);
			ithMuStruct.inTrk_eta              = theEvent_->normMuons_inTrk_eta_.at(muIdx);
			ithMuStruct.inTrk_phi              = theEvent_->normMuons_inTrk_phi_.at(muIdx);
			ithMuStruct.inTrk_charge           = theEvent_->normMuons_inTrk_charge_.at(muIdx);
			ithMuStruct.inTrk_numValidPixHits  = theEvent_->normMuons_inTrk_numValidPixHits_.at(muIdx);
			ithMuStruct.inTrk_numValidTrkrHits = theEvent_->normMuons_inTrk_numValidTrkrHits_.at(muIdx);
			ithMuStruct.inTrk_dxyVsOrigin      = theEvent_->normMuons_inTrk_dxyVsOrigin_.at(muIdx);
	   	// ... and from the outer track ...
			ithMuStruct.outTrk_exists = theEvent_->normMuons_outTrk_exists_.at(muIdx);
			ithMuStruct.outTrk_pT     = theEvent_->normMuons_outTrk_pT_.at(muIdx);
			ithMuStruct.outTrk_eta    = theEvent_->normMuons_outTrk_eta_.at(muIdx);
			ithMuStruct.outTrk_phi    = theEvent_->normMuons_outTrk_phi_.at(muIdx);
			ithMuStruct.outTrk_charge = theEvent_->normMuons_outTrk_charge_.at(muIdx);

			Muon ithMuon(ithMuStruct);
			vecOfMuons.push_back(ithMuon);
		}
		MuonCollection theMuColln(&vecOfMuons, "normal muons");
		return theMuColln;
	}
}

#endif
