
#include "NTupler/BstdZeeNTupler/interface/tswEvent.h"
#include <iostream>

namespace tsw{

	Event::Event() :
		runNum_(0), lumiSec_(0), evtNum_(0),
		normMuons_charge_(), normMuons_isGlobalMuon_(), normMuons_isTrackerMuon_(), normMuons_isStandAloneMuon_(), normMuons_numMatchedMuonStns_(), normMuons_isolR03_sumPt_(),
		normMuons_globTrk_exists_(), normMuons_globTrk_pT_(), normMuons_globTrk_eta_(), normMuons_globTrk_phi_(), normMuons_globTrk_charge_(), normMuons_globTrk_numberOfValidMuonHits_(), normMuons_globTrk_normalisedChi2_(),
		normMuons_inTrk_exists_(), normMuons_inTrk_pT_(), normMuons_inTrk_eta_(), normMuons_inTrk_phi_(), normMuons_inTrk_charge_(), normMuons_inTrk_numValidPixHits_(), normMuons_inTrk_numValidTrkrHits_(), normMuons_inTrk_dxyVsOrigin_(),
		normMuons_outTrk_exists_(), normMuons_outTrk_pT_(), normMuons_outTrk_eta_(), normMuons_outTrk_phi_(), normMuons_outTrk_charge_()
	{

	}
	Event::~Event() {}

	void Event::SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber){
		runNum_  = runNumber;
		lumiSec_ = lumiSection;
		evtNum_  = eventNumber;
	}
	void Event::PrintBasicEventInformation(){
		std::cout << "Run " << runNum_ << ", LumiSec " << lumiSec_ << ", event " << evtNum_ << std::endl;
	}

	void Event::SetEMuTriggerInfo(const std::string eMuPathName, const bool eMuPathDecision)
	{
		trg_emuPath_name_ = eMuPathName;
		trg_emuPath_decision_ = eMuPathDecision;
	}

	void Event::AddNormMuon(tsw::MuStruct* theMuon){
		normMuons_p4_                .push_back(theMuon->p4);
		normMuons_charge_            .push_back(theMuon->charge);
		normMuons_isGlobalMuon_      .push_back(theMuon->isGlobalMuon);
		normMuons_isTrackerMuon_     .push_back(theMuon->isTrackerMuon);
		normMuons_isStandAloneMuon_  .push_back(theMuon->isStandAloneMuon);
		normMuons_numMatchedMuonStns_.push_back(theMuon->numMatchedMuonStns);
		normMuons_isolR03_sumPt_     .push_back(theMuon->isolR03_sumPt);

		normMuons_globTrk_exists_.push_back(theMuon->globTrk_exists);
		normMuons_globTrk_pT_    .push_back(theMuon->globTrk_pT);
		normMuons_globTrk_eta_   .push_back(theMuon->globTrk_eta);
		normMuons_globTrk_phi_   .push_back(theMuon->globTrk_phi);
		normMuons_globTrk_charge_.push_back(theMuon->globTrk_charge);
		normMuons_globTrk_numberOfValidMuonHits_.push_back(theMuon->globTrk_numberOfValidMuonHits);
		normMuons_globTrk_normalisedChi2_.push_back(theMuon->globTrk_normalisedChi2);

		normMuons_inTrk_exists_         .push_back(theMuon->inTrk_exists);
		normMuons_inTrk_pT_              .push_back(theMuon->inTrk_pT);
		normMuons_inTrk_eta_             .push_back(theMuon->inTrk_eta);
		normMuons_inTrk_phi_             .push_back(theMuon->inTrk_phi);
		normMuons_inTrk_charge_          .push_back(theMuon->inTrk_charge);
		normMuons_inTrk_numValidPixHits_ .push_back(theMuon->inTrk_numValidPixHits);
		normMuons_inTrk_numValidTrkrHits_.push_back(theMuon->inTrk_numValidTrkrHits);
		normMuons_inTrk_dxyVsOrigin_     .push_back(theMuon->inTrk_dxyVsOrigin);

		normMuons_outTrk_exists_.push_back(theMuon->outTrk_exists);
		normMuons_outTrk_pT_    .push_back(theMuon->outTrk_pT);
		normMuons_outTrk_eta_   .push_back(theMuon->outTrk_eta);
		normMuons_outTrk_phi_   .push_back(theMuon->outTrk_phi);
		normMuons_outTrk_charge_.push_back(theMuon->outTrk_charge);
	}
}


#if !defined(__CINT__)
  ClassImp(tsw::Event)
#endif
