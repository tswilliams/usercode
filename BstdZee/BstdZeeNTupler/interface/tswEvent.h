
#ifndef tswEvent_h
#define tswEvent_h

#include "TObject.h"

#include <vector>

#include "NTupler/BstdZeeNTupler/interface/tswMuStruct.h"
//#include "BstdZeeFirst/Analyser/interface/tswEventHelper.h"

namespace tsw{
	class Event : public TObject {
		friend class EventHelper;
		public:
			Event();
			~Event();

			//Methods to set information ...
			void SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber);
			void PrintBasicEventInformation();
			void AddNormMuon(tsw::MuStruct* theMuon);

		private:
			//General event information ...
			unsigned int runNum_;
			unsigned int lumiSec_;
			unsigned int evtNum_;

			//Information about the standard muons ...
			std::vector<ROOT::Math::XYZTVector> normMuons_p4_;
			std::vector<Int_t>    normMuons_charge_;
			std::vector<Bool_t>   normMuons_isGlobalMuon_;
			std::vector<Bool_t>   normMuons_isTrackerMuon_;
			std::vector<Bool_t>   normMuons_isStandAloneMuon_;
			std::vector<Int_t>    normMuons_numMatchedMuonStns_;
			std::vector<Double_t> normMuons_isolR03_sumPt_;

			std::vector<Bool_t>   normMuons_globTrk_exists_;
			std::vector<Double_t> normMuons_globTrk_pT_;
			std::vector<Double_t> normMuons_globTrk_eta_;
			std::vector<Double_t> normMuons_globTrk_phi_;
			std::vector<Int_t>    normMuons_globTrk_charge_;
			std::vector<Int_t>    normMuons_globTrk_numberOfValidMuonHits_;
			std::vector<Double_t> normMuons_globTrk_normalisedChi2_;

			std::vector<Bool_t>   normMuons_inTrk_exists_;
			std::vector<Double_t> normMuons_inTrk_pT_;
			std::vector<Double_t> normMuons_inTrk_eta_;
			std::vector<Double_t> normMuons_inTrk_phi_;
			std::vector<Int_t>    normMuons_inTrk_charge_;
			std::vector<Int_t>    normMuons_inTrk_numValidPixHits_;
			std::vector<Int_t>    normMuons_inTrk_numValidTrkrHits_;
			std::vector<Double_t> normMuons_inTrk_dxyVsOrigin_;

			std::vector<Bool_t>   normMuons_outTrk_exists_;
			std::vector<Double_t> normMuons_outTrk_pT_;
			std::vector<Double_t> normMuons_outTrk_eta_;
			std::vector<Double_t> normMuons_outTrk_phi_;
			std::vector<Int_t>    normMuons_outTrk_charge_;

			ClassDef(tsw::Event,5);
	};
}

#endif
