
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

			// Enum for size of inner vetos
			enum InnerVetoSize{
				xSmallVeto = 0,
				smallVeto = 1,
				mediumVeto = 2,
				largeVeto = 3,
				xLargeVeto = 4
			};

			//Methods to set information ...
			void SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber);
			void SetMCGenWeight(float genWeight);
			void SetLHEZbosonInfo(ROOT::Math::XYZTVector p4_Zboson);
			void SetMCPUInfo(float nVtxPoissonMean, unsigned int nVertices, std::vector<float> vtxZposns, double puWeight_1D);
			void SetRecoVtxInfo(unsigned int nRecoVtxs, unsigned int nGoodRecoVtxs);
			void SetPURho(const double rhoValue);
			void SetEMuTriggerInfo(const std::string eMuPathName, const bool eMuPathDecision);
			void AddStdEleInfo_isoDep_std(double ,double ,double );
			void AddStdEleInfo_isoDep_inrVeto(double, double, double );
			void AddStdEleInfo_inrVetoModIso(tsw::Event::InnerVetoSize, double, double, double );
			void AddStdEleInfo_genHadronsDr04(unsigned int nGenHadrons_dR04, double ptSumGenHadrons_dR04);

			void AddNormMuon(tsw::MuStruct* theMuon);

			/// Prints out basic event information to screen
			void PrintBasicEventInformation();
			void PrintMCGenWeight();
			void PrintLHEZbosonInfo();
			void PrintPUVtxInfo();
			void PrintPURho();
			void PrintStdEleInfo_isoValues(unsigned int iEle);

		private:
			//General event information ...
			unsigned int runNum_;
			unsigned int lumiSec_;
			unsigned int evtNum_;

			// Generator information from MC ...
			float mc_genWeight_;
			ROOT::Math::XYZTVector mcLHE_ZbosonP4_;

			// Pile up information from MC ...
			float mc_bx0_nPUVtxPoissonMean_;
			unsigned int mc_bx0_nPUVtx_;
			std::vector<float> mc_bx0_PUVtxZPosns_;
			double mc_puWeight_1D_;

			// Reconstructed vertex info ...
			unsigned int recoVtx_totalNum_;
			unsigned int recoVtx_numGoodVtxs_;

			// Isol value PU correction info ...
			double pu_rho_;

			// Information about e-mu trigger ...
			std::string trg_emuPath_name_;
			bool trg_emuPath_decision_;

			// Information about the standard electrons ...
			std::vector<Double_t> stdEles_isoDeps_stdTrkIso_;
			std::vector<Double_t> stdEles_isoDeps_stdEcalIso_;
			std::vector<Double_t> stdEles_isoDeps_stdHcalD1Iso_;
			std::vector<Double_t> stdEles_isoDeps_inrVetoModTrkIso_;
			std::vector<Double_t> stdEles_isoDeps_inrVetoModEcalIso_;
			std::vector<Double_t> stdEles_isoDeps_inrVetoModHcalD1Iso_;

			std::vector<Double_t> stdEles_inrVetoModIso_Trk_;
			std::vector<Double_t> stdEles_inrVetoModIso_Ecal_;
			std::vector<Double_t> stdEles_inrVetoModIso_HcalD1_;
			std::vector<Double_t> stdEles_inrVetoXSModIso_Trk_;
			std::vector<Double_t> stdEles_inrVetoXSModIso_Ecal_;
			std::vector<Double_t> stdEles_inrVetoXSModIso_HcalD1_;
			std::vector<Double_t> stdEles_inrVetoSModIso_Trk_;
			std::vector<Double_t> stdEles_inrVetoSModIso_Ecal_;
			std::vector<Double_t> stdEles_inrVetoSModIso_HcalD1_;
			std::vector<Double_t> stdEles_inrVetoLModIso_Trk_;
			std::vector<Double_t> stdEles_inrVetoLModIso_Ecal_;
			std::vector<Double_t> stdEles_inrVetoLModIso_HcalD1_;
			std::vector<Double_t> stdEles_inrVetoXLModIso_Trk_;
			std::vector<Double_t> stdEles_inrVetoXLModIso_Ecal_;
			std::vector<Double_t> stdEles_inrVetoXLModIso_HcalD1_;

			std::vector<UInt_t>   stdEles_nGenHadronsDr04_;
			std::vector<Double_t> stdEles_ptSumGenHadronsDr04_;

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

			ClassDef(tsw::Event,12);
	};
}

#endif
