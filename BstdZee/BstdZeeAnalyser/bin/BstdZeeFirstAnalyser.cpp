//Include files...
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TMath.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h" // Included for the typedef of ROOT::Math::LorentzVector<PxPyPzE4D<double> > to ROOT::Math::XYZTVector
//#include "Math/PxPyPzE4D.h"
//#include "Math/GenVector/GenVector_exception.h"
#include "TLorentzVector.h"
#include <cmath>
//Must include TROOT.h in order for the gROOT->ProcessLine("#inlcude <vector>") line to work.
#include "TROOT.h" 
//    gROOT->ProcessLine("#include <vector>");

//#include "NTupler/BstdZeeNTupler/interface/tswEvent.h"
#include "BstdZeeFirst/Analyser/interface/tswEventHelper.h"

#include "NTupler/BstdZeeNTupler/interface/tswHEEPEle.h"
#include "BstdZeeFirst/Analyser/interface/tswReconValidationHistos.h"
#include "BstdZeeFirst/Analyser/interface/tswUsefulFunctions.h"
#include "BstdZeeFirst/Analyser/interface/tswHEEPDiEle.h"
#include "BstdZeeFirst/Analyser/interface/tswDiEleDistns.h"
#include "BstdZeeFirst/Analyser/interface/tswMCParticle.h"
#include "BstdZeeFirst/Analyser/interface/tswMCParticleDistns.h"
#include "BstdZeeFirst/Analyser/interface/tswMuonDistns.h"

#include "BstdZeeFirst/Analyser/interface/tswEleMuObject.h"

//#include "linkdef.h"
//#pragma link C++ class vector<float> + ;
//#pragma link C++ class vector<double> + ;
//#pragma link C++ class vector<int> + ;
//#pragma link C++ class vector<bool> + ;

#ifdef __MAKECINT__
#pragma link C++ class vector<Int_t>+;
#endif

//Functions...


//=====================================================================================================
//-----------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------------------
//=========================================== Analyser class ==========================================

class BstdZeeFirstAnalyser{
	public:
		BstdZeeFirstAnalyser(int runMode, unsigned int numEvts, bool isMC, const TString& inFileName, const TString& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax);
		~BstdZeeFirstAnalyser();
		void DoAnalysis(const Double_t evtWeight);
		
	private:
		void AnalyseEvent(const Double_t );
		void SetupEleClassVectors();
		void SetupMuonCollection();
		void PrintOutBranchVariables();
		void FinishOffAnalysis(); //Method to be called after all events analysed - i.e. for normalisation of histograms, calculating errors on each bins, etc

		//Methods for ordering reconstructed objects in terms of kinematic variables...
		std::vector <Int_t> OrderByEt(const std::vector<Double_t>& etValues);
		void TwoHighestEtObjects(const std::vector<Double_t>& etValues, int& idx_HighestEtObject, int& idx_2ndHighestEtObject);

		//Methods for application of HEEP cuts ...
		std::vector<bool> HEEPCuts(int eleType); //Apply to normGSFEles if eleType=0, bstdGSFEles if eleType=1, and return empty vector (and error message) otherwise
		std::vector<bool> HEEPCutsWithoutIso(int eleType);

		//Methods for filling sets of histograms...
		void FillReconValidationHistos(const Double_t );
		void FillDiEleHistos(const Double_t );
		void FillHistograms(const Double_t );
		void DoEMuAnalysis(const Double_t );

		//Methods for storing histograms in output file...
		void SaveReconValidationHistos(TFile* histosFile);

		//Methods called by the constructor (e.g. Initialise methods)...
		void SetMemberVariablesToDefaultValues();
		void InitialiseReconValidationHistos();
		void SetupBranchLinks(const TFile*);
		void SetupEMuMethodTrees();

		//Methods called by the destructor (i.e. deletion of heap-based objects)...
		void DeleteReconValidationHistos();
	
	//Member variables...
	private:
		int runMode_;
		int vFlg_;
		unsigned int numEvts_;
		bool isMC_;

		//Input and output files...
		TString inputFile_name_;
		TTree* inputFile_tree_;
		TString outputFile_name_;

		//------------//
		//Variables for storing contents of event information branches...
		tsw::Event* event_;
		tsw::EventHelper eventHelper_;
		unsigned int evt_runNum_;
		unsigned int evt_lumiSec_;
		unsigned int evt_evtNum_;

		//------------//
		//Variables for storing contents of trigger information branches...
		Bool_t trg_PathA_decision_;
		std::string* trg_PathA_name_ptr_;
		std::string* trg_PathA_nameOfLastFilter_ptr_;
		float trg_PathA_highestTrigObjEt_;

		//------------//
		//Variables for storing contents of MC final state particle branches...
		unsigned int mc_numFinalStateEles_;

		Int_t mcEles_HighestEt_charge_;
		Int_t mcEles_HighestEt_PDGid_;
		Int_t mcEles_HighestEt_status_;
		Double_t mcEles_HighestEt_pt_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcEles_HighestEt_OLDp4_ptr_; //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcEles_HighestEt_p4_;

		Int_t mcEles_2ndHighestEt_charge_;
		Int_t mcEles_2ndHighestEt_PDGid_;
		Int_t mcEles_2ndHighestEt_status_;
		Double_t mcEles_2ndHighestEt_pt_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcEles_2ndHighestEt_OLDp4_ptr_; //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcEles_2ndHighestEt_p4_;

		Double_t mcZcandidate_pt_;
		Double_t mcZcandidate_eta_;
		Double_t mcZcandidate_phi_;
		Double_t mcZcandidate_mass_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcZcandidate_OLDp4_ptr_; //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcZcandidate_p4_;
		Double_t mcZcandidate_dEtaEles_;
		Double_t mcZcandidate_dPhiEles_;
		Double_t mcZcandidate_dREles_;
		Double_t mcZcandidate_openingAngle_;

		Int_t mcZ_pdgId_;
		Int_t mcZ_status_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcZ_p4_ptr_; //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		Int_t mcZ_numDaughters_;
		Double_t mcZ_daughters_dR_;
		Double_t mcZ_daughters_dEta_;
		Double_t mcZ_daughters_dPhi_;
		Double_t mcZ_daughters_openingAngle_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcZ_daughterA_OLDp4_ptr_; //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcZ_daughterA_p4_;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcZ_daughterB_OLDp4_ptr_; //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcZ_daughterB_p4_;
		tsw::MCParticle mcZboson_;

		//------------//
		//Variables for storing contents of normal GSF electron info branches...
		unsigned int normGsfEles_number_;
		std::vector<ROOT::Math::XYZTVector>* normGsfEles_OLDp4_;
		std::vector<TLorentzVector> normGsfEles_p4_;
		std::vector<Int_t>* normGsfEles_charge_;

		std::vector<Double_t>* normGsfEles_Et_;
		std::vector<Double_t>* normGsfEles_HEEP_Et_;
		std::vector<Double_t>* normGsfEles_Eta_;
		std::vector<Double_t>* normGsfEles_scEta_;
		std::vector<bool>*     normGsfEles_ecalDriven_;
		std::vector<bool>*     normGsfEles_ecalDrivenSeed_;

		std::vector<Double_t>* normGsfEles_HEEP_dEtaIn_;
		std::vector<Double_t>* normGsfEles_HEEP_dPhiIn_;
		std::vector<Double_t>* normGsfEles_HoverE_;
		std::vector<Double_t>* normGsfEles_sigmaIetaIeta_;
		std::vector<Double_t>* normGsfEles_scSigmaIetaIeta_;

		std::vector<Double_t>* normGsfEles_dr03EmIsoEt_;
		std::vector<Double_t>* normGsfEles_dr03HadDepth1IsoEt_;
		std::vector<Double_t>* normGsfEles_dr03HadDepth2IsoEt_;
		std::vector<Double_t>* normGsfEles_dr03TkIsoPt_;

		std::vector<Double_t>* normGsfEles_e2x5Max_;
		std::vector<Double_t>* normGsfEles_e5x5_;

		//kinematic and geometric methods
		std::vector<float>* normHEEPEles_et_;
		std::vector<float>* normHEEPEles_gsfEt_;
		std::vector<float>* normHEEPEles_scEt_;
		std::vector<float>* normHEEPEles_energy_;
		std::vector<float>* normHEEPEles_gsfEnergy_;
		std::vector<float>* normHEEPEles_caloEnergy_;
		std::vector<float>* normHEEPEles_eta_;
		std::vector<float>* normHEEPEles_scEta_;
		std::vector<float>* normHEEPEles_detEta_;
		std::vector<float>* normHEEPEles_detEtaAbs_;
		std::vector<float>* normHEEPEles_phi_;
		std::vector<float>* normHEEPEles_scPhi_;
		std::vector<float>* normHEEPEles_detPhi_;
		std::vector<float>* normHEEPEles_zVtx_;
		std::vector<ROOT::Math::XYZTVector>* normHEEPEles_p4ptr_;
		std::vector<ROOT::Math::XYZTVector>* normHEEPEles_gsfP4ptr_;

		//classification (couldnt they have just named it 'type')
		std::vector<Int_t>*  normHEEPEles_classification_;
		std::vector<Bool_t>* normHEEPEles_isEcalDriven_;
		std::vector<Bool_t>* normHEEPEles_isTrackerDriven_;
		std::vector<Bool_t>* normHEEPEles_isEB_;
		std::vector<Bool_t>* normHEEPEles_isEE_;

		//track methods
		std::vector<Int_t>* normHEEPEles_charge_;
		std::vector<Int_t>* normHEEPEles_trkCharge_;
		std::vector<float>* normHEEPEles_pVtx_;
		std::vector<float>* normHEEPEles_pCalo_;
		std::vector<float>* normHEEPEles_ptVtx_;
		std::vector<float>* normHEEPEles_ptCalo_;

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		std::vector<float>* normHEEPEles_hOverE_;
		std::vector<float>* normHEEPEles_dEtaIn_;
		std::vector<float>* normHEEPEles_dPhiIn_;
		std::vector<float>* normHEEPEles_dPhiOut_;
		std::vector<float>* normHEEPEles_epIn_;
		std::vector<float>* normHEEPEles_epOut_;
		std::vector<float>* normHEEPEles_fbrem_;
		std::vector<float>* normHEEPEles_bremFrac_;
		std::vector<float>* normHEEPEles_invEOverInvP_;

		//shower shape variables
		std::vector<float>* normHEEPEles_sigmaEtaEta_;
		std::vector<float>* normHEEPEles_sigmaEtaEtaUnCorr_;
		std::vector<float>* normHEEPEles_sigmaIEtaIEta_;
		std::vector<float>* normHEEPEles_e1x5_;
		std::vector<float>* normHEEPEles_e2x5Max_;
		std::vector<float>* normHEEPEles_e5x5_;
		std::vector<float>* normHEEPEles_e1x5Over5x5_;
		std::vector<float>* normHEEPEles_e2x5MaxOver5x5_;

		//isolation, we use cone of 0.3
		std::vector<float>* normHEEPEles_isolEm_;
		std::vector<float>* normHEEPEles_isolHad_;
		std::vector<float>* normHEEPEles_isolHadDepth1_;
		std::vector<float>* normHEEPEles_isolHadDepth2_;
		std::vector<float>* normHEEPEles_isolPtTrks_;
		std::vector<float>* normHEEPEles_isolEmHadDepth1_;

		std::vector<tsw::HEEPEle>  normEles_;

		//------------//
		//Variables for storing contents of special boosted reco'n GSF electron branches...
		unsigned int bstdGsfEles_number_;
		std::vector<ROOT::Math::XYZTVector>* bstdGsfEles_OLDp4_;

		std::vector<TLorentzVector> bstdGsfEles_p4_;
		std::vector<Int_t>* bstdGsfEles_charge_;

		std::vector<Double_t>* bstdGsfEles_Et_;
		std::vector<Double_t>* bstdGsfEles_HEEP_Et_;
		std::vector<Double_t>* bstdGsfEles_Eta_;
		std::vector<Double_t>* bstdGsfEles_scEta_;
		std::vector<bool>*     bstdGsfEles_ecalDriven_;
		std::vector<bool>*     bstdGsfEles_ecalDrivenSeed_;

		std::vector<Double_t>* bstdGsfEles_HEEP_dEtaIn_;
		std::vector<Double_t>* bstdGsfEles_HEEP_dPhiIn_;
		std::vector<Double_t>* bstdGsfEles_HoverE_;
		std::vector<Double_t>* bstdGsfEles_sigmaIetaIeta_;
		std::vector<Double_t>* bstdGsfEles_scSigmaIetaIeta_;

		std::vector<Double_t>* bstdGsfEles_dr03EmIsoEt_;
		std::vector<Double_t>* bstdGsfEles_dr03HadDepth1IsoEt_;
		std::vector<Double_t>* bstdGsfEles_dr03HadDepth2IsoEt_;
		std::vector<Double_t>* bstdGsfEles_dr03TkIsoPt_;

		std::vector<Double_t>* bstdGsfEles_e2x5Max_;
		std::vector<Double_t>* bstdGsfEles_e5x5_;

		// heep::Ele method values ...
		//kinematic and geometric methods
		std::vector<float>* bstdHEEPEles_et_;
		std::vector<float>* bstdHEEPEles_gsfEt_;
		std::vector<float>* bstdHEEPEles_scEt_;
		std::vector<float>* bstdHEEPEles_energy_;
		std::vector<float>* bstdHEEPEles_gsfEnergy_;
		std::vector<float>* bstdHEEPEles_caloEnergy_;
		std::vector<float>* bstdHEEPEles_eta_;
		std::vector<float>* bstdHEEPEles_scEta_;
		std::vector<float>* bstdHEEPEles_detEta_;
		std::vector<float>* bstdHEEPEles_detEtaAbs_;
		std::vector<float>* bstdHEEPEles_phi_;
		std::vector<float>* bstdHEEPEles_scPhi_;
		std::vector<float>* bstdHEEPEles_detPhi_;
		std::vector<float>* bstdHEEPEles_zVtx_;
		std::vector<ROOT::Math::XYZTVector>* bstdHEEPEles_p4ptr_;
		std::vector<ROOT::Math::XYZTVector>* bstdHEEPEles_gsfP4ptr_;

		//classification (couldnt they have just named it 'type')
		std::vector<Int_t>*  bstdHEEPEles_classification_;
		std::vector<Bool_t>* bstdHEEPEles_isEcalDriven_;
		std::vector<Bool_t>* bstdHEEPEles_isTrackerDriven_;
		std::vector<Bool_t>* bstdHEEPEles_isEB_;
		std::vector<Bool_t>* bstdHEEPEles_isEE_;

		//track methods
		std::vector<Int_t>* bstdHEEPEles_charge_;
		std::vector<Int_t>* bstdHEEPEles_trkCharge_;
		std::vector<float>* bstdHEEPEles_pVtx_;
		std::vector<float>* bstdHEEPEles_pCalo_;
		std::vector<float>* bstdHEEPEles_ptVtx_;
		std::vector<float>* bstdHEEPEles_ptCalo_;

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		std::vector<float>* bstdHEEPEles_hOverE_;
		std::vector<float>* bstdHEEPEles_dEtaIn_;
		std::vector<float>* bstdHEEPEles_dPhiIn_;
		std::vector<float>* bstdHEEPEles_dPhiOut_;
		std::vector<float>* bstdHEEPEles_epIn_;
		std::vector<float>* bstdHEEPEles_epOut_;
		std::vector<float>* bstdHEEPEles_fbrem_;
		std::vector<float>* bstdHEEPEles_bremFrac_;
		std::vector<float>* bstdHEEPEles_invEOverInvP_;

		//shower shape variables
		std::vector<float>* bstdHEEPEles_sigmaEtaEta_;
		std::vector<float>* bstdHEEPEles_sigmaEtaEtaUnCorr_;
		std::vector<float>* bstdHEEPEles_sigmaIEtaIEta_;
		std::vector<float>* bstdHEEPEles_e1x5_;
		std::vector<float>* bstdHEEPEles_e2x5Max_;
		std::vector<float>* bstdHEEPEles_e5x5_;
		std::vector<float>* bstdHEEPEles_e1x5Over5x5_;
		std::vector<float>* bstdHEEPEles_e2x5MaxOver5x5_;

		//isolation, we use cone of 0.3
		std::vector<float>* bstdHEEPEles_isolEm_;
		std::vector<float>* bstdHEEPEles_isolHad_;
		std::vector<float>* bstdHEEPEles_isolHadDepth1_;
		std::vector<float>* bstdHEEPEles_isolHadDepth2_;
		std::vector<float>* bstdHEEPEles_isolPtTrks_;
		std::vector<float>* bstdHEEPEles_isolEmHadDepth1_;

		std::vector<tsw::HEEPEle> bstdEles_;

		//----------------------------------------------//
		// Muon collections ...
		tsw::MuonCollection normMuons_;
		tsw::MuonCollection normMuons_tight_;

		// Muon histograms ...
		tsw::MuonDistns normMuons_1stpT_Histos_;
		tsw::MuonDistns normMuons_tight_1stpT_Histos_;

		// emu method: emu object & output tree stuff
		tsw::HEEPEle normEles_HEEPNoIso_1stpT_;
		bool normEles_HEEPNoIso_1stpT_exists_;
		tsw::HEEPDiEle normDiEle_HEEPNoIso_;
		bool normDiEle_HEEPNoIso_exists_;
		tsw::EleMuObject eleMu_HEEPNoIso_tight_;

		TTree* eleMuTree_;
		Double_t eleMuTreeVar_mass_;
		Double_t eleMuTreeVar_weight_;

		TTree* diEleTree_;
		Double_t diEleTreeVar_mass_;
		Double_t diEleTreeVar_weight_;

		//----------------------------------------------//
		//Histograms...
		//Special reco'n validation histos:
		tsw::ReconValidationHistos normEles_reconValidationHistos_;
		tsw::ReconValidationHistos normEles_simpleCuts_reconValHistos_;
		tsw::ReconValidationHistos normEles_HEEPCuts_reconValHistos_;
		tsw::ReconValidationHistos bstdEles_reconValidationHistos_;
		tsw::ReconValidationHistos bstdEles_simpleCuts_reconValHistos_;
		tsw::ReconValidationHistos bstdEles_HEEPCuts_reconValHistos_;
		tsw::DiEleDistns normDiEle_AllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_HEEPNoIso_Histos_;
		tsw::DiEleDistns bstdDiEle_AllHEEP_Histos_;
		tsw::DiEleDistns bstdDiEle_M0_AllHEEP_Histos_;
		tsw::DiEleDistns bstdDiEle_Mleq0_AllHEEP_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_AllHEEP_Histos_;
		tsw::DiEleDistns bstdDiEle_HEEPNoIso_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_Histos_;

		// DiEle distributions for various different signal triggger Et thresholds ...
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_5e32_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_1e33_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_2e33_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_TrgEt80_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_TrgEt100_Histos_;
		tsw::DiEleDistns bstdDiEle_diffSC_HEEPNoIso_TrgEt120_Histos_;

		// Continuing previous set ...
		tsw::DiEleDistns normDiEle_AllHEEP_EBEB_Histos_;
		tsw::DiEleDistns bstdDiEle_AllHEEP_EBEB_Histos_;
		tsw::DiEleDistns bstdDiEle_HEEPNoIso_EBEB_Histos_;

		// MC particle distributions ...
		tsw::MCParticleDistns mcZboson_allEvts_Histos_;
		tsw::MCParticleDistns mcZboson_withinCMS_Histos_;
		tsw::MCParticleDistns mcZboson_withinAcc_Histos_;
		tsw::MCParticleDistns mcZboson_withinAcc_EBEB_Histos_;
		tsw::MCParticleDistns mcZboson_withinAcc_EEEE_Histos_;
		tsw::MCParticleDistns mcZboson_stdReconNoIsoEvts_Histos_;
		tsw::MCParticleDistns mcZboson_stdReconNoIsoEvts_EBEB_Histos_;
		tsw::MCParticleDistns mcZboson_stdReconNoIsoEvts_EEEE_Histos_;
		tsw::MCParticleDistns mcZboson_modReconNoIsoEvts_Histos_;
		tsw::MCParticleDistns mcZboson_modReconNoIsoEvts_EBEB_Histos_;
		tsw::MCParticleDistns mcZboson_modReconNoIsoEvts_EEEE_Histos_;

};


//-------------------------------------------------------------------//
//--------------------- Public methods ------------------------------//

BstdZeeFirstAnalyser::BstdZeeFirstAnalyser(int runMode, unsigned int numEvts, bool isMC, const TString& inFileName, const TString& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax):
	normMuons_1stpT_Histos_(            "h_normMuons_1stpT_",       "normal", "",           1, 1, 50, 1000.0, 50, 1000.0),
	normMuons_tight_1stpT_Histos_(      "h_normMuons_tight_1stpT_", "normal", "tight cuts", 2, 1, 50, 1000.0, 50, 1000.0),
	normEles_reconValidationHistos_(    "h_normEles_",      "standard",  "",              1,  false),
	normEles_simpleCuts_reconValHistos_("h_normEles_sCuts_", "standard", ", simple cuts", 1,  false), //was 2
	normEles_HEEPCuts_reconValHistos_(  "h_normEles_HEEP_",  "standard", ", HEEP cuts",   1,  true),  //was 8
	bstdEles_reconValidationHistos_(    "h_bstdEles_",       "special",   "",             4,  false),//was 4
	bstdEles_simpleCuts_reconValHistos_("h_bstdEles_sCuts_", "special",  ", simple cuts", 4, false), //was 38
	bstdEles_HEEPCuts_reconValHistos_(  "h_bstdEles_HEEP_",  "special",  ", HEEP cuts",   4,  true),  //was 6
	normDiEle_AllHEEP_Histos_(          "h_normDiEle_AllHEEP_",  "standard", "all HEEP cuts",     1, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_HEEPNoIso_Histos_(        "h_normDiEle_HEEPNoIso_","standard", "all HEEP cuts",     2, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_AllHEEP_Histos_(          "h_bstdDiEle_AllHEEP_",   "special", "all HEEP cuts",    4, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_M0_AllHEEP_Histos_(       "h_bstdDiEle_M0_AllHEEP_","special", "all HEEP cuts, M!=0", 2, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_Mleq0_AllHEEP_Histos_(       "h_bstdDiEle_Mleq0_AllHEEP_","special", "all HEEP cuts, M==0", 8, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_AllHEEP_Histos_(   "h_bstdDiEle_diffSC_AllHEEP_", "special", "all HEEP cuts, diff SC", 4, 1, nbins_mass, nbins_pt, ptmax), //was 8, 7
	bstdDiEle_HEEPNoIso_Histos_(        "h_bstdDiEle_HEEPNoIso_", "special",  "HEEP cuts no iso", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_", "special", "HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax), //was 2, 2
	bstdDiEle_diffSC_HEEPNoIso_5e32_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_5e32_", "special", "5e32 trigger, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_1e33_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_1e33_", "special", "1e33 trigger, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_2e33_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_2e33_", "special", "2e33 trigger, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_TrgEt80_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_TrgEt80_", "special", "trigger E_{T}^{thr}=80GeVc^{-1}, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_TrgEt100_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_TrgEt100_", "special", "trigger E_{T}^{thr}=100GeVc^{-1}, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_diffSC_HEEPNoIso_TrgEt120_Histos_( "h_bstdDiEle_diffSC_HEEPNoIso_TrgEt120_", "special", "trigger E_{T}^{thr}=120GeVc^{-1}, HEEP cuts no iso diff SC", 6, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_AllHEEP_EBEB_Histos_(          "h_normDiEle_AllHEEP_EBEB_",  "standard", "all HEEP cuts",     1, 7, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_AllHEEP_EBEB_Histos_(          "h_bstdDiEle_AllHEEP_EBEB_",   "special",  "all HEEP cuts",    4, 7, nbins_mass, nbins_pt, ptmax),
	bstdDiEle_HEEPNoIso_EBEB_Histos_(        "h_bstdDiEle_HEEPNoIso_EBEB_", "special",  "HEEP cuts no iso", 6, 7, nbins_mass, nbins_pt, ptmax),
	mcZboson_allEvts_Histos_(        "h_mcZboson_allEvts_",           "Z bosons", "all events",                          1, 5, nbins_pt, ptmax),
	mcZboson_withinCMS_Histos_(      "h_mcZboson_withinCMS_",         "Z bosons", "both eles |#eta|<2.5",                1, 1, nbins_pt, ptmax),
	mcZboson_withinAcc_Histos_(      "h_mcZboson_withinAcc_",         "Z bosons", "both eles in (p_{T}, #eta) acc",      9, 1, nbins_pt, ptmax),
	mcZboson_withinAcc_EBEB_Histos_( "h_mcZboson_withinAcc_EBEB_",    "Z bosons", "both eles in (p_{T}, #eta (EB)) acc", 9, 1, nbins_pt, ptmax),
	mcZboson_withinAcc_EEEE_Histos_( "h_mcZboson_withinAcc_EEEE_",    "Z bosons", "both eles in (p_{T}, #eta (EE)) acc", 9, 1, nbins_pt, ptmax),
	mcZboson_stdReconNoIsoEvts_Histos_(      "h_mcZboson_stdReconNoIsoEvts_",      "Z bosons", "std recon events",       4, 1, nbins_pt, ptmax),
	mcZboson_stdReconNoIsoEvts_EBEB_Histos_( "h_mcZboson_stdReconNoIsoEvts_EBEB_", "Z bosons", "std recon EB-EB events", 4, 1, nbins_pt, ptmax),
	mcZboson_stdReconNoIsoEvts_EEEE_Histos_( "h_mcZboson_stdReconNoIsoEvts_EEEE_", "Z bosons", "std recon EE-EE events", 4, 1, nbins_pt, ptmax),
	mcZboson_modReconNoIsoEvts_Histos_(      "h_mcZboson_modReconNoIsoEvts_",      "Z bosons", "mod recon events",       2, 1, nbins_pt, ptmax),
	mcZboson_modReconNoIsoEvts_EBEB_Histos_( "h_mcZboson_modReconNoIsoEvts_EBEB_", "Z bosons", "mod recon EB-EB events", 2, 1, nbins_pt, ptmax),
	mcZboson_modReconNoIsoEvts_EEEE_Histos_( "h_mcZboson_modReconNoIsoEvts_EEEE_", "Z bosons", "mod recon EE-EE events", 2, 1, nbins_pt, ptmax)
{
//BstdZeeFirstAnalyser::BstdZeeFirstAnalyser(){
	//TODO - initialisation list for efficiency/good practice

	//Initialise the member variables that are arguments of the CTOR...
	runMode_ = runMode;
	numEvts_ = numEvts;
	isMC_ = isMC;
	inputFile_name_ = inFileName;
	outputFile_name_ = outFileName;
	vFlg_ = vFlg;
	
	//Set default (typically clearly 'incorrect') values for non-histo variables...
	SetMemberVariablesToDefaultValues();

	//Initialise the histograms that will be filled in analysis...
	InitialiseReconValidationHistos();

	//Setting up the emu method trees ...
	SetupEMuMethodTrees();
}

BstdZeeFirstAnalyser::~BstdZeeFirstAnalyser(){
	//Clear any heap variables here...
	DeleteReconValidationHistos();
}

void BstdZeeFirstAnalyser::DoAnalysis(const Double_t evtWeight){

	//Opening the datafile(s)
	std::cout << " Opening the ntuple datafile..." << std::endl;
	TFile* inputFile = TFile::Open(inputFile_name_,"READ");
	/*if(!inputFile_ptr)
		std::cout << "ERROR in opening file." << std::endl;
	else
		std::cout << "   file opened successfully." << std::endl;*/

	//Changing into the correct directory within the input file ..
	inputFile->cd("demo");
	//Getting a pointer to the ntuple tree in the file
	inputFile_tree_ = (TTree*)gDirectory->Get("EventDataTree");

	//Setup the branch links...
	SetupBranchLinks(inputFile);

	//TODO - put in code so that a negative value for numEvts_ results in all events being run over.

	//Call AnalyseEvent method for each event...
	for(unsigned int evtIdx=0; evtIdx<numEvts_; evtIdx++){
		if(vFlg_>0){std::cout << std::endl << " Analysing event no. " << evtIdx << std::endl;}
		//Load in data for the evtIdx'th event...
		Long64_t dataTreeEntry = inputFile_tree_->LoadTree(evtIdx);
		inputFile_tree_->GetEntry(dataTreeEntry);

		if( (vFlg_>-2) && (evtIdx%10000==0) ){std::cout << " *** Event no. " << evtIdx << " reached." << std::endl;}
		//TODO - Put setting up of TLorentzVector 4 momenta here ...

		//Analyse this data...
		AnalyseEvent(evtWeight);
	}

	//Close the input file ...
	inputFile->Close();

	FinishOffAnalysis(); //Output file is opened in here ...
}

//-------------------------------------------------------------------//
//---------------------- Private methods ----------------------------//
//-------------------------------------------------------------------//

void BstdZeeFirstAnalyser::AnalyseEvent(const Double_t eventWeight){
	//Setting up the TLorentzVector 4momenta...

	//mcEles_HighestEt_p4.SetPxPyPzE(mcEles_HighestEt_OLDp4_ptr->Px(),mcEles_HighestEt_OLDp4_ptr->Py(),mcEles_HighestEt_OLDp4_ptr->Pz(),mcEles_HighestEt_OLDp4_ptr->E());
	if(isMC_){
		mcEles_HighestEt_p4_ = ConvertToTLorentzVector(mcEles_HighestEt_OLDp4_ptr_);
		mcEles_2ndHighestEt_p4_ = ConvertToTLorentzVector(mcEles_2ndHighestEt_OLDp4_ptr_);
		mcZcandidate_p4_ = ConvertToTLorentzVector(mcZcandidate_OLDp4_ptr_);
		mcZ_daughterA_p4_ = ConvertToTLorentzVector(mcZ_daughterA_OLDp4_ptr_);
		mcZ_daughterB_p4_ = ConvertToTLorentzVector(mcZ_daughterB_OLDp4_ptr_);
	}

	normGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	ROOT::Math::XYZTVector tmpXYZTVector;
	for(unsigned int iEle = 0; iEle<normGsfEles_number_; iEle++){
		tmpXYZTVector = normGsfEles_OLDp4_->at(iEle);
		normGsfEles_p4_.push_back(  ConvertToTLorentzVector(&tmpXYZTVector)  );
	}
	bstdGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	for(unsigned int iEle = 0; iEle<bstdGsfEles_number_; iEle++){
		tmpXYZTVector = bstdGsfEles_OLDp4_->at(iEle);
		bstdGsfEles_p4_.push_back(  ConvertToTLorentzVector(&tmpXYZTVector)  );
	}

	SetupEleClassVectors();

	mcZboson_ = tsw::MCParticle(mcZ_pdgId_, mcZ_status_, -999, *mcZ_p4_ptr_);
	//mcZboson_ = tsw::MCParticle();
	//mcZ_daughterA_p4_.SetPxPyPzE(99999.9, 0.0, 0.0, 99999.9);
	//mcZ_daughterB_p4_.SetPxPyPzE(99999.9, 0.0, 0.0, 99999.9);

	if(vFlg_>0){PrintOutBranchVariables();}
	normEles_ = OrderHEEPElesByEt(normEles_);
	bstdEles_ = OrderHEEPElesByEt(bstdEles_);

	SetupMuonCollection();

	///////////////////////////////////////
	// Now, can actually do some physics!

	FillHistograms(eventWeight);
	DoEMuAnalysis(eventWeight);
}

void BstdZeeFirstAnalyser::SetupEleClassVectors(){
	tsw::EleStruct ithEleStruct;
	SetDefaultValuesInEleStruct(ithEleStruct);
	tsw::HEEPEle ithHEEPEle;

	if(vFlg_>0){std::cout << " Resetting the Ele class vectors ..." << std::endl;}
	normEles_.clear();
	bstdEles_.clear();
	if(vFlg_>0){std::cout << " ... and setting up the new Ele class vectors for this event..." << std::endl;}

	for(unsigned int iEle = 0; iEle < normGsfEles_number_; iEle++){
		if(vFlg_>0){std::cout << "   iEle=" << iEle << std::endl;}
		SetDefaultValuesInEleStruct(ithEleStruct);
		ithHEEPEle = tsw::HEEPEle::HEEPEle(ithEleStruct);

		//kinematic and geometric methods
		ithEleStruct.et_         = normHEEPEles_et_->at(iEle);
		ithEleStruct.gsfEt_      = normHEEPEles_gsfEt_->at(iEle);
		ithEleStruct.scEt_       = normHEEPEles_scEt_->at(iEle);
		ithEleStruct.energy_     = normHEEPEles_energy_->at(iEle);
		ithEleStruct.gsfEnergy_  = normHEEPEles_gsfEnergy_->at(iEle);
		ithEleStruct.caloEnergy_ = normHEEPEles_caloEnergy_->at(iEle);
		ithEleStruct.eta_        = normHEEPEles_eta_->at(iEle);
		ithEleStruct.scEta_      = normHEEPEles_scEta_->at(iEle);
		ithEleStruct.detEta_     = normHEEPEles_detEta_->at(iEle);
		ithEleStruct.detEtaAbs_  = normHEEPEles_detEtaAbs_->at(iEle);
		ithEleStruct.phi_        = normHEEPEles_phi_->at(iEle);
		ithEleStruct.scPhi_      = normHEEPEles_scPhi_->at(iEle);
		ithEleStruct.detPhi_     = normHEEPEles_detPhi_->at(iEle);
		ithEleStruct.zVtx_       = normHEEPEles_zVtx_->at(iEle);
		ithEleStruct.p4_         = normHEEPEles_p4ptr_->at(iEle);
		ithEleStruct.gsfP4_      = normHEEPEles_gsfP4ptr_->at(iEle);

		//classification (couldnt they have just named it 'type')
		ithEleStruct.classification_ =  normHEEPEles_classification_->at(iEle);
		ithEleStruct.isEcalDriven_ = normHEEPEles_isEcalDriven_->at(iEle);
		ithEleStruct.isTrackerDriven_ = normHEEPEles_isTrackerDriven_->at(iEle);
		ithEleStruct.isEB_ = normHEEPEles_isEB_->at(iEle);
		ithEleStruct.isEE_ = normHEEPEles_isEE_->at(iEle);

		//track methods
		ithEleStruct.charge_ = normHEEPEles_charge_->at(iEle);
		ithEleStruct.trkCharge_ = normHEEPEles_trkCharge_->at(iEle);
		ithEleStruct.pVtx_ = normHEEPEles_pVtx_->at(iEle);
		ithEleStruct.pCalo_ = normHEEPEles_pCalo_->at(iEle);
		ithEleStruct.ptVtx_ = normHEEPEles_ptVtx_->at(iEle);
		ithEleStruct.ptCalo_ = normHEEPEles_ptCalo_->at(iEle);

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		ithEleStruct.hOverE_ = normHEEPEles_hOverE_->at(iEle);
		ithEleStruct.dEtaIn_ = normHEEPEles_dEtaIn_->at(iEle);
		ithEleStruct.dPhiIn_ = normHEEPEles_dPhiIn_->at(iEle);
		ithEleStruct.dPhiOut_ = normHEEPEles_dPhiOut_->at(iEle);
		ithEleStruct.epIn_ = normHEEPEles_epIn_->at(iEle);
		ithEleStruct.epOut_ = normHEEPEles_epOut_->at(iEle);
		ithEleStruct.fbrem_ = normHEEPEles_fbrem_->at(iEle);
		ithEleStruct.bremFrac_ = normHEEPEles_bremFrac_->at(iEle);
		ithEleStruct.invEOverInvP_ = normHEEPEles_invEOverInvP_->at(iEle);

		//shower shape variables
		ithEleStruct.sigmaEtaEta_ = normHEEPEles_sigmaEtaEta_->at(iEle);
		ithEleStruct.sigmaEtaEtaUnCorr_ = normHEEPEles_sigmaEtaEtaUnCorr_->at(iEle);
		ithEleStruct.sigmaIEtaIEta_ = normHEEPEles_sigmaIEtaIEta_->at(iEle);
		ithEleStruct.e1x5_ = normHEEPEles_e1x5_->at(iEle);
		ithEleStruct.e2x5Max_ = normHEEPEles_e2x5Max_->at(iEle);
		ithEleStruct.e5x5_ = normHEEPEles_e5x5_->at(iEle);
		ithEleStruct.e1x5Over5x5_ = normHEEPEles_e1x5Over5x5_->at(iEle);
		ithEleStruct.e2x5MaxOver5x5_ = normHEEPEles_e2x5MaxOver5x5_->at(iEle);

		//isolation, we use cone of 0.3
		ithEleStruct.isolEm_ = normHEEPEles_isolEm_->at(iEle);
		ithEleStruct.isolHad_ = normHEEPEles_isolHad_->at(iEle);
		ithEleStruct.isolHadDepth1_ = normHEEPEles_isolHadDepth1_->at(iEle);
		ithEleStruct.isolHadDepth2_ = normHEEPEles_isolHadDepth2_->at(iEle);
		ithEleStruct.isolPtTrks_ = normHEEPEles_isolPtTrks_->at(iEle);
		ithEleStruct.isolEmHadDepth1_ = normHEEPEles_isolEmHadDepth1_->at(iEle);

		ithHEEPEle = tsw::HEEPEle::HEEPEle(ithEleStruct);
		normEles_.push_back(ithHEEPEle);
	}

	for(unsigned int iEle = 0; iEle < bstdGsfEles_number_; iEle++){
		if(vFlg_>0){std::cout << "   iEle=" << iEle << std::endl;}
		SetDefaultValuesInEleStruct(ithEleStruct);
		ithHEEPEle = tsw::HEEPEle::HEEPEle(ithEleStruct);

		//kinematic and geometric methods
		ithEleStruct.et_         = bstdHEEPEles_et_->at(iEle);
		ithEleStruct.gsfEt_      = bstdHEEPEles_gsfEt_->at(iEle);
		ithEleStruct.scEt_       = bstdHEEPEles_scEt_->at(iEle);
		ithEleStruct.energy_     = bstdHEEPEles_energy_->at(iEle);
		ithEleStruct.gsfEnergy_  = bstdHEEPEles_gsfEnergy_->at(iEle);
		ithEleStruct.caloEnergy_ = bstdHEEPEles_caloEnergy_->at(iEle);
		ithEleStruct.eta_        = bstdHEEPEles_eta_->at(iEle);
		ithEleStruct.scEta_      = bstdHEEPEles_scEta_->at(iEle);
		ithEleStruct.detEta_     = bstdHEEPEles_detEta_->at(iEle);
		ithEleStruct.detEtaAbs_  = bstdHEEPEles_detEtaAbs_->at(iEle);
		ithEleStruct.phi_        = bstdHEEPEles_phi_->at(iEle);
		ithEleStruct.scPhi_      = bstdHEEPEles_scPhi_->at(iEle);
		ithEleStruct.detPhi_     = bstdHEEPEles_detPhi_->at(iEle);
		ithEleStruct.zVtx_       = bstdHEEPEles_zVtx_->at(iEle);
		ithEleStruct.p4_         = bstdHEEPEles_p4ptr_->at(iEle);
		ithEleStruct.gsfP4_      = bstdHEEPEles_gsfP4ptr_->at(iEle);

		//classification (couldnt they have just named it 'type')
		ithEleStruct.classification_ =  bstdHEEPEles_classification_->at(iEle);
		ithEleStruct.isEcalDriven_ = bstdHEEPEles_isEcalDriven_->at(iEle);
		ithEleStruct.isTrackerDriven_ = bstdHEEPEles_isTrackerDriven_->at(iEle);
		ithEleStruct.isEB_ = bstdHEEPEles_isEB_->at(iEle);
		ithEleStruct.isEE_ = bstdHEEPEles_isEE_->at(iEle);

		//track methods
		ithEleStruct.charge_ = bstdHEEPEles_charge_->at(iEle);
		ithEleStruct.trkCharge_ = bstdHEEPEles_trkCharge_->at(iEle);
		ithEleStruct.pVtx_ = bstdHEEPEles_pVtx_->at(iEle);
		ithEleStruct.pCalo_ = bstdHEEPEles_pCalo_->at(iEle);
		ithEleStruct.ptVtx_ = bstdHEEPEles_ptVtx_->at(iEle);
		ithEleStruct.ptCalo_ = bstdHEEPEles_ptCalo_->at(iEle);

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		ithEleStruct.hOverE_ = bstdHEEPEles_hOverE_->at(iEle);
		ithEleStruct.dEtaIn_ = bstdHEEPEles_dEtaIn_->at(iEle);
		ithEleStruct.dPhiIn_ = bstdHEEPEles_dPhiIn_->at(iEle);
		ithEleStruct.dPhiOut_ = bstdHEEPEles_dPhiOut_->at(iEle);
		ithEleStruct.epIn_ = bstdHEEPEles_epIn_->at(iEle);
		ithEleStruct.epOut_ = bstdHEEPEles_epOut_->at(iEle);
		ithEleStruct.fbrem_ = bstdHEEPEles_fbrem_->at(iEle);
		ithEleStruct.bremFrac_ = bstdHEEPEles_bremFrac_->at(iEle);
		ithEleStruct.invEOverInvP_ = bstdHEEPEles_invEOverInvP_->at(iEle);

		//shower shape variables
		ithEleStruct.sigmaEtaEta_ = bstdHEEPEles_sigmaEtaEta_->at(iEle);
		ithEleStruct.sigmaEtaEtaUnCorr_ = bstdHEEPEles_sigmaEtaEtaUnCorr_->at(iEle);
		ithEleStruct.sigmaIEtaIEta_ = bstdHEEPEles_sigmaIEtaIEta_->at(iEle);
		ithEleStruct.e1x5_ = bstdHEEPEles_e1x5_->at(iEle);
		ithEleStruct.e2x5Max_ = bstdHEEPEles_e2x5Max_->at(iEle);
		ithEleStruct.e5x5_ = bstdHEEPEles_e5x5_->at(iEle);
		ithEleStruct.e1x5Over5x5_ = bstdHEEPEles_e1x5Over5x5_->at(iEle);
		ithEleStruct.e2x5MaxOver5x5_ = bstdHEEPEles_e2x5MaxOver5x5_->at(iEle);

		//isolation, we use cone of 0.3
		ithEleStruct.isolEm_ = bstdHEEPEles_isolEm_->at(iEle);
		ithEleStruct.isolHad_ = bstdHEEPEles_isolHad_->at(iEle);
		ithEleStruct.isolHadDepth1_ = bstdHEEPEles_isolHadDepth1_->at(iEle);
		ithEleStruct.isolHadDepth2_ = bstdHEEPEles_isolHadDepth2_->at(iEle);
		ithEleStruct.isolPtTrks_ = bstdHEEPEles_isolPtTrks_->at(iEle);
		ithEleStruct.isolEmHadDepth1_ = bstdHEEPEles_isolEmHadDepth1_->at(iEle);

		ithHEEPEle = tsw::HEEPEle::HEEPEle(ithEleStruct);
		bstdEles_.push_back(ithHEEPEle);
	}

	if(vFlg_>0){std::cout << "   done." << std::endl;}


}

void BstdZeeFirstAnalyser::SetupMuonCollection(){
	eventHelper_ = tsw::EventHelper(event_);

	normMuons_ = eventHelper_.GetNormMuons();
	if(vFlg_>0)
		normMuons_.Print();

	// Find two highest pT muons before ordering the collection by pT as a test of these methods ...
	tsw::Muon highestPtMuon;
	if(normMuons_.NumOfMuons()>0)
		highestPtMuon = normMuons_.HighestPtMuon();
	tsw::Muon secondHighestPtMuon;
	if(normMuons_.NumOfMuons()>1)
		secondHighestPtMuon = normMuons_.SecondHighestPtMuon();

	if(vFlg_>0){
		std::cout << "   ***-------------------------------------***" << std::endl;
		std::cout << "   * (From BEFORE ordering muons by pT ... " << std::endl;
		std::cout << "   * Highest pT muon: " << std::endl;
		std::cout << "   *     (pT, eta, phi) = (" << highestPtMuon.pT() << ", " << highestPtMuon.eta() << ", " << highestPtMuon.phi() << std::endl;
		std::cout << "   * 2nd highest pT muon" << std::endl;
		std::cout << "   *     (pT, eta, phi) = (" << secondHighestPtMuon.pT() << ", " << secondHighestPtMuon.eta() << ", " << secondHighestPtMuon.phi() << std::endl;
		std::cout << "   ***-------------------------------------***" << std::endl;
	}
	// Now, ordering the collection by pT, and then finding the two highest pT muons (to check that both methods => same result)...
	normMuons_.OrderBypT();
	if(normMuons_.NumOfMuons()>0)
		highestPtMuon = normMuons_.HighestPtMuon();
	if(normMuons_.NumOfMuons()>1)
		secondHighestPtMuon = normMuons_.SecondHighestPtMuon();

	if(vFlg_>0){
		std::cout << "   ***-------------------------------------***" << std::endl;
		std::cout << "   * (From AFTER ordering muons by pT ... " << std::endl;
		std::cout << "   * Highest pT muon: " << std::endl;
		std::cout << "   *     (pT, eta, phi) = (" << highestPtMuon.pT() << ", " << highestPtMuon.eta() << ", " << highestPtMuon.phi() << std::endl;
		std::cout << "   * 2nd highest pT muon" << std::endl;
		std::cout << "   *     (pT, eta, phi) = (" << secondHighestPtMuon.pT() << ", " << secondHighestPtMuon.eta() << ", " << secondHighestPtMuon.phi() << std::endl;
		std::cout << "   ***-------------------------------------***" << std::endl;
	}

	//Now, applying the tight cuts to the muons, printing out the new collection, ordering it by (descending) pT, and then re-printing out the collection
	if(vFlg_>0){
		std::cout << "   ***---***---***---***---***---***---***---***" << std::endl;
		std::cout << "   ***---***   (TIGHT CUTS APPLIED!)   ***---***" << std::endl;
	}
	normMuons_tight_ = normMuons_.GetTightMuons();

	normMuons_tight_.OrderBypT();
	if(vFlg_>0)
		normMuons_tight_.Print();
}

void BstdZeeFirstAnalyser::PrintOutBranchVariables(){
	//Printing the event information to screen ...
	std::cout << "  ->Evt info (from event branch):"; event_->PrintBasicEventInformation();
	std::cout << "  ->Evt info: run no. " << evt_runNum_ << ", lumi sec. " << evt_lumiSec_<< ", event num. " << evt_evtNum_ << std::endl;
	//Printing the trigger decision to screen ...
	if(trg_PathA_decision_)
		std::cout << "  ->This event passed HLT trigger path A (i.e. " << *trg_PathA_name_ptr_ << ")" << std::endl;
	else
		std::cout << "  ->This event did *NOT* pass HLT trigger path A (i.e. " << *trg_PathA_name_ptr_ << ")" << std::endl;
	std::cout << "      Et of highest Et obj passing this path was: " << trg_PathA_highestTrigObjEt_ << std::endl;
	std::cout << "      (Last filter for this path was: " << *trg_PathA_nameOfLastFilter_ptr_ << ")" << std::endl;

	if(isMC_){
		//Printing the properties of the 2 highest Et MC eles and the resulting Z candidate to screen ...
		std::cout << "  ->*** genParticle information ... ***" << std::endl;
		std::cout << "    ->For the highest Et ele in the event..." << std::endl;
		std::cout << "        charge=" << mcEles_HighestEt_charge_ << "; PDGid=" << mcEles_HighestEt_PDGid_ << "; status=" << mcEles_HighestEt_status_ << std::endl;
		std::cout << "        P=" << mcEles_HighestEt_OLDp4_ptr_->P() << "=" << mcEles_HighestEt_p4_.P() << "?; pt=" << mcEles_HighestEt_pt_ << "=" << mcEles_HighestEt_OLDp4_ptr_->Pt() << "=" << mcEles_HighestEt_p4_.Pt() <<"?" << std::endl;
		std::cout << "        eta=" << mcEles_HighestEt_OLDp4_ptr_->Eta() << "=" << mcEles_HighestEt_p4_.Eta() << "?; phi=" << mcEles_HighestEt_OLDp4_ptr_->Phi() << "=" << mcEles_HighestEt_p4_.Phi() << "?" << std::endl;

		std::cout << "    ->For the 2nd highest Et ele in the event..." << std::endl;
		std::cout << "        charge=" << mcEles_2ndHighestEt_charge_ << "; PDGid=" << mcEles_2ndHighestEt_PDGid_ << "; status=" << mcEles_2ndHighestEt_status_ << std::endl;
		std::cout << "        P=" << mcEles_2ndHighestEt_OLDp4_ptr_->P() << "=" << mcEles_2ndHighestEt_p4_.P() << "?; pt=" << mcEles_2ndHighestEt_pt_ << "=" << mcEles_2ndHighestEt_OLDp4_ptr_->Pt() << "=" << mcEles_2ndHighestEt_p4_.Pt() <<"?" << std::endl;
		std::cout << "        eta=" << mcEles_2ndHighestEt_OLDp4_ptr_->Eta() << "=" << mcEles_2ndHighestEt_p4_.Eta() << "?; phi=" << mcEles_2ndHighestEt_OLDp4_ptr_->Phi() << "=" << mcEles_2ndHighestEt_p4_.Phi() << "?" << std::endl;

		std::cout << "    ->For the Z candidate formed from the highest two Et final state electrons in this event..." << std::endl;
		std::cout << "        P=" << mcZcandidate_OLDp4_ptr_->P() << "=" << mcZcandidate_p4_.P() << "?; Pt=" << mcZcandidate_pt_ << "=" << mcZcandidate_p4_.Pt() << "; inv mass=" << mcZcandidate_p4_.M() << std::endl;
		std::cout << "        eta=" << mcZcandidate_eta_ << "=" << mcZcandidate_p4_.Eta() << "?; phi=" << mcZcandidate_phi_ << "=" << mcZcandidate_p4_.Phi() << std::endl;
		std::cout << "        dEta=" << mcZcandidate_dEtaEles_ << "; dPhi=" << mcZcandidate_dPhiEles_ << "; dR=" << mcZcandidate_dREles_ << std::endl;
		std::cout << "        lab frame opening angle = " << mcZcandidate_openingAngle_ << "=" << mcEles_HighestEt_p4_.Angle(mcEles_2ndHighestEt_p4_.Vect()) << std::endl;

		std::cout << "    ->For the Z ('status 3') Z boson ..." << std::endl;
		std::cout << "        pdgId=" << mcZboson_.pdgId() << "; status=" << mcZboson_.hepMCstatus() << "; charge=" << mcZboson_.charge() << "" << std::endl;
		std::cout << "        (pT, eta, phi) = (" << mcZboson_.pT() << ", " << mcZboson_.eta() << ", " << mcZboson_.phi() << ")" << std::endl;
		std::cout << "        Daughters: number=" << mcZ_numDaughters_ << "; dR=" << mcZ_daughters_dR_ << "; dEta=" << mcZ_daughters_dEta_ << "; dPhi=" << mcZ_daughters_dPhi_ << "; ang sep'n=" << mcZ_daughters_openingAngle_ << std::endl;
		std::cout << "           Daughter A: (pT, eta, phi) = (" << mcZ_daughterA_p4_.Pt() << ", " << mcZ_daughterA_p4_.Eta() << ", " << mcZ_daughterA_p4_.Phi() << ")" << std::endl;
		std::cout << "           Daughter B: (pT, eta, phi) = (" << mcZ_daughterB_p4_.Pt() << ", " << mcZ_daughterB_p4_.Eta() << ", " << mcZ_daughterB_p4_.Phi() << ")" << std::endl;
	}

	std::cout << "  ->*** Reconstructed particle information ... ***" << std::endl;
	//Printing the information from the 'normal' GSF electron branches to screen ...
	std::cout << "  ->There are " << normGsfEles_number_ << "=" << normGsfEles_charge_->size() << "=" << normGsfEles_OLDp4_->size() << "=" << normGsfEles_HEEP_Et_->size() << "=" << normGsfEles_ecalDriven_->size() << "?? 'normal' GSF electrons in this event." << std::endl;
	if(normGsfEles_number_!=normEles_.size()){
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
		std::cout << "         ***  ERROR: normGsfEles_number_ does NOT equal normEles_.size()" << std::endl;
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
	}
	else if(vFlg_>1)
		std::cout << "     normEles_.size()=" << normEles_.size() << std::endl;

	for(unsigned int iEle = 0; iEle < normGsfEles_number_; iEle++){
		std::cout << "    ->Normal GSF ele no. " << iEle << ":";
		std::cout << " Charge = " << normGsfEles_charge_->at(iEle) << std::endl;

		std::cout << "       p4: Px=" << normGsfEles_OLDp4_->at(iEle).Px() << "; Py=" << normGsfEles_OLDp4_->at(iEle).Py() << "; Pz=" << normGsfEles_OLDp4_->at(iEle).Pz() << std::endl;
		std::cout << "           Px=" << normGsfEles_p4_.at(iEle).Px() << "; Py="  << normGsfEles_p4_.at(iEle).Py()  << "; Pz="  << normGsfEles_p4_.at(iEle).Pz() << std::endl;
		std::cout << "           Pt=" << normGsfEles_p4_.at(iEle).Pt() << "; Eta=" << normGsfEles_p4_.at(iEle).Eta() << "; Phi=" << normGsfEles_p4_.at(iEle).Phi() << std::endl;

		std::cout << "       Et=" << normGsfEles_Et_->at(iEle);
		std::cout << "; HEEP_Et=" << normGsfEles_HEEP_Et_->at(iEle) << std::endl;
		std::cout << "       Eta=" << normGsfEles_Eta_->at(iEle);
		std::cout << "; scEta=" << normGsfEles_scEta_->at(iEle) << std::endl;
		std::cout << "       ecalDriven=" << normGsfEles_ecalDriven_->at(iEle);
		std::cout << "; ecalDrivenSeed=" << normGsfEles_ecalDrivenSeed_->at(iEle) << std::endl;

		std::cout << "       HEEP_dEtaIn=" << normGsfEles_HEEP_dEtaIn_->at(iEle);
		std::cout << "; HEEP_dPhiIn=" << normGsfEles_HEEP_dPhiIn_->at(iEle) << std::endl;
		std::cout << "       HoverE=" << normGsfEles_HoverE_->at(iEle);
		std::cout << "; sigmaIetaIeta=" << normGsfEles_sigmaIetaIeta_->at(iEle);
		std::cout << "; scSigmaIetaIeta=" << normGsfEles_scSigmaIetaIeta_->at(iEle) << std::endl;

		std::cout << "       dr03EmIsoEt=" << normGsfEles_dr03EmIsoEt_->at(iEle);
		std::cout << "; dr03HadDepth1IsoEt=" << normGsfEles_dr03HadDepth1IsoEt_->at(iEle) << std::endl;
		std::cout << "       dr03HadDepth2IsoEt=" << normGsfEles_dr03HadDepth2IsoEt_->at(iEle);
		std::cout << "; dr03TkIsoPt=" << normGsfEles_dr03TkIsoPt_->at(iEle) << std::endl;

		std::cout << "       e2x5Max=" << normGsfEles_e2x5Max_->at(iEle);
		std::cout << "; e5x5=" << normGsfEles_e5x5_->at(iEle) << std::endl;

		std::cout << "     ->HEEP variables ..." << std::endl;
		normEles_.at(iEle).PrintOutVariables();
	}

	//Printing the information from the special boosted Z(ee) reconstruction GSF electron branches to screen ...
	std::cout << "  ->There are " << bstdGsfEles_number_ << "=" << bstdGsfEles_charge_->size() << "=" << bstdGsfEles_OLDp4_->size() << "=" << bstdGsfEles_HEEP_Et_->size() << "=" << bstdGsfEles_ecalDriven_->size() << " 'special reconstruction' boosted electrons in this event." << std::endl;

	if(bstdGsfEles_number_!=bstdEles_.size()){
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
		std::cout << "         ***  ERROR: bstdGsfEles_number_ does NOT equal bstdEles_.size()" << std::endl;
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
	}
	else if(vFlg_>1)
		std::cout << "     bstdEles_.size()=" << bstdEles_.size() << std::endl;

	for(unsigned int iEle = 0; iEle < bstdGsfEles_number_; iEle++){
		std::cout << "    ->Special reco'n GSF ele no. " << iEle << ":";
		std::cout << " Charge = " << bstdGsfEles_charge_->at(iEle) << std::endl;

		std::cout << "       p4: Px=" << bstdGsfEles_OLDp4_->at(iEle).Px() << "; Py=" << bstdGsfEles_OLDp4_->at(iEle).Py() << "; Pz=" << bstdGsfEles_OLDp4_->at(iEle).Pz() << std::endl;
		std::cout << "           Px=" << bstdGsfEles_p4_.at(iEle).Px() << "; Py="  << bstdGsfEles_p4_.at(iEle).Py()  << "; Pz="  << bstdGsfEles_p4_.at(iEle).Pz() << std::endl;
		std::cout << "           Pt=" << bstdGsfEles_p4_.at(iEle).Pt() << "; Eta=" << bstdGsfEles_p4_.at(iEle).Eta() << "; Phi=" << bstdGsfEles_p4_.at(iEle).Phi() << std::endl;

		std::cout << "       Et=" << bstdGsfEles_Et_->at(iEle);
		std::cout << "; HEEP_Et=" << bstdGsfEles_HEEP_Et_->at(iEle) << std::endl;
		std::cout << "       Eta=" << bstdGsfEles_Eta_->at(iEle);
		std::cout << "; scEta=" << bstdGsfEles_scEta_->at(iEle) << std::endl;
		std::cout << "       ecalDriven=" << bstdGsfEles_ecalDriven_->at(iEle);
		std::cout << "; ecalDrivenSeed=" << bstdGsfEles_ecalDrivenSeed_->at(iEle) << std::endl;

		std::cout << "       HEEP_dEtaIn=" << bstdGsfEles_HEEP_dEtaIn_->at(iEle);
		std::cout << "; HEEP_dPhiIn=" << bstdGsfEles_HEEP_dPhiIn_->at(iEle) << std::endl;
		std::cout << "       HoverE=" << bstdGsfEles_HoverE_->at(iEle);
		std::cout << "; sigmaIetaIeta=" << bstdGsfEles_sigmaIetaIeta_->at(iEle);
		std::cout << "; scSigmaIetaIeta=" << bstdGsfEles_scSigmaIetaIeta_->at(iEle) << std::endl;

		std::cout << "       dr03EmIsoEt=" << bstdGsfEles_dr03EmIsoEt_->at(iEle);
		std::cout << "; dr03HadDepth1IsoEt=" << bstdGsfEles_dr03HadDepth1IsoEt_->at(iEle) << std::endl;
		std::cout << "       dr03HadDepth2IsoEt=" << bstdGsfEles_dr03HadDepth2IsoEt_->at(iEle);
		std::cout << "; dr03TkIsoPt=" << bstdGsfEles_dr03TkIsoPt_->at(iEle) << std::endl;

		std::cout << "       e2x5Max=" << bstdGsfEles_e2x5Max_->at(iEle);
		std::cout << "; e5x5=" << bstdGsfEles_e5x5_->at(iEle) << std::endl;

		std::cout << "     ->HEEP variables ..." << std::endl;
		bstdEles_.at(iEle).PrintOutVariables();
	}
}

void BstdZeeFirstAnalyser::FinishOffAnalysis(){

	//Calculate errors on histogram bins ...



	//Open the output file ...
	TFile f_histos(outputFile_name_+".root","RECREATE");
	f_histos.Write();

	//Save all of the histograms from analysis ...
	SaveReconValidationHistos(&f_histos);

	//Close the output file ...
	f_histos.Close();

	// Open up eMu method files, save the appropriate tree, and then
	TFile f_eleMuTree("eleMuTree/" + outputFile_name_ + "_sample.root","RECREATE");
	eleMuTree_->Write();
	f_eleMuTree.Close();
	delete eleMuTree_;

	TFile f_diEleTree("diEleTree/" + outputFile_name_ + "_sample.root","RECREATE");
	diEleTree_->Write();
	f_diEleTree.Close();
	delete diEleTree_;
}


//Methods for ordering reconstructed objects in terms of kinematic variables...
std::vector<Int_t> BstdZeeFirstAnalyser::OrderByEt(const std::vector<Double_t>& etValues){
	//Currently a placeholder
	std::vector<Int_t> tmpAnswer;
	tmpAnswer.clear(); tmpAnswer.push_back(-999);
	return tmpAnswer;
}
void BstdZeeFirstAnalyser::TwoHighestEtObjects(const std::vector<Double_t>& etValues, int& idx_HighestEtObject, int& idx_2ndHighestEtObject){
	//Currently a placeholder ...
}

//Methods for application of HEEP cuts ...
//
std::vector<bool> BstdZeeFirstAnalyser::HEEPCuts(int eleType){ //Apply to normGSFEles if eleType=0, bstdGSFEles if eleType=1, and return empty vector (and error message) otherwise
	//Currently a placeholder ...
	std::vector<bool> tmpAnswer;
	tmpAnswer.clear();
	return tmpAnswer;
}
std::vector<bool> BstdZeeFirstAnalyser::HEEPCutsWithoutIso(int eleType){
	//Currently a placeholder ...
	std::vector<bool> tmpAnswer;
	tmpAnswer.clear();
	return tmpAnswer;
}

//--------------------------------------------------------//
//---- Methods for filling sets of histograms ...
//--------------------------------------------------------//
void BstdZeeFirstAnalyser::FillReconValidationHistos(const Double_t weight){

	//hist_normEles_number_->Fill(normGsfEles_number_, weight);
	/*for(unsigned int eleIdx=0; eleIdx<normGsfEles_number_; eleIdx++){
		hist_normEles_charge_->Fill(normGsfEles_charge_->at(eleIdx), weight);
		hist_normEles_Et_->Fill(    normGsfEles_Et_->at(eleIdx),     weight);
		hist_normEles_heepEt_->Fill(normGsfEles_HEEP_Et_, weight);
		hist_normEles_Eta_->Fill(normGsfEles_Eta_->at(eleIdx), weight);
		hist_normEles_scEta_->Fill(normGsfEles_scEta_->at(eleIdx), weight);
		hist_normEles_ecalDriven_->Fill(normGsfEles_ecalDriven_->at(eleIdx), weight);
		hist_normEles_ecalDrivenSeed_->Fill(normGsfEles_ecalDrivenSeed_->at(eleIdx), weight);

		hist_normEles_heepdEtaIn_->Fill(normGsfEles_HEEP_dEtaIn_->at(eleIdx), weight);
		hist_normEles_heepdPhiIn_->Fill(normGsfEles_HEEP_dPhiIn_->at(eleIdx), weight);
		hist_normEles_HoverE_->Fill(normGsfEles_HoverE_->at(eleIdx), weight);
		hist_normEles_sigmaIetaIeta_->Fill(normGsfEles_sigmaIetaIeta_->at(eleIdx), weight);
		hist_normEles_scSigmaIetaIeta_->Fill(normGsfEles_scSigmaIetaIeta_->at(eleIdx), weight);

		hist_normEles_dr03EmIsoEt_  ->Fill(normGsfEles_dr03EmIsoEt_->at(eleIdx), weight);
		hist_normEles_dr03Had1IsoEt_->Fill(normGsfEles_dr03HadDepth1IsoEt_->at(eleIdx), weight);
		hist_normEles_dr03Had2IsoEt_->Fill(normGsfEles_dr03HadDepth2IsoEt_->at(eleIdx), weight);
		hist_normEles_dr03TkIsoPt_  ->Fill(normGsfEles_dr03TkIsoPt_->at(eleIdx), weight);
		hist_normEles_e2x5Max_->Fill(      normGsfEles_e2x5Max_->at(eleIdx), weight);
		hist_normEles_e5x5_->Fill(         normGsfEles_e5x5_->at(eleIdx), weight);
	}*/


}
void BstdZeeFirstAnalyser::FillDiEleHistos(const Double_t weight){
	//Currently a placeholder ...
}
void BstdZeeFirstAnalyser::FillHistograms(const Double_t evtWeight){

	normDiEle_HEEPNoIso_exists_ = false;
	normEles_HEEPNoIso_1stpT_exists_ = false;

	//Currently a placeholder ...
	//FillDiEleHistos(evtWeight);
	unsigned int idx_highestEtEle = 0;
	float highestEleEt = -999.9;
	std::vector<bool> normEles_sCutsFlags;         normEles_sCutsFlags.clear();
	std::vector<bool> normEles_HEEPCutsFlags;      normEles_HEEPCutsFlags.clear();
	std::vector<bool> normEles_HEEPCutsNoIsoFlags; normEles_HEEPCutsNoIsoFlags.clear();
	std::vector<bool> bstdEles_sCutsFlags;         bstdEles_sCutsFlags.clear();
	std::vector<bool> bstdEles_HEEPCutsFlags;      bstdEles_HEEPCutsFlags.clear();
	std::vector<bool> bstdEles_HEEPCutsNoIsoFlags; bstdEles_HEEPCutsNoIsoFlags.clear();

	mcZboson_allEvts_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
	if( fabs(mcZ_daughterA_p4_.Eta())<=2.5 && fabs(mcZ_daughterB_p4_.Eta())<=2.5 ){
		mcZboson_withinCMS_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);

		bool bothElesWithinEtaAcceptance = true;
		bothElesWithinEtaAcceptance = bothElesWithinEtaAcceptance && ( fabs(mcZ_daughterA_p4_.Eta())>=1.56 || fabs(mcZ_daughterA_p4_.Eta())<=1.442 );
		bothElesWithinEtaAcceptance = bothElesWithinEtaAcceptance && ( fabs(mcZ_daughterB_p4_.Eta())>=1.56 || fabs(mcZ_daughterB_p4_.Eta())<=1.442 );

		if( (mcZ_daughterA_p4_.Pt()>=25.0 && mcZ_daughterB_p4_.Pt()>=25.0) && bothElesWithinEtaAcceptance ){
			mcZboson_withinAcc_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			if( fabs(mcZ_daughterA_p4_.Eta())<=1.5 && fabs(mcZ_daughterB_p4_.Eta())<=1.5 )
				mcZboson_withinAcc_EBEB_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			else if( fabs(mcZ_daughterA_p4_.Eta())>=1.5 && fabs(mcZ_daughterB_p4_.Eta())>=1.5 )
				mcZboson_withinAcc_EEEE_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
		}
	}

	normEles_reconValidationHistos_.FillNumberElesHist(normGsfEles_number_ , evtWeight);
	for(unsigned int iEle=0; iEle<normGsfEles_number_; iEle++){
//		normEles_reconValidationHistos_.FillChargeHist(normGsfEles_charge_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillEtHists(normGsfEles_Et_->at(iEle), normGsfEles_HEEP_Et_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillEtaHists(normGsfEles_Eta_->at(iEle), normGsfEles_scEta_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillEcalDrivenHists(normGsfEles_ecalDriven_->at(iEle), normGsfEles_ecalDrivenSeed_->at(iEle), weight);
//		normEles_reconValidationHistos_.FilldEtadPhiInHists(normGsfEles_HEEP_dEtaIn_->at(iEle), normGsfEles_HEEP_dPhiIn_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillHoverEHist(normGsfEles_HoverE_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillSigmaIetaIetaHists(normGsfEles_sigmaIetaIeta_->at(iEle), normGsfEles_scSigmaIetaIeta_->at(iEle), weight);
//		normEles_reconValidationHistos_.FillIsoHists(normGsfEles_dr03EmIsoEt_->at(iEle), normGsfEles_dr03HadDepth1IsoEt_->at(iEle), normGsfEles_dr03HadDepth2IsoEt_->at(iEle), normGsfEles_dr03TkIsoPt_->at(iEle), weight);
//		normEles_reconValidationHistos_.FilleAxBHists(normGsfEles_e2x5Max_->at(iEle), normGsfEles_e5x5_->at(iEle), weight);
		normEles_reconValidationHistos_.FillHistos( normEles_.at(iEle), evtWeight);
		if ( normEles_.at(iEle).ApplySimpleCuts() ){
			normEles_simpleCuts_reconValHistos_.FillHistos( normEles_.at(iEle), evtWeight);}
		//if ( normEles_.at(iEle).ApplyAllHEEPCuts() ){
		//	normEles_HEEPCuts_reconValHistos_.FillHistos( normEles_.at(iEle) );}

		normEles_sCutsFlags.push_back(   normEles_.at(iEle).ApplySimpleCuts());
		normEles_HEEPCutsFlags.push_back(normEles_.at(iEle).ApplyAllHEEPCuts()  );
		normEles_HEEPCutsNoIsoFlags.push_back(normEles_.at(iEle).ApplyHEEPCutsNoIso()  );

		//if(normEles_.at(iEle).et()>highestEleEt){
		//	highestEleEt = normEles_.at(iEle).et();
		//	idx_highestEtEle = iEle;
		//}
	}
	normEles_HEEPCuts_reconValHistos_.FillHistos(normEles_, normEles_HEEPCutsFlags, evtWeight);
	//if(normGsfEles_number_>0)
	//	normEles_reconValidationHistos_.FillHighestEtEleHistos( normEles_.at(idx_highestEtEle) );

	idx_highestEtEle = 0;
	highestEleEt = -999.9;

	bstdEles_reconValidationHistos_.FillNumberElesHist(bstdGsfEles_number_ , evtWeight);
	for(unsigned int iEle=0; iEle<bstdGsfEles_number_; iEle++){
//		bstdEles_reconValidationHistos_.FillChargeHist(bstdGsfEles_charge_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillEtHists(bstdGsfEles_Et_->at(iEle), bstdGsfEles_HEEP_Et_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillEtaHists(bstdGsfEles_Eta_->at(iEle), bstdGsfEles_scEta_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillEcalDrivenHists(bstdGsfEles_ecalDriven_->at(iEle), bstdGsfEles_ecalDrivenSeed_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FilldEtadPhiInHists(bstdGsfEles_HEEP_dEtaIn_->at(iEle), bstdGsfEles_HEEP_dPhiIn_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillHoverEHist(bstdGsfEles_HoverE_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillSigmaIetaIetaHists(bstdGsfEles_sigmaIetaIeta_->at(iEle), bstdGsfEles_scSigmaIetaIeta_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FillIsoHists(bstdGsfEles_dr03EmIsoEt_->at(iEle), bstdGsfEles_dr03HadDepth1IsoEt_->at(iEle), bstdGsfEles_dr03HadDepth2IsoEt_->at(iEle), bstdGsfEles_dr03TkIsoPt_->at(iEle), weight);
//		bstdEles_reconValidationHistos_.FilleAxBHists(bstdGsfEles_e2x5Max_->at(iEle), bstdGsfEles_e5x5_->at(iEle), weight);
		bstdEles_reconValidationHistos_.FillHistos( bstdEles_.at(iEle), evtWeight);
		if ( bstdEles_.at(iEle).ApplySimpleCuts() ){
			bstdEles_simpleCuts_reconValHistos_.FillHistos( bstdEles_.at(iEle), evtWeight);}
		//if ( bstdEles_.at(iEle).ApplyAllHEEPCuts() ){
		//	bstdEles_HEEPCuts_reconValHistos_.FillHistos( bstdEles_.at(iEle) );}

		bstdEles_sCutsFlags.push_back(         bstdEles_.at(iEle).ApplySimpleCuts()    );
		bstdEles_HEEPCutsFlags.push_back(      bstdEles_.at(iEle).ApplyAllHEEPCuts()      );
		bstdEles_HEEPCutsNoIsoFlags.push_back( bstdEles_.at(iEle).ApplyHEEPCutsNoIso() );

		//if(bstdEles_.at(iEle).et()>highestEleEt){
		//	highestEleEt = bstdEles_.at(iEle).et();
		//	idx_highestEtEle = iEle;
		//}
	}
	//if(bstdGsfEles_number_>0)
	//	bstdEles_reconValidationHistos_.FillHighestEtEleHistos( bstdEles_.at(idx_highestEtEle) );
	bstdEles_HEEPCuts_reconValHistos_.FillHistos(bstdEles_, bstdEles_HEEPCutsFlags, evtWeight);


	if(tsw::NumPassingCuts(normEles_HEEPCutsNoIsoFlags)>0){
		std::vector<unsigned int> idxs_elesPassedCuts = IndicesOfElesPassingCuts(normEles_, normEles_HEEPCutsNoIsoFlags, 0);
		normEles_HEEPNoIso_1stpT_ = normEles_.at( idxs_elesPassedCuts.at(0) );
		normEles_HEEPNoIso_1stpT_exists_ = true;
	}

	/////////////////////////////////////////////////
	// Form the di-electron pairs ...
	tsw::HEEPDiEle normDiEle_AllHEEP;
	tsw::HEEPDiEle normDiEle_HEEPNoIso;
	tsw::HEEPDiEle bstdDiEle_AllHEEP;
	tsw::HEEPDiEle bstdDiEle_HEEPNoIso;

	if(tsw::NumPassingCuts(normEles_HEEPCutsFlags)>1){
		normDiEle_AllHEEP = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_HEEPCutsFlags );
		normDiEle_AllHEEP_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);
		if( normDiEle_AllHEEP.eleA().isEB() && normDiEle_AllHEEP.eleB().isEB() )
			normDiEle_AllHEEP_EBEB_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);
	}

	if(tsw::NumPassingCuts(normEles_HEEPCutsNoIsoFlags)>1){
		// Set some of the emu method diele tree branch variablles here ...
		normDiEle_HEEPNoIso = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_HEEPCutsNoIsoFlags );
		normDiEle_HEEPNoIso_ = normDiEle_HEEPNoIso;
		normDiEle_HEEPNoIso_exists_ =  true;

		if(normDiEle_HEEPNoIso.isInZMassRange()){
			normDiEle_HEEPNoIso_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);

			mcZboson_stdReconNoIsoEvts_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);

			if( fabs(mcZ_daughterA_p4_.Eta())<=1.5 && fabs(mcZ_daughterB_p4_.Eta())<=1.5 )
				mcZboson_stdReconNoIsoEvts_EBEB_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			else if( fabs(mcZ_daughterA_p4_.Eta())>=1.56 && fabs(mcZ_daughterB_p4_.Eta())>=1.56 )
				mcZboson_stdReconNoIsoEvts_EEEE_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
		}
	}

	//Modfied reconstruction di-electron pairs ...

	if(tsw::NumPassingCuts(bstdEles_HEEPCutsFlags)>1){
		bstdDiEle_AllHEEP = tsw::HEEPDiEle::HEEPDiEle( bstdEles_, bstdEles_HEEPCutsFlags );
		bstdDiEle_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
		if( bstdDiEle_AllHEEP.invMass()>0.1 ){
			bstdDiEle_M0_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);}
		else{
			bstdDiEle_Mleq0_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);}
		if( fabs(bstdDiEle_AllHEEP.scDeltaEta())>0.001 || fabs(bstdDiEle_AllHEEP.scDeltaPhi()>0.001) )
			bstdDiEle_diffSC_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
		if( bstdDiEle_AllHEEP.eleA().isEB() && bstdDiEle_AllHEEP.eleB().isEB() )
			bstdDiEle_AllHEEP_EBEB_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
	}

	if(tsw::NumPassingCuts(bstdEles_HEEPCutsNoIsoFlags)>1){
		bstdDiEle_HEEPNoIso = tsw::HEEPDiEle::HEEPDiEle( bstdEles_, bstdEles_HEEPCutsNoIsoFlags );
		bstdDiEle_HEEPNoIso_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
		if( bstdDiEle_HEEPNoIso.eleA().isEB() && bstdDiEle_HEEPNoIso.eleB().isEB() )
			bstdDiEle_HEEPNoIso_EBEB_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
		if( ( fabs(bstdDiEle_HEEPNoIso.scDeltaEta())>0.001 || fabs(bstdDiEle_HEEPNoIso.scDeltaPhi()>0.001) ) && bstdDiEle_HEEPNoIso.isInZMassRange() ){
			bstdDiEle_diffSC_HEEPNoIso_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);

			mcZboson_modReconNoIsoEvts_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			if( fabs(mcZ_daughterA_p4_.Eta())<=1.442 && fabs(mcZ_daughterB_p4_.Eta())<=1.442 )
				mcZboson_modReconNoIsoEvts_EBEB_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			else if( fabs(mcZ_daughterA_p4_.Eta())>=1.56 && fabs(mcZ_daughterB_p4_.Eta())>=1.56 )
				mcZboson_modReconNoIsoEvts_EEEE_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);

			if(trg_PathA_decision_){
				bstdDiEle_diffSC_HEEPNoIso_5e32_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
				if(trg_PathA_highestTrigObjEt_ > 52.0)
					bstdDiEle_diffSC_HEEPNoIso_1e33_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
				if(trg_PathA_highestTrigObjEt_ > 65.0)
					bstdDiEle_diffSC_HEEPNoIso_2e33_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
				if(trg_PathA_highestTrigObjEt_ > 80.0)
					bstdDiEle_diffSC_HEEPNoIso_TrgEt80_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
				if(trg_PathA_highestTrigObjEt_ > 100.0)
					bstdDiEle_diffSC_HEEPNoIso_TrgEt100_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
				if(trg_PathA_highestTrigObjEt_ > 120.0)
					bstdDiEle_diffSC_HEEPNoIso_TrgEt120_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);

			}
		}
	}
}

//--------------------------------------------------------//
//---- Method for implementing emu analysis ...           //
void BstdZeeFirstAnalyser::DoEMuAnalysis(const Double_t evtWeight){
	//  Filling histograms with the properties of the highest pT muon (no cuts), and the highest pT tight Muon ...
	if(normMuons_.NumOfMuons()>0)
		normMuons_1stpT_Histos_.FillHistos( normMuons_.HighestPtMuon(), 1.0);
	if(normMuons_tight_.NumOfMuons()>0)
		normMuons_tight_1stpT_Histos_.FillHistos( normMuons_tight_.HighestPtMuon(), 1.0);

	// Forming the emu object (IFF there is at least 1 HEEP non-isol electron, and at least one tight muon) ...
	if(normMuons_tight_.NumOfMuons()>0 && normEles_HEEPNoIso_1stpT_exists_){
		eleMu_HEEPNoIso_tight_ = tsw::EleMuObject(normEles_HEEPNoIso_1stpT_, normMuons_tight_.HighestPtMuon());
		//Now, setting the bracnch variables, and adding a new event to the tree ...
		eleMuTreeVar_mass_   = eleMu_HEEPNoIso_tight_.mass();
		eleMuTreeVar_weight_ = evtWeight;
		eleMuTree_->Fill();
	}

	//Settting the branch variables for the diEle tree - IFF a di-ele has been constructed ...
	if(normDiEle_HEEPNoIso_exists_){
		diEleTreeVar_mass_ = normDiEle_HEEPNoIso_.invMass();
		diEleTreeVar_weight_ = evtWeight;
		diEleTree_->Fill();
	}
}

//--------------------------------------------------------//
//---- Methods for storing histograms in output file...   //
//--------------------------------------------------------//
void BstdZeeFirstAnalyser::SaveReconValidationHistos(TFile* histosFile){
	/*hist_normEles_number_->Write();
	hist_normEles_number_->Write();
	hist_normEles_charge_->Write();
	hist_normEles_Et_->Write();
	hist_normEles_heepEt_->Write();
	hist_normEles_Eta_->Write();
	hist_normEles_scEta_->Write();
	hist_normEles_ecalDriven_->Write();
	hist_normEles_ecalDrivenSeed_->Write();

	hist_normEles_heepdEtaIn_->Write();
	hist_normEles_heepdPhiIn_->Write();
	hist_normEles_HoverE_->Write();
	hist_normEles_sigmaIetaIeta_->Write();
	hist_normEles_scSigmaIetaIeta_->Write();

	hist_normEles_dr03EmIsoEt_->Write();
	hist_normEles_dr03Had1IsoEt_->Write();
	hist_normEles_dr03Had2IsoEt_->Write();
	hist_normEles_dr03TkIsoPt_->Write();
	hist_normEles_e2x5Max_->Write();
	hist_normEles_e5x5_->Write();*/

	normEles_reconValidationHistos_.WriteHistos(histosFile);
	normEles_simpleCuts_reconValHistos_.WriteHistos(histosFile);
	normEles_HEEPCuts_reconValHistos_.WriteHistos(histosFile);
	bstdEles_reconValidationHistos_.WriteHistos(histosFile);
	bstdEles_simpleCuts_reconValHistos_.WriteHistos(histosFile);
	bstdEles_HEEPCuts_reconValHistos_.WriteHistos(histosFile);

	normDiEle_AllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_Histos_.WriteHistos(histosFile);
	bstdDiEle_AllHEEP_Histos_.WriteHistos(histosFile);
	bstdDiEle_HEEPNoIso_Histos_.WriteHistos(histosFile);

	normDiEle_AllHEEP_EBEB_Histos_.WriteHistos(histosFile);
	bstdDiEle_AllHEEP_EBEB_Histos_.WriteHistos(histosFile);
	bstdDiEle_HEEPNoIso_EBEB_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_5e32_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_1e33_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_2e33_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_TrgEt80_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_TrgEt100_Histos_.WriteHistos(histosFile);
	bstdDiEle_diffSC_HEEPNoIso_TrgEt120_Histos_.WriteHistos(histosFile);

	bstdDiEle_diffSC_AllHEEP_Histos_.WriteHistos(histosFile);
	bstdDiEle_M0_AllHEEP_Histos_.WriteHistos(histosFile);
	bstdDiEle_Mleq0_AllHEEP_Histos_.WriteHistos(histosFile);

	mcZboson_allEvts_Histos_.WriteHistos(histosFile);
	mcZboson_withinCMS_Histos_.WriteHistos(histosFile);
	mcZboson_withinAcc_Histos_.WriteHistos(histosFile);
	mcZboson_withinAcc_EBEB_Histos_.WriteHistos(histosFile);
	mcZboson_withinAcc_EEEE_Histos_.WriteHistos(histosFile);
	mcZboson_stdReconNoIsoEvts_Histos_.WriteHistos(histosFile);
	mcZboson_stdReconNoIsoEvts_EBEB_Histos_.WriteHistos(histosFile);
	mcZboson_stdReconNoIsoEvts_EEEE_Histos_.WriteHistos(histosFile);
	mcZboson_modReconNoIsoEvts_Histos_.WriteHistos(histosFile);
	mcZboson_modReconNoIsoEvts_EBEB_Histos_.WriteHistos(histosFile);
	mcZboson_modReconNoIsoEvts_EEEE_Histos_.WriteHistos(histosFile);
}

//----------------------------------------------------------------//
//Methods called by the constructor (e.g. Initialise methods)...
//----------------------------------------------------------------//
void BstdZeeFirstAnalyser::SetMemberVariablesToDefaultValues(){
	std::vector<Int_t> defaultVectInts;
	defaultVectInts.clear(); defaultVectInts.push_back(-999);

	std::vector<Double_t> defaultVectDoubs;
	defaultVectDoubs.clear(); defaultVectDoubs.push_back(-999.9);

	//Variables for storing contents of event information branches...
	evt_runNum_ = 0;
	evt_lumiSec_= 0;
	evt_evtNum_ = 0;

	//Variables for storing contents of trigger branches ...
	trg_PathA_decision_ = false;
	trg_PathA_name_ptr_ = 0;
	trg_PathA_nameOfLastFilter_ptr_ = 0;
	trg_PathA_highestTrigObjEt_ = -999.9;


	//Variables for storing contents of MC final state particle branches...
	mc_numFinalStateEles_ = 0; //unsigned int

   mcEles_HighestEt_charge_ = -999; //Int_t
	mcEles_HighestEt_PDGid_  = -999; //Int_t
	mcEles_HighestEt_status_ = -999; //Int_t
	mcEles_HighestEt_pt_   = -999.9; //Double_t
	mcEles_HighestEt_OLDp4_ptr_ = 0; //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
	mcEles_HighestEt_p4_ = TLorentzVector(0.0,0.0,0.0,0.0); //TLorentzVector

	mcEles_2ndHighestEt_charge_ = -999; //Int_t
	mcEles_2ndHighestEt_PDGid_  = -999; //Int_t
	mcEles_2ndHighestEt_status_ = -999; //Int_t
	mcEles_2ndHighestEt_pt_   = -999.9; //Double_t
	mcEles_2ndHighestEt_OLDp4_ptr_ = 0; //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
	mcEles_2ndHighestEt_p4_ = TLorentzVector(0.0,0.0,0.0,0.0); //TLorentzVector

	mcZcandidate_pt_   = -999.9; //Double_t
	mcZcandidate_eta_  = -999.9; //Double_t
	mcZcandidate_phi_  = -999.9; //Double_t
	mcZcandidate_mass_ = -999.9; //Double_t
	mcZcandidate_OLDp4_ptr_ = 0; //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* //ROOT::Math::XYZTVector
										 //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
	mcZcandidate_p4_ = TLorentzVector(0.0,0.0,0.0,0.0); //TLorentzVector
	mcZcandidate_dEtaEles_ = -999.9; //Double_t
	mcZcandidate_dPhiEles_ = -999.9; //Double_t
	mcZcandidate_dREles_   = -999.9; //Double_t
	mcZcandidate_openingAngle_ = -999.9; //Double_t

	mcZ_pdgId_  = -999;
	mcZ_status_ = -999;
	mcZ_p4_ptr_ = 0; //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
	mcZ_numDaughters_ = 0;
	mcZ_daughters_dR_   = -999.9;
	mcZ_daughters_dEta_ = -999.9;
	mcZ_daughters_dPhi_ = -999.9;
	mcZ_daughters_openingAngle_ = -999.9;
	mcZ_daughterA_OLDp4_ptr_ = 0;
	mcZ_daughterA_p4_ = TLorentzVector(0.0,0.0,0.0,0.0);
	mcZ_daughterB_OLDp4_ptr_ = 0;
	mcZ_daughterB_p4_ = TLorentzVector(0.0,0.0,0.0,0.0);

	//Variables for storing contents of 'normal' GSF electron branches ...
	normGsfEles_number_ = 0;

	normGsfEles_OLDp4_  = 0; // std::vector<ROOT::Math::XYZTVector>*
	normGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	normGsfEles_charge_ = 0; //std::vector<Int_t>*

	normGsfEles_Et_             = 0; // std::vector<Double_t>*
	normGsfEles_HEEP_Et_        = 0; // std::vector<Double_t>*
	normGsfEles_Eta_            = 0; // std::vector<Double_t>*
	normGsfEles_scEta_          = 0; // std::vector<Double_t>*
	normGsfEles_ecalDriven_     = 0; // std::vector<bool>*
	normGsfEles_ecalDrivenSeed_ = 0; // std::vector<bool>*

	normGsfEles_HEEP_dEtaIn_     = 0; // std::vector<Double_t>*
	normGsfEles_HEEP_dPhiIn_     = 0; // std::vector<Double_t>*
	normGsfEles_HoverE_          = 0; // std::vector<Double_t>*
	normGsfEles_sigmaIetaIeta_   = 0; // std::vector<Double_t>*
	normGsfEles_scSigmaIetaIeta_ = 0; // std::vector<Double_t>*

	normGsfEles_dr03EmIsoEt_        = 0; // std::vector<Double_t>*
	normGsfEles_dr03HadDepth1IsoEt_ = 0; // std::vector<Double_t>*
	normGsfEles_dr03HadDepth2IsoEt_ = 0; // std::vector<Double_t>*
	normGsfEles_dr03TkIsoPt_        = 0; // std::vector<Double_t>*

	normGsfEles_e2x5Max_ = 0; // std::vector<Double_t>*
	normGsfEles_e5x5_    = 0; // std::vector<Double_t>*

	//*** Member variables for heep::Ele variable branches ***
	//kinematic and geometric methods
	normHEEPEles_et_ = 0;
	normHEEPEles_gsfEt_ = 0;
	normHEEPEles_scEt_ = 0;
	normHEEPEles_energy_ = 0;
	normHEEPEles_gsfEnergy_ = 0;
	normHEEPEles_caloEnergy_ = 0;
	normHEEPEles_eta_ = 0;
	normHEEPEles_scEta_ = 0;
	normHEEPEles_detEta_ = 0;
	normHEEPEles_detEtaAbs_ = 0;
	normHEEPEles_phi_ = 0;
	normHEEPEles_scPhi_ = 0;
	normHEEPEles_detPhi_ = 0;
	normHEEPEles_zVtx_ = 0;
	normHEEPEles_p4ptr_ = 0;
	normHEEPEles_gsfP4ptr_ = 0;

	//classification (couldnt they have just named it 'type')
	 normHEEPEles_classification_ = 0;
	normHEEPEles_isEcalDriven_ = 0;
	normHEEPEles_isTrackerDriven_ = 0;
	normHEEPEles_isEB_ = 0;
	normHEEPEles_isEE_ = 0;

	//track methods
	normHEEPEles_charge_ = 0;
	normHEEPEles_trkCharge_ = 0;
	normHEEPEles_pVtx_ = 0;
	normHEEPEles_pCalo_ = 0;
	normHEEPEles_ptVtx_ = 0;
	normHEEPEles_ptCalo_ = 0;

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	normHEEPEles_hOverE_ = 0;
	normHEEPEles_dEtaIn_ = 0;
	normHEEPEles_dPhiIn_ = 0;
	normHEEPEles_dPhiOut_ = 0;
	normHEEPEles_epIn_ = 0;
	normHEEPEles_epOut_ = 0;
	normHEEPEles_fbrem_ = 0;
	normHEEPEles_bremFrac_ = 0;
	normHEEPEles_invEOverInvP_ = 0;

	//shower shape variables
	normHEEPEles_sigmaEtaEta_ = 0;
	normHEEPEles_sigmaEtaEtaUnCorr_ = 0;
	normHEEPEles_sigmaIEtaIEta_ = 0;
	normHEEPEles_e1x5_ = 0;
	normHEEPEles_e2x5Max_ = 0;
	normHEEPEles_e5x5_ = 0;
	normHEEPEles_e1x5Over5x5_ = 0;
	normHEEPEles_e2x5MaxOver5x5_ = 0;

	//isolation, we use cone of 0.3
	normHEEPEles_isolEm_ = 0;
	normHEEPEles_isolHad_ = 0;
	normHEEPEles_isolHadDepth1_ = 0;
	normHEEPEles_isolHadDepth2_ = 0;
	normHEEPEles_isolPtTrks_ = 0;
	normHEEPEles_isolEmHadDepth1_ = 0;

	normEles_.clear();

	//Variables for storing contents of special reconstruction GSF electron branches ...
	bstdGsfEles_number_ = 0;

	bstdGsfEles_OLDp4_  = 0; // std::vector<ROOT::Math::XYZTVector>*
	bstdGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	bstdGsfEles_charge_ = 0; //std::vector<Int_t>*

	bstdGsfEles_Et_             = 0; // std::vector<Double_t>*
	bstdGsfEles_HEEP_Et_        = 0; // std::vector<Double_t>*
	bstdGsfEles_Eta_            = 0; // std::vector<Double_t>*
	bstdGsfEles_scEta_          = 0; // std::vector<Double_t>*
	bstdGsfEles_ecalDriven_     = 0; // std::vector<bool>*
	bstdGsfEles_ecalDrivenSeed_ = 0; // std::vector<bool>*

	bstdGsfEles_HEEP_dEtaIn_     = 0; // std::vector<Double_t>*
	bstdGsfEles_HEEP_dPhiIn_     = 0; // std::vector<Double_t>*
	bstdGsfEles_HoverE_          = 0; // std::vector<Double_t>*
	bstdGsfEles_sigmaIetaIeta_   = 0; // std::vector<Double_t>*
	bstdGsfEles_scSigmaIetaIeta_ = 0; // std::vector<Double_t>*

	bstdGsfEles_dr03EmIsoEt_        = 0; // std::vector<Double_t>*
	bstdGsfEles_dr03HadDepth1IsoEt_ = 0; // std::vector<Double_t>*
	bstdGsfEles_dr03HadDepth2IsoEt_ = 0; // std::vector<Double_t>*
	bstdGsfEles_dr03TkIsoPt_        = 0; // std::vector<Double_t>*

	bstdGsfEles_e2x5Max_ = 0; // std::vector<Double_t>*
	bstdGsfEles_e5x5_    = 0; // std::vector<Double_t>*

	//*** Member variables for heep::Ele variable branches ***
	//kinematic and geometric methods
	bstdHEEPEles_et_ = 0;
	bstdHEEPEles_gsfEt_ = 0;
	bstdHEEPEles_scEt_ = 0;
	bstdHEEPEles_energy_ = 0;
	bstdHEEPEles_gsfEnergy_ = 0;
	bstdHEEPEles_caloEnergy_ = 0;
	bstdHEEPEles_eta_ = 0;
	bstdHEEPEles_scEta_ = 0;
	bstdHEEPEles_detEta_ = 0;
	bstdHEEPEles_detEtaAbs_ = 0;
	bstdHEEPEles_phi_ = 0;
	bstdHEEPEles_scPhi_ = 0;
	bstdHEEPEles_detPhi_ = 0;
	bstdHEEPEles_zVtx_ = 0;
	bstdHEEPEles_p4ptr_ = 0;
	bstdHEEPEles_gsfP4ptr_ = 0;

	//classification (couldnt they have just named it 'type')
	 bstdHEEPEles_classification_ = 0;
	bstdHEEPEles_isEcalDriven_ = 0;
	bstdHEEPEles_isTrackerDriven_ = 0;
	bstdHEEPEles_isEB_ = 0;
	bstdHEEPEles_isEE_ = 0;

	//track methods
	bstdHEEPEles_charge_ = 0;
	bstdHEEPEles_trkCharge_ = 0;
	bstdHEEPEles_pVtx_ = 0;
	bstdHEEPEles_pCalo_ = 0;
	bstdHEEPEles_ptVtx_ = 0;
	bstdHEEPEles_ptCalo_ = 0;

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	bstdHEEPEles_hOverE_ = 0;
	bstdHEEPEles_dEtaIn_ = 0;
	bstdHEEPEles_dPhiIn_ = 0;
	bstdHEEPEles_dPhiOut_ = 0;
	bstdHEEPEles_epIn_ = 0;
	bstdHEEPEles_epOut_ = 0;
	bstdHEEPEles_fbrem_ = 0;
	bstdHEEPEles_bremFrac_ = 0;
	bstdHEEPEles_invEOverInvP_ = 0;

	//shower shape variables
	bstdHEEPEles_sigmaEtaEta_ = 0;
	bstdHEEPEles_sigmaEtaEtaUnCorr_ = 0;
	bstdHEEPEles_sigmaIEtaIEta_ = 0;
	bstdHEEPEles_e1x5_ = 0;
	bstdHEEPEles_e2x5Max_ = 0;
	bstdHEEPEles_e5x5_ = 0;
	bstdHEEPEles_e1x5Over5x5_ = 0;
	bstdHEEPEles_e2x5MaxOver5x5_ = 0;

	//isolation, we use cone of 0.3
	bstdHEEPEles_isolEm_ = 0;
	bstdHEEPEles_isolHad_ = 0;
	bstdHEEPEles_isolHadDepth1_ = 0;
	bstdHEEPEles_isolHadDepth2_ = 0;
	bstdHEEPEles_isolPtTrks_ = 0;
	bstdHEEPEles_isolEmHadDepth1_ = 0;

	bstdEles_.clear();

}

//================================================================//
void BstdZeeFirstAnalyser::InitialiseReconValidationHistos(){
	//normEles_reconValidationHistos_ = ReconValidationHistos("h_normEles_","standard","");
	/*hist_normEles_number_ = new TH1D("h_normEles_number", "Histogram of number of standard GSF eles in each event; Number of standard GSF electrons in event; Number of events", 6, -0.5, 5.5);
	hMin_eleCharge = ;
	hist_normEles_charge_ = new TH1D("h_normEles_charge", "title; x-axis; y-axis", numbins, min, max); //hMin_eleEt
	hist_normEles_Et_ = new TH1D("h_normEles_Et", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_heepEt_ = new TH1D("h_normEles_heepEt", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_Eta_ = new TH1D("h_normEles_Eta", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_scEta_ = new TH1D("h_normEles_scEta", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_ecalDriven_ = new TH1D("h_normEles_ecalDriven", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_ecalDrivenSeed_ = new TH1D("h_normEles_ecalDrivenSeed", "title; x-axis; y-axis", numbins, min, max);

	hist_normEles_heepdEtaIn_ = new TH1D("h_normEles_heepdEtaIn", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_heepdPhiIn_ = new TH1D("h_normEles_heepdPhiIn", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_HoverE_ = new TH1D("h_normEles_HoverE", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_sigmaIetaIeta_ = new TH1D("h_normEles_sigmaIetaIeta", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_scSigmaIetaIeta_ = new TH1D("h_normEles_scSigmaIetaIeta", "title; x-axis; y-axis", numbins, min, max);

	hist_normEles_dr03EmIsoEt_ = new TH1D("h_normEles_dr03EmIsoEt", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_dr03Had1IsoEt_ = new TH1D("h_normEles_dr03Had1IsoPt", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_dr03Had2IsoEt_ = new TH1D("h_normEles_dr03Had2IsoPt", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_dr03TkIsoPt_ = new TH1D("h_normEles_dr03TkIsoPt", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_e2x5Max_ = new TH1D("h_normEles_e2x5Max", "title; x-axis; y-axis", numbins, min, max);
	hist_normEles_e5x5_ = new TH1D("h_normEles_e5x5", "title; x-axis; y-axis", numbins, min, max);*/
	//Currently a placeholder ...
}

//================================================================//
void BstdZeeFirstAnalyser::SetupBranchLinks(const TFile* inFile_ptr){
	inputFile_tree_->SetBranchAddress("event", &event_);
	//Setting up the pointer links for the event information branches ...
	inputFile_tree_->SetBranchAddress("evt_runNum",&evt_runNum_);   // unsigned int
	inputFile_tree_->SetBranchAddress("evt_lumiSec",&evt_lumiSec_); // unsigned int
	inputFile_tree_->SetBranchAddress("evt_evtNum",&evt_evtNum_);   // unsigned int

	//Setting up the pointer links for the trigger branches ...
	inputFile_tree_->SetBranchAddress("trg_PathA_decision", &trg_PathA_decision_);
	inputFile_tree_->SetBranchAddress("trg_PathA_name",     &trg_PathA_name_ptr_);
	inputFile_tree_->SetBranchAddress("trg_PathA_nameOfLastFilter", &trg_PathA_nameOfLastFilter_ptr_);
	inputFile_tree_->SetBranchAddress("trg_PathA_highestTrigObjEt", &trg_PathA_highestTrigObjEt_);


	if(isMC_){
		//Setting up the pointer links for the MC electron branches ...
		inputFile_tree_->SetBranchAddress("mcEles_number",          &mc_numFinalStateEles_      );
		inputFile_tree_->SetBranchAddress("mcEles_HighestEt_charge",&mcEles_HighestEt_charge_   );
		inputFile_tree_->SetBranchAddress("mcEles_HighestEt_PDGid", &mcEles_HighestEt_PDGid_    );
		inputFile_tree_->SetBranchAddress("mcEles_HighestEt_status",&mcEles_HighestEt_status_   );
		inputFile_tree_->SetBranchAddress("mcEles_HighestEt_pt",    &mcEles_HighestEt_pt_       );
		inputFile_tree_->SetBranchAddress("mcEles_HighestEt_p4",    &mcEles_HighestEt_OLDp4_ptr_);

		inputFile_tree_->SetBranchAddress("mcEles_2ndHighestEt_charge", &mcEles_2ndHighestEt_charge_   );
		inputFile_tree_->SetBranchAddress("mcEles_2ndHighestEt_PDGid",  &mcEles_2ndHighestEt_PDGid_    );
		inputFile_tree_->SetBranchAddress("mcEles_2ndHighestEt_status", &mcEles_2ndHighestEt_status_   );
		inputFile_tree_->SetBranchAddress("mcEles_2ndHighestEt_pt",     &mcEles_2ndHighestEt_pt_       );
		inputFile_tree_->SetBranchAddress("mcEles_2ndHighestEt_p4",     &mcEles_2ndHighestEt_OLDp4_ptr_);

		inputFile_tree_->SetBranchAddress("mcZcandidate_pt",      &mcZcandidate_pt_  ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_eta",     &mcZcandidate_eta_ ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_phi",     &mcZcandidate_phi_ ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_mass",    &mcZcandidate_mass_); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_p4",           &mcZcandidate_OLDp4_ptr_   ); //ROOT::Math::XYZTLorentzVector
		inputFile_tree_->SetBranchAddress("mcZcandidate_dEtaEles",     &mcZcandidate_dEtaEles_    ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_dPhiEles",     &mcZcandidate_dPhiEles_    ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_dREles",       &mcZcandidate_dREles_      ); //Double_t
		inputFile_tree_->SetBranchAddress("mcZcandidate_openingAngle", &mcZcandidate_openingAngle_); //Double_t

		inputFile_tree_->SetBranchAddress("mcZboson_pdgId",  &mcZ_pdgId_ );
		inputFile_tree_->SetBranchAddress("mcZboson_status", &mcZ_status_);
		inputFile_tree_->SetBranchAddress("mcZboson_p4",     &mcZ_p4_ptr_); //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		inputFile_tree_->SetBranchAddress("mcZboson_numDaughters",   &mcZ_numDaughters_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughters_dR",   &mcZ_daughters_dR_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughters_dEta", &mcZ_daughters_dEta_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughters_dPhi", &mcZ_daughters_dPhi_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughters_openingAngle", &mcZ_daughters_openingAngle_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughterA_p4", &mcZ_daughterA_OLDp4_ptr_);
		inputFile_tree_->SetBranchAddress("mcZboson_daughterB_p4", &mcZ_daughterB_OLDp4_ptr_);
	}

	//Setting up the pointer links for the 'normal' GSF electron branches ...
	inputFile_tree_->SetBranchAddress("normGsfEles_number", &normGsfEles_number_); //unsigned int
	inputFile_tree_->SetBranchAddress("normGsfEles_p4ptr_", &normGsfEles_OLDp4_  ); // std::vector<ROOT::Math::XYZTVector>*
	inputFile_tree_->SetBranchAddress("normGsfEles_charge", &normGsfEles_charge_); //std::vector<Int_t>*

	inputFile_tree_->SetBranchAddress("normGsfEles_Et",      &normGsfEles_Et_     );           // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_HEEP_Et", &normGsfEles_HEEP_Et_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_Eta",     &normGsfEles_Eta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_scEta",   &normGsfEles_scEta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_ecalDriven",     &normGsfEles_ecalDriven_);         // std::vector<bool>*
	inputFile_tree_->SetBranchAddress("normGsfEles_ecalDrivenSeed", &normGsfEles_ecalDrivenSeed_); // std::vector<bool>*

	inputFile_tree_->SetBranchAddress("normGsfEles_HEEP_dEtaIn", & normGsfEles_HEEP_dEtaIn_);    // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_HEEP_dPhiIn", &normGsfEles_HEEP_dPhiIn_);     // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_HoverE", &normGsfEles_HoverE_);               // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_sigmaIetaIeta", &normGsfEles_sigmaIetaIeta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_scSigmaIetaIeta", &normGsfEles_scSigmaIetaIeta_); // std::vector<Double_t>*

	inputFile_tree_->SetBranchAddress("normGsfEles_dr03EmIsoEt", &normGsfEles_dr03EmIsoEt_);             // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_dr03HadDepth1IsoEt", &normGsfEles_dr03HadDepth1IsoEt_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_dr03HadDepth2IsoEt", &normGsfEles_dr03HadDepth2IsoEt_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_dr03TkIsoPt", &normGsfEles_dr03TkIsoPt_);  // std::vector<Double_t>*

	inputFile_tree_->SetBranchAddress("normGsfEles_e2x5Max", &normGsfEles_e2x5Max_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("normGsfEles_e5x5", &normGsfEles_e5x5_); // std::vector<Double_t>*

	//*** Branches containing heep::Ele method values for norm eles ***
	//kinematic and geometric methods
	inputFile_tree_->SetBranchAddress("normHEEPEles_et", &normHEEPEles_et_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfEt", &normHEEPEles_gsfEt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scEt", &normHEEPEles_scEt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_energy", &normHEEPEles_energy_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfEnergy", &normHEEPEles_gsfEnergy_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_caloEnergy", &normHEEPEles_caloEnergy_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_eta", &normHEEPEles_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scEta", &normHEEPEles_scEta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_detEta", &normHEEPEles_detEta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_detEtaAbs", &normHEEPEles_detEtaAbs_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_phi", &normHEEPEles_phi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scPhi", &normHEEPEles_scPhi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_detPhi", &normHEEPEles_detPhi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_zVtx", &normHEEPEles_zVtx_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_p4ptr", &normHEEPEles_p4ptr_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfP4ptr", &normHEEPEles_gsfP4ptr_);

	//classification (couldnt they have just named it 'type')
	inputFile_tree_->SetBranchAddress("normHEEPEles_classification", &normHEEPEles_classification_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEcalDriven", &normHEEPEles_isEcalDriven_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isTrackerDriven", &normHEEPEles_isTrackerDriven_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEB", &normHEEPEles_isEB_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEE", &normHEEPEles_isEE_);

	//track methods
	inputFile_tree_->SetBranchAddress("normHEEPEles_charge", &normHEEPEles_charge_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_trkCharge", &normHEEPEles_trkCharge_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_pVtx", &normHEEPEles_pVtx_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_pCalo", &normHEEPEles_pCalo_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_ptVtx", &normHEEPEles_ptVtx_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_ptCalo", &normHEEPEles_ptCalo_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	inputFile_tree_->SetBranchAddress("normHEEPEles_hOverE", &normHEEPEles_hOverE_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_dEtaIn", &normHEEPEles_dEtaIn_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_dPhiIn", &normHEEPEles_dPhiIn_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_dPhiOut", &normHEEPEles_dPhiOut_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_epIn", &normHEEPEles_epIn_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_epOut", &normHEEPEles_epOut_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_fbrem", &normHEEPEles_fbrem_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_bremFrac", &normHEEPEles_bremFrac_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_invEOverInvP", &normHEEPEles_invEOverInvP_);

	//shower shape variables
	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaEtaEta", &normHEEPEles_sigmaEtaEta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaEtaEtaUnCorr", &normHEEPEles_sigmaEtaEtaUnCorr_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaIEtaIEta", &normHEEPEles_sigmaIEtaIEta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e1x5", &normHEEPEles_e1x5_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e2x5Max", &normHEEPEles_e2x5Max_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e5x5", &normHEEPEles_e5x5_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e1x5Over5x5", &normHEEPEles_e1x5Over5x5_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e2x5MaxOver5x5", &normHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolEm", &normHEEPEles_isolEm_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHad", &normHEEPEles_isolHad_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHadDepth1", &normHEEPEles_isolHadDepth1_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHadDepth2", &normHEEPEles_isolHadDepth2_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolPtTrks", &normHEEPEles_isolPtTrks_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolEmHadDepth1", &normHEEPEles_isolEmHadDepth1_);

	/////////////////////////////////////////////////////////////////////////////////////////
	//Setting up the pointer links for the special reconstruction GSF electron branches ...
	inputFile_tree_->SetBranchAddress("bstdGsfEles_number", &bstdGsfEles_number_); //unsigned int
	inputFile_tree_->SetBranchAddress("bstdGsfEles_p4ptr_", &bstdGsfEles_OLDp4_  ); // std::vector<ROOT::Math::XYZTVector>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_charge", &bstdGsfEles_charge_); //std::vector<Int_t>*

	inputFile_tree_->SetBranchAddress("bstdGsfEles_Et",      &bstdGsfEles_Et_     );           // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_HEEP_Et", &bstdGsfEles_HEEP_Et_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_Eta",     &bstdGsfEles_Eta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_scEta",   &bstdGsfEles_scEta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_ecalDriven",     &bstdGsfEles_ecalDriven_);         // std::vector<bool>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_ecalDrivenSeed", &bstdGsfEles_ecalDrivenSeed_); // std::vector<bool>*

	inputFile_tree_->SetBranchAddress("bstdGsfEles_HEEP_dEtaIn", & bstdGsfEles_HEEP_dEtaIn_);    // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_HEEP_dPhiIn", &bstdGsfEles_HEEP_dPhiIn_);     // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_HoverE", &bstdGsfEles_HoverE_);               // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_sigmaIetaIeta", &bstdGsfEles_sigmaIetaIeta_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_scSigmaIetaIeta", &bstdGsfEles_scSigmaIetaIeta_); // std::vector<Double_t>*

	inputFile_tree_->SetBranchAddress("bstdGsfEles_dr03EmIsoEt", &bstdGsfEles_dr03EmIsoEt_);             // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_dr03HadDepth1IsoEt", &bstdGsfEles_dr03HadDepth1IsoEt_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_dr03HadDepth2IsoEt", &bstdGsfEles_dr03HadDepth2IsoEt_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_dr03TkIsoPt", &bstdGsfEles_dr03TkIsoPt_);  // std::vector<Double_t>*

	inputFile_tree_->SetBranchAddress("bstdGsfEles_e2x5Max", &bstdGsfEles_e2x5Max_); // std::vector<Double_t>*
	inputFile_tree_->SetBranchAddress("bstdGsfEles_e5x5", &bstdGsfEles_e5x5_); // std::vector<Double_t>*

	//*** Branches containing heep::Ele method values for special reco'n eles ***
	//kinematic and geometric methods
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_et", &bstdHEEPEles_et_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_gsfEt", &bstdHEEPEles_gsfEt_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_scEt", &bstdHEEPEles_scEt_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_energy", &bstdHEEPEles_energy_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_gsfEnergy", &bstdHEEPEles_gsfEnergy_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_caloEnergy", &bstdHEEPEles_caloEnergy_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_eta", &bstdHEEPEles_eta_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_scEta", &bstdHEEPEles_scEta_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_detEta", &bstdHEEPEles_detEta_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_detEtaAbs", &bstdHEEPEles_detEtaAbs_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_phi", &bstdHEEPEles_phi_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_scPhi", &bstdHEEPEles_scPhi_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_detPhi", &bstdHEEPEles_detPhi_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_zVtx", &bstdHEEPEles_zVtx_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_p4ptr", &bstdHEEPEles_p4ptr_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_gsfP4ptr", &bstdHEEPEles_gsfP4ptr_);

	//classification (couldnt they have just named it 'type')
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_classification", &bstdHEEPEles_classification_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isEcalDriven", &bstdHEEPEles_isEcalDriven_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isTrackerDriven", &bstdHEEPEles_isTrackerDriven_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isEB", &bstdHEEPEles_isEB_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isEE", &bstdHEEPEles_isEE_);

	//track methods
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_charge", &bstdHEEPEles_charge_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_trkCharge", &bstdHEEPEles_trkCharge_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_pVtx", &bstdHEEPEles_pVtx_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_pCalo", &bstdHEEPEles_pCalo_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_ptVtx", &bstdHEEPEles_ptVtx_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_ptCalo", &bstdHEEPEles_ptCalo_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_hOverE", &bstdHEEPEles_hOverE_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_dEtaIn", &bstdHEEPEles_dEtaIn_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_dPhiIn", &bstdHEEPEles_dPhiIn_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_dPhiOut", &bstdHEEPEles_dPhiOut_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_epIn", &bstdHEEPEles_epIn_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_epOut", &bstdHEEPEles_epOut_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_fbrem", &bstdHEEPEles_fbrem_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_bremFrac", &bstdHEEPEles_bremFrac_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_invEOverInvP", &bstdHEEPEles_invEOverInvP_);

	//shower shape variables
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_sigmaEtaEta", &bstdHEEPEles_sigmaEtaEta_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_sigmaEtaEtaUnCorr", &bstdHEEPEles_sigmaEtaEtaUnCorr_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_sigmaIEtaIEta", &bstdHEEPEles_sigmaIEtaIEta_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_e1x5", &bstdHEEPEles_e1x5_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_e2x5Max", &bstdHEEPEles_e2x5Max_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_e5x5", &bstdHEEPEles_e5x5_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_e1x5Over5x5", &bstdHEEPEles_e1x5Over5x5_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_e2x5MaxOver5x5", &bstdHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolEm", &bstdHEEPEles_isolEm_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolHad", &bstdHEEPEles_isolHad_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolHadDepth1", &bstdHEEPEles_isolHadDepth1_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolHadDepth2", &bstdHEEPEles_isolHadDepth2_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolPtTrks", &bstdHEEPEles_isolPtTrks_);
	inputFile_tree_->SetBranchAddress("bstdHEEPEles_isolEmHadDepth1", &bstdHEEPEles_isolEmHadDepth1_);

	std::cout << "   pointer links for branch read-out now set up." << std::endl;
}

//==============================================================
void BstdZeeFirstAnalyser::SetupEMuMethodTrees(){
	eleMuTree_ = new TTree("myTree","eleMu data");
	eleMuTree_->Branch("mass",   &eleMuTreeVar_mass_,   "mass/D");
	eleMuTree_->Branch("weight", &eleMuTreeVar_weight_, "weight/D");

	diEleTree_ = new TTree("myTree","diEle data");
	diEleTree_->Branch("mass",   &diEleTreeVar_mass_,   "mass/D");
	diEleTree_->Branch("weight", &diEleTreeVar_weight_, "weight/D");
}


//-----------------------------------------------------------------------------------------//
//----- Methods called by the destructor (i.e. deletion of heap-based objects)...     -----//
//-----------------------------------------------------------------------------------------//
void BstdZeeFirstAnalyser::DeleteReconValidationHistos(){
	/*delete hist_normEles_number_;
	delete hist_normEles_charge_;
	delete hist_normEles_Et_;
	delete hist_normEles_heepEt_;
	delete hist_normEles_Eta_;
	delete hist_normEles_scEta_;
	delete hist_normEles_ecalDriven_;
	delete hist_normEles_ecalDrivenSeed_;

	delete hist_normEles_heepdEtaIn_;
	delete hist_normEles_heepdPhiIn_;
	delete hist_normEles_HoverE_;
	delete hist_normEles_sigmaIetaIeta_;
	delete hist_normEles_scSigmaIetaIeta_;

	delete hist_normEles_dr03EmIsoEt_;
	delete hist_normEles_dr03Had1IsoEt_;
	delete hist_normEles_dr03Had2IsoEt_;
	delete hist_normEles_dr03TkIsoPt_;
	delete hist_normEles_e2x5Max_;
	delete hist_normEles_e5x5_;*/
}

//=====================================================================================================
//-----------------------------------------------------------------------------------------------------

int main()
{
	bool testClass = true;
	
	if(testClass){
		gROOT->ProcessLine("#include <vector>"); //Without this line, the std::vector<bool> branches are not linked correctly(???), giving 'dictionary for class...' errors
				//upon running of ZEventAnalyser program, and accessing these branches results in nonsensical results (e.g. methods such as size()).

		/*TString outFileTag = "2011-06-23";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/sampleDataFile_NTuple-100evts_2011-Apr-24.root";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/sampleSignalmcFile_NTuple-100evts_2011-Apr-25.root";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/Data/samplemcFile_Fall10_NTuple_2011-Apr-25.root";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/testNTuple_2-00TeVu.root";
		TString myOutputFileName = "histos_DYJetsToEE_" + outFileTag + ".root";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/dataNTuple_170kEvts_2011-05-13.root";
		//TString myInputFileName  = "/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/bkgdMC_41X-Ntuple_DYJetsToLL_eeOnly_66795Evts_2011-06-22.root";
		//TString myOutputFileName = "histos_DYJetsToLL_" + outFileTag + ".root";

		std::cout << std::endl << std::endl;
		std::cout << "  ***-------------------------------***" << std::endl;
		std::cout << "  *** Data analysis ...." << std::endl;
		//BstdZeeFirstAnalyser myAnalyser(0, 2530516, false, myInputFileName, myOutputFileName, -1, 180, 120, 1200.0);
		BstdZeeFirstAnalyser* myAnalyser = new BstdZeeFirstAnalyser(0, 20, true, "/home/ppd/nnd85574/Work/BstdZee/CMSSW_4_2_4_patch1/src/NTupler/BstdZeeNTupler/testNTuple.root", "testOutputHistos.root", 6, 180, 120, 1200.0);
		std::cout << " Running the DoAnalysis method ..." << std::endl;
		myAnalyser->DoAnalysis( 1.0 ); //myAnalyser.DoAnalysis( 1.0*(2321000.0/2530516.0) );
		delete myAnalyser;

		TString sigInputFilename;  //sigInputFilenames.clear();
		TString sigOutputFilename; //sigOutputFilenames.clear();
		//std::vector<BstdZeeFirstAnalyser> sigAnalysers; sigAnalysers.clear();

		//sigInputFilename =  "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/signalmcNTuple_0-75TeVu_20kEvts_2011-05-13.root" ;
		sigInputFilename =  "/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/sigMC_41X-Ntuple_0-75TeVu_20kEvts_2011-06-23.root";
		sigOutputFilename = "histos_0-75TeVu_" + outFileTag + ".root" ;
		std::cout << std::endl << std::endl;
		std::cout << "  ***-------------------------------***" << std::endl;
		std::cout << "  *** Signal analysis A ..." << std::endl;
		//BstdZeeFirstAnalyser sigAnalyserA(0, 20000, true, sigInputFilename, sigOutputFilename, -1, 180, 30, 900.0);
		std::cout << " Running the DoAnalysis method ..." << std::endl;
		//sigAnalyserA.DoAnalysis( 1.0 );  //sigAnalyserA.DoAnalysis( 1.0*(925.0/20000.0) );

		//sigInputFilename = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/signalmcNTuple_1-00TeVu_20kEvts_2011-05-13.root" ;
		sigInputFilename =  "/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/sigMC_41X-Ntuple_1-00TeVu_20kEvts_2011-06-23.root";
		sigOutputFilename = "histos_1-00TeVu_" + outFileTag + ".root" ;
		std::cout << std::endl << std::endl;
		std::cout << "  ***-------------------------------***" << std::endl;
		std::cout << "  *** Signal analysis B ..." << std::endl;
		//BstdZeeFirstAnalyser sigAnalyserB(0, 20000, true, sigInputFilename, sigOutputFilename, -1, 180, 30, 900.0);
		std::cout << " Running the DoAnalysis method ..." << std::endl;
		//sigAnalyserB.DoAnalysis( 1.0 ); //sigAnalyserB.DoAnalysis( 1.0*(110.0/20000.0) );

		//sigInputFilename = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/signalmcNTuple_2-00TeVu_20kEvts_2011-05-13.root" ;
		sigInputFilename =  "/opt/ppd/scratch/williams/Datafiles/NTuples/CMSSW_41X/sigMC_41X-Ntuple_2-00TeVu_20kEvts_2011-06-23.root";
		sigOutputFilename = "histos_2-00TeVu_" + outFileTag + ".root" ;
		std::cout << std::endl << std::endl;
		std::cout << "  ***-------------------------------***" << std::endl;
		std::cout << "  *** Signal analysis C ..." << std::endl;
		//BstdZeeFirstAnalyser sigAnalyserC(0, 20000, true, sigInputFilename, sigOutputFilename, -1, 180, 30, 1800.0);
		std::cout << " Running the DoAnalysis method ..." << std::endl;
		//sigAnalyserC.DoAnalysis( 1.0 ); //sigAnalyserC.DoAnalysis( 1.0*(1.2/20000.0) );*/

		TString inFilePrefix = "/home/ppd/nnd85574/Work/BstdZee/CMSSW_4_2_4_patch1/src/NTupler/BstdZeeNTupler/";
		std::vector<BstdZeeFirstAnalyser> myAnalysers; myAnalysers.clear();
		std::vector<TString> inFile; inFile.clear();
		std::vector<TString> outFile; outFile.clear();
		std::vector<unsigned int> nEvents; nEvents.clear();
		std::vector<Double_t> intLumiPerEvent; intLumiPerEvent.clear();
		Double_t desiredIntLumi = 0.001;

		// ttbar
		inFile.push_back("ttbarNTuple_2011-07-14.root");
		outFile.push_back("ttbarOutput_2011-07-14");
		nEvents.push_back( 5 );
		intLumiPerEvent.push_back(11.6/1090000.0); // in inv fb

		// DYToTauTau
		inFile.push_back("DYToTauTauNTuple_2011-07-14.root");
		outFile.push_back("DYToTauTauOutput_2011-07-14");
		nEvents.push_back( 5 );
		intLumiPerEvent.push_back( 1.6/2032536.0 ); // in inv fb

		// data fake A - ttbar
		inFile.push_back("dataFakeANTuple_2011-07-14.root");
		outFile.push_back("dataFakeAOutput_2011-07-14");
		nEvents.push_back(5);
		intLumiPerEvent.push_back( 11.6/1090000.0 ); // in inv fb

		// data fake B - DYToTauTau
		inFile.push_back("dataFakeBNTuple_2011-07-14.root");
		outFile.push_back("dataFakeBOutput_2011-07-14");
		nEvents.push_back( 9000 );
		intLumiPerEvent.push_back( 1.6/2032536.0 ); // in inv fb


		// Doing the analysis ...
		for(unsigned int idx=0; idx<inFile.size(); idx++){
			//TStopwatch watch;
			//watch.Start();
			std::cout << std::endl << "***--- NEW FILE (no " << idx+1 << "/" << inFile.size() << ") ---***" << std::endl;
			std::cout << " * Input:  " << inFile.at(idx) << std::endl;
			std::cout << " * Output: " << outFile.at(idx) << std::endl;
			BstdZeeFirstAnalyser myAnalysers = BstdZeeFirstAnalyser(0, nEvents.at(idx), true, inFilePrefix+inFile.at(idx), outFile.at(idx), -1, 180, 120, 1200.0);
			std::cout << " Running the DoAnalysis method ..." << std::endl;
			myAnalysers.DoAnalysis(  desiredIntLumi/( static_cast<double>(nEvents.at(idx))*intLumiPerEvent.at(idx))  ); //myAnalyser.DoAnalysis( 1.0*(2321000.0/2530516.0) );

			//watch.Stop();
			//std::cout << " * " << nEvents.at(idx) << " events analysed in " << watch.Print() << " seconds." << std::endl;
			//delete myAnalysers.at(idx);
		}

		/*for(unsigned int i=0; i<sigInputFilenames.size()-1; i++){
			std::cout << std::endl << std::endl;
			std::cout << "  ***-------------------------------***" << std::endl;
			std::cout << "  *** Signal analysis #" << i << " ..." << std::endl;
			sigAnalysers.push_back( BstdZeeFirstAnalyser(0, 20000, false, sigInputFilenames.at(i), sigOutputFilenames.at(i), -1) );
			std::cout << " Running the DoAnalysis method ..." << std::endl;
			sigAnalysers.at(i).DoAnalysis(1.0);
		}*/
	}
	else
	{
	/*	std::cout << "BstdZeeFirstAnalyser program is starting..." << std::endl;
		gROOT->ProcessLine("#include <vector>"); //Without this line, the std::vector<bool> branches are not linked correctly(???), giving 'dictionary for class...' errors
		//upon running of ZEventAnalyser program, and accessing these branches results in nonsensical results (e.g. methods such as size()).

		//Variable declarations
		int vFlg = 1;
		int num_evts = 20;
		Double_t evtWeight = 1.0;
		//Double_t crossSection = -999.9; //in *fb*
		int sampleType = 3; //1 for 0.75TeV u*, 2 for 1.00TeVu*, and 3 for 2.00TeVu*

		TString inputFileName("tmpFileName");
		TString outputFileName("tmpFileName");
	 */
		/*if(sampleType==1){
			inputFileName = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_3_11_0/src/MCNTupler/BstdZeeMCTruthNTupler/MCTruthEventData_0-75TeVu_2011-03-16_1000evts.root";
			outputFileName = "BstdZeeMCTruthHistos_0-75TeVu_2011-03-16.root";
			//crossSection = 2.5*1000.0*1000.0; //in fb
		}
		else if(sampleType==2){
			inputFileName = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_3_11_0/src/MCNTupler/BstdZeeMCTruthNTupler/MCTruthEventData_1-00TeVu_2011-03-16_1000evts.root";
			outputFileName = "BstdZeeMCTruthHistos_1-00TeVu_2011-03-16.root";
			//crossSection = *1000.0*1000.0; //in fb
		}
		else if(sampleType==3){
			inputFileName = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_3_11_0/src/MCNTupler/BstdZeeMCTruthNTupler/MCTruthEventData_2-00TeVu_2011-03-16_1000evts.root";
			outputFileName = "BstdZeeMCTruthHistos_2-00TeVu_2011-03-16.root";
			//crossSection = 2.5*1000.0*1000.0; //in fb
		}
		else
			return 0;*/
	/*
		inputFileName  = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/NTupler/BstdZeeNTupler/eventData_ntuple.root";
		outputFileName = "/opt/ppd/scratch/williams/EDAnalysers/BstdZee/CMSSW_4_1_4_patch2/src/BstdZeeFirst/Analyser/testOutputHists.root";
		//evtWeight = (crossSection*1.0)/static_cast<double>(num_evts);
		evtWeight = 1.0/static_cast<double>(num_evts);
	
		//------------------------------------------------
		//Variables that will be filled with branch data...
		unsigned int mc_numFinalStateEles = 0;
	
		Int_t mcEles_HighestEt_charge = -999;
		Int_t mcEles_HighestEt_PDGid  = -999;
		Int_t mcEles_HighestEt_status = -999;
		Double_t mcEles_HighestEt_pt  = -999.9;
		//ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcEles_HighestEt_OLDp4_ptr = 0; //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcEles_HighestEt_p4(0.0,0.0,0.0,0.0);
	
		Int_t mcEles_2ndHighestEt_charge = -999;
		Int_t mcEles_2ndHighestEt_PDGid  = -999;
		Int_t mcEles_2ndHighestEt_status = -999;
		Double_t mcEles_2ndHighestEt_pt  = -999.9;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcEles_2ndHighestEt_OLDp4_ptr = 0; //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcEles_2ndHighestEt_p4(0.0,0.0,0.0,0.0);
	
		Double_t mcZcandidate_pt  = -999.9;
		Double_t mcZcandidate_eta = -999.9;
		Double_t mcZcandidate_phi = -999.9;
		Double_t mcZcandidate_mass= -999.9;
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* mcZcandidate_OLDp4_ptr = 0; //ROOT::Math::XYZTVector
																															//NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		TLorentzVector mcZcandidate_p4(0.0,0.0,0.0,0.0);
		Double_t mcZcandidate_dEtaEles = -999.9;
		Double_t mcZcandidate_dPhiEles = -999.9;
		Double_t mcZcandidate_dREles   = -999.9;
		Double_t mcZcandidate_openingAngle = -999.9;
	
		Bool_t hltPathA_decision = false; // Bool_t
		std::string *hltPathA_name_ptr = 0;     // std::string hltPathA_
	
		//std::vector<Double_t> *ele_SCs_Et = 0;
		//std::vector<bool> *ele_ecalDrivenFlags = 0;
		//ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > diele_p4(0.0,0.0,0.0,0.0); // N.B. ROOT::Math::XYZTVector is typedef for ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
	
		Int_t numBins_PtZ = 24;
		Double_t min_PtZ = 0.0;
		Double_t max_PtZ = 1200.0;
		
		Double_t ptBinLims_Z[21];
		Int_t numBins_LogPtZ = 20;
		//std::cout << std::endl <<"ptBinLims_Z:" << std::endl;
		for(int idx=0; idx<=numBins_LogPtZ; idx++){
			ptBinLims_Z[idx] = 0.1*pow(10.0,static_cast<double>(idx)/5.0);
			//std::cout << "    " << ptBinLims_Z[idx] << std::endl;
		}
		
		Int_t numBins_PtEle = 25;
		Double_t min_PtEle  = 0.0;
		Double_t max_PtEle  = 1000.0;

		Double_t ptBinLims_ele[21];
		Int_t numBins_LogPtEle = 20;
		//std::cout << std::endl <<"ptBinLims_Z:" << std::endl;
		for(int idx=0; idx<=numBins_LogPtEle; idx++){
			ptBinLims_ele[idx] = 5.0*pow(10.0,static_cast<double>(idx)/10.0);
			//std::cout << "    " << ptBinLims_Z[idx] << std::endl;
		}
		double value_pi = 3.14159265;

		//Histogram creation...
		TH1D* mcHighestEtEle_PtHist            = new TH1D("mcHighestEtEle_PtHist", "p_{T} distribution of highest p_{T} electron; electron p_{T} /GeVc^{-1}; Fraction per 40GeV", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mcHighestEtEle_LogPtHist         = new TH1D("mcHighestEtEle_LogPtHist", "p_{T} distribution of highest p_{T} electron; electron p_{T} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T}/GeVc^{-1})", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mcHighestEtEle_PostTrigPtHist    = new TH1D("mcHighestEtEle_PostTrigPtHist", "p_{T} distribution of highest p_{T} electron in triggered events; electron p_{T} /GeVc^{-1}; Fraction per 40GeV", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mcHighestEtEle_PostTrigLogPtHist = new TH1D("mcHighestEtEle_PostTrigLogPtHist", "p_{T} distribution of highest p_{T} electron in triggered events; electron p_{T} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T}/GeVc^{-1})", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mcHighestEtEle_EtaHist         = new TH1D("mcHighestEtEle_EtaHist"        , "#eta distribution of the highest p_{T} electron; Electron pseudorapidity, #eta; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mcHighestEtEle_PostTrigEtaHist = new TH1D("mcHighestEtEle_PostTrigEtaHist", "#eta distribution of the highest p_{T} electron in triggered events; Electron pseudorapidity, #eta; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mcHighestEtEle_PhiHist         = new TH1D("mcHighestEtEle_PhiHist"        , "#phi distribution of the highest p_{T} electron; #phi /radians; Fraction per 0.3", 20, -value_pi, +value_pi);
		TH1D* mcHighestEtEle_PostTrigPhiHist = new TH1D("mcHighestEtEle_PostTrigPhiHist", "#phi distribution of the highest p_{T} electron in triggered events; #phi /radians; Fraction per 0.3", 20, -value_pi, +value_pi);
		//TrigEffi hists...
		TH1D* mcHighestEtEle_PtTrigEffi            = new TH1D("mcHighestEtEle_PtTrigEffi"   , "Trigger efficiency in p_{T} distribution of the highest p_{T} electron; electron p_{T} /GeVc^{-1}; Trigger efficiency", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mcHighestEtEle_LogPtTrigEffi         = new TH1D("mcHighestEtEle_LogPtTrigEffi", "Trigger efficiency in p_{T} of the highest p_{T} electron; electron p_{T} /GeVc^{-1}; Trigger efficiency", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mcHighestEtEle_EtaTrigEffi           = new TH1D("mcHighestEtEle_EtaTrigEffi"  , "Trigger efficiency in #eta of the highest p_{T} electron; Electron pseudorapidity, #eta; Trigger efficiency", 20, -3.0, +3.0);
		TH1D* mcHighestEtEle_PhiTrigEffi           = new TH1D("mcHighestEtEle_PhiTrigEffi"  , "Trigger efficiency in #phi of the highest p_{T} electron; #phi /radians; Trigger efficiency", 20, -value_pi, +value_pi);

		//
		// 2nd highest Et ele hists...
		TH1D* mc2ndHighestEtEle_PtHist            = new TH1D("mc2ndHighestEtEle_PtHist",         "p_{T} distribution of 2nd highest p_{T} electron; electron p_{T} /GeVc^{-1}; Fraction per 40GeV", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mc2ndHighestEtEle_LogPtHist         = new TH1D("mc2ndHighestEtEle_LogPtHist",      "p_{T} distribution of 2nd highest p_{T} electron; electron p_{T} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T}/GeVc^{-1})", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mc2ndHighestEtEle_PostTrigPtHist    = new TH1D("mc2ndHighestEtEle_PostTrigPtHist", "p_{T} distribution of 2nd highest p_{T} electron in triggered events; electron p_{T} /GeVc^{-1}; Fraction per 40GeV", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mc2ndHighestEtEle_PostTrigLogPtHist = new TH1D("mc2ndHighestEtEle_PostTrigLogPtHist", "p_{T} distribution of 2nd highest p_{T} electron in triggered events; electron p_{T} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T}/GeVc^{-1})", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mc2ndHighestEtEle_EtaHist         = new TH1D("mc2ndHighestEtEle_EtaHist"        , "#eta distribution of the 2nd highest p_{T} electron; Electron pseudorapidity, #eta; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mc2ndHighestEtEle_PostTrigEtaHist = new TH1D("mc2ndHighestEtEle_PostTrigEtaHist", "#eta distribution of the 2nd highest p_{T} electron in triggered events; Electron pseudorapidity, #eta; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mc2ndHighestEtEle_PhiHist         = new TH1D("mc2ndHighestEtEle_PhiHist"        , "#phi distribution of the 2nd highest p_{T} electron; #phi /radians; Fraction per 0.3", 20, -value_pi, +value_pi);
		TH1D* mc2ndHighestEtEle_PostTrigPhiHist = new TH1D("mc2ndHighestEtEle_PostTrigPhiHist", "#phi distribution of the 2nd highest p_{T} electron in triggered events; #phi /radians; Fraction per 0.3", 20, -value_pi, +value_pi);
		//TrigEffi hists...
		TH1D* mc2ndHighestEtEle_PtTrigEffi            = new TH1D("mc2ndHighestEtEle_PtTrigEffi"   , "Trigger efficiency in p_{T} distribution of the 2nd highest p_{T} electron; electron p_{T} /GeVc^{-1}; Trigger efficiency", numBins_PtEle, min_PtEle, max_PtEle);
		TH1D* mc2ndHighestEtEle_LogPtTrigEffi         = new TH1D("mc2ndHighestEtEle_LogPtTrigEffi", "Trigger efficiency in p_{T} of the 2nd highest p_{T} electron; electron p_{T} /GeVc^{-1}; Trigger efficiency", numBins_LogPtEle, ptBinLims_ele);
		TH1D* mc2ndHighestEtEle_EtaTrigEffi           = new TH1D("mc2ndHighestEtEle_EtaTrigEffi"  , "Trigger efficiency in #eta of the 2nd highest p_{T} electron; Electron pseudorapidity, #eta; Trigger efficiency", 20, -3.0, +3.0);
		TH1D* mc2ndHighestEtEle_PhiTrigEffi           = new TH1D("mc2ndHighestEtEle_PhiTrigEffi"  , "Trigger efficiency in #phi of the 2nd highest p_{T} electron; #phi /radians; Trigger efficiency", 20, -value_pi, +value_pi);

		//
		// Z candidate hists...
		TH1D* mcZcandidate_PtHist      = new TH1D("mcZcandidate_PtHist",   "p_{T} distribution of the Z candidate; transverse momentum, p_{T, Z} /GeVc^{-1}; Fraction per 50GeVc^{-1}", numBins_PtZ, min_PtZ, max_PtZ);
		TH1D* mcZcandidate_LogPtHist   = new TH1D("mcZcandidate_LogPtHist",   "p_{T} distribution of the Z candidate; transverse momentum, p_{T, Z} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T})", numBins_LogPtZ, ptBinLims_Z);
		TH1D* mcZcandidate_EtaHist     = new TH1D("mcZcandidate_EtaHist",  "#eta distribution of the Z candidate; #eta_{Z} /radians; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mcZcandidate_PhiHist     = new TH1D("mcZcandidate_PhiHist",  "#phi distribution of the Z candidate; #phi_{Z} /radians; Fraction per #pi/20", 20, -value_pi, +value_pi);
		TH1D* mcZcandidate_MassHist    = new TH1D("mcZcandidate_MassHist", "Invariant mass distribution of the Z candidate; Invariant mass of dielectron system, M_{ee} /GeVc^{-2}; Fraction per 5GeVc^{-2}", 20, 40.0, 140.0);
		TH1D* mcZcandidate_dEtaHist    = new TH1D("mcZcandidate_dEtaHist", "#Delta#eta_{ee} distribution of the Z candidate; #Delta#eta_{ee}; Fraction per 0.2", 20, 0.0, 4.0);
		TH1D* mcZcandidate_dPhiHist    = new TH1D("mcZcandidate_dPhiHist", "#Delta#phi_{ee} distribution of the Z candidate; #Delta#phi_{ee} /radians; Fraction per 0.1", 20, 0.0, 2.0);
		TH1D* mcZcandidate_dRHist      = new TH1D("mcZcandidate_dRHist",   "#DeltaR_{ee} distribution of the Z candidate; #DeltaR_{ee}; Fraction per 0.1", 20, 0.0, 2.0);
		TH1D* mcZcandidate_openingAngleHist = new TH1D("mcZcandidate_openingAngleHist", "Opening angle distribution of the Z candidate; dielectron opening angle, #theta_{ee} /radians; Fraction per 0.1", 20, 0.0, 2.0);

		TH1D* mcZcandidate_PostTrigPtHist      = new TH1D("mcZcandidate_PostTrigPtHist",   "p_{T} distribution of the Z candidate in triggered events; transverse momentum, p_{T, Z} /GeVc^{-1}; Fraction per 50GeVc^{-1}", numBins_PtZ, min_PtZ, max_PtZ);
		TH1D* mcZcandidate_PostTrigLogPtHist   = new TH1D("mcZcandidate_PostTrigLogPtHist",   "p_{T} distribution of the Z candidate in triggered events; transverse momentum, p_{T, Z} /GeVc^{-1}; Fraction per 0.2 in log_{10}(p_{T})", numBins_LogPtZ, ptBinLims_Z);
		TH1D* mcZcandidate_PostTrigEtaHist     = new TH1D("mcZcandidate_PostTrigEtaHist",  "#eta distribution of the Z candidate in triggered events; #eta_{Z} /radians; Fraction per 0.3", 20, -3.0, +3.0);
		TH1D* mcZcandidate_PostTrigPhiHist     = new TH1D("mcZcandidate_PostTrigPhiHist",  "#phi distribution of the Z candidate in triggered events; #phi_{Z} /radians; Fraction per #pi/20", 20, -value_pi, +value_pi);
		TH1D* mcZcandidate_PostTrigMassHist    = new TH1D("mcZcandidate_PostTrigMassHist", "Invariant mass distribution of the Z candidate in triggered events; Invariant mass of dielectron system, M_{ee} /GeVc^{-2}; Fraction per 5GeVc^{-2}", 20, 40.0, 140.0);
		TH1D* mcZcandidate_PostTrigdEtaHist    = new TH1D("mcZcandidate_PostTrigdEtaHist", "#Delta#eta_{ee} distribution of the Z candidate in triggered events; #Delta#eta_{ee}; Fraction per 0.2", 20, 0.0, 4.0);
		TH1D* mcZcandidate_PostTrigdPhiHist    = new TH1D("mcZcandidate_PostTrigdPhiHist", "#Delta#phi_{ee} distribution of the Z candidate in triggered events; #Delta#phi_{ee} /radians; Fraction per 0.1", 20, 0.0, 2.0);
		TH1D* mcZcandidate_PostTrigdRHist      = new TH1D("mcZcandidate_PostTrigdRHist",   "#DeltaR_{ee} distribution of the Z candidate in triggered events; #DeltaR_{ee}; Fraction per 0.1", 20, 0.0, 2.0);
		TH1D* mcZcandidate_PostTrigopeningAngleHist = new TH1D("mcZcandidate_PostTrigopeningAngleHist", "Opening angle distribution of the Z candidate in triggered events; dielectron opening angle, #theta_{ee} /radians; Fraction per 0.1", 20, 0.0, 2.0);

		TH1D* mcZcandidate_PtTrigEffi      = new TH1D("mcZcandidate_PtTrigEffi",   "Trigger efficiency in p_{T} of the Z candidate; transverse momentum, p_{T, Z} /GeVc^{-1}; Trigger efficiency", numBins_PtZ, min_PtZ, max_PtZ);
		TH1D* mcZcandidate_LogPtTrigEffi   = new TH1D("mcZcandidate_LogPtTrigEffi","Trigger efficiency in p_{T} of the Z candidate; transverse momentum, p_{T, Z} /GeVc^{-1}; Trigger efficiency", numBins_LogPtZ, ptBinLims_Z);
		TH1D* mcZcandidate_EtaTrigEffi     = new TH1D("mcZcandidate_EtaTrigEffi",  "Trigger efficiency in #eta of the Z candidate; #eta_{Z} /radians; Trigger efficiency", 20, -3.0, +3.0);
		TH1D* mcZcandidate_PhiTrigEffi     = new TH1D("mcZcandidate_PhiTrigEffi",  "Trigger efficiency in #phi of the Z candidate; #phi_{Z} /radians; Trigger efficiency", 20, -value_pi, +value_pi);
		TH1D* mcZcandidate_MassTrigEffi    = new TH1D("mcZcandidate_MassTrigEffi", "Trigger efficiency in the Z candidate invariant mass; Invariant mass of dielectron system, M_{ee} /GeVc^{-2}; Trigger efficiency", 20, 40.0, 140.0);
		TH1D* mcZcandidate_dEtaTrigEffi    = new TH1D("mcZcandidate_dEtaTrigEffi", "Trigger efficiency in #Delta#eta_{ee} of the Z candidate; #Delta#eta_{ee}; Trigger efficiency", 20, 0.0, 4.0);
		TH1D* mcZcandidate_dPhiTrigEffi    = new TH1D("mcZcandidate_dPhiTrigEffi", "Trigger efficiency in #Delta#phi_{ee} of the Z candidate; #Delta#phi_{ee} /radians; Trigger efficiency", 20, 0.0, 2.0);
		TH1D* mcZcandidate_dRTrigEffi      = new TH1D("mcZcandidate_dRTrigEffi",   "Trigger efficiency in #DeltaR_{ee} of the Z candidate; #DeltaR_{ee}; Trigger efficiency", 20, 0.0, 2.0);
		TH1D* mcZcandidate_openingAngleTrigEffi = new TH1D("mcZcandidate_openingAngleTrigEffi", "Trigger efficiency in the opening angle of the Z candidate; dielectron opening angle, #theta_{ee} /radians; Trigger efficiency", 20, 0.0, 2.0);

		//-----------------------
		//Opening the datafile(s)
		std::cout << "Opening the ntuple datafile..." << std::endl;
		TFile *f = TFile::Open(inputFileName,"READ");
		if(!f)
			return 0;
		else
			std::cout << "   file opened successfully." << std::endl;

		//Setting up reading of the tree data from file...
		std::cout << "Beginning to analyse the event data..." << std::endl;
		f->cd("demo");
		std::cout << "   changed 'folder' within the ntuple successfully" << std::endl;
		TTree *dataTree = (TTree*)gDirectory->Get("EventDataTree");
		std::cout << "   got the 'EventDataTree' object successfully" << std::endl;

		dataTree->SetBranchAddress("mcEles_number",          &mc_numFinalStateEles      );
		dataTree->SetBranchAddress("mcEles_HighestEt_charge",&mcEles_HighestEt_charge   );
		dataTree->SetBranchAddress("mcEles_HighestEt_PDGid", &mcEles_HighestEt_PDGid    );
		dataTree->SetBranchAddress("mcEles_HighestEt_status",&mcEles_HighestEt_status   );
		dataTree->SetBranchAddress("mcEles_HighestEt_pt",    &mcEles_HighestEt_pt       );
		dataTree->SetBranchAddress("mcEles_HighestEt_p4",    &mcEles_HighestEt_OLDp4_ptr);
		std::cout << "   dealt with the highest Et ele branches successfully" << std::endl;

		dataTree->SetBranchAddress("mcEles_2ndHighestEt_charge", &mcEles_2ndHighestEt_charge   );
		dataTree->SetBranchAddress("mcEles_2ndHighestEt_PDGid",  &mcEles_2ndHighestEt_PDGid    );
		dataTree->SetBranchAddress("mcEles_2ndHighestEt_status", &mcEles_2ndHighestEt_status   );
		dataTree->SetBranchAddress("mcEles_2ndHighestEt_pt",     &mcEles_2ndHighestEt_pt       );
		dataTree->SetBranchAddress("mcEles_2ndHighestEt_p4",     &mcEles_2ndHighestEt_OLDp4_ptr);
		//std::cout << "   dealt with the 2nd highest Et ele branches successfully" << std::endl;

		dataTree->SetBranchAddress("mcZcandidate_pt",      &mcZcandidate_pt  ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_eta",     &mcZcandidate_eta ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_phi",     &mcZcandidate_phi ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_mass",    &mcZcandidate_mass); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_p4",           &mcZcandidate_OLDp4_ptr   ); //ROOT::Math::XYZTLorentzVector
		dataTree->SetBranchAddress("mcZcandidate_dEtaEles",     &mcZcandidate_dEtaEles    ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_dPhiEles",     &mcZcandidate_dPhiEles    ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_dREles",       &mcZcandidate_dREles      ); //Double_t
		dataTree->SetBranchAddress("mcZcandidate_openingAngle", &mcZcandidate_openingAngle); //Double_t
		std::cout << "   dealt with the Z candidate branches successfully" << std::endl;

		dataTree->SetBranchAddress("trg_PathA_decision", &hltPathA_decision); // Bool_t
		dataTree->SetBranchAddress("trg_PathA_name",     &hltPathA_name_ptr);     // std::string hltPathA_
		std::cout << "   dealt with the trigger branches successfully" << std::endl;

		//---------------------------------------
		//Running over the data event-by-event...
		for(Int_t i = 0; i < num_evts; i++){
			//Loading this event's information from the TTree
			Long64_t dataTreeEntry = dataTree->LoadTree(i);
			dataTree->GetEntry(dataTreeEntry);

			if(vFlg>1){std::cout << " Analysing event no. " << i << " (run .., lumi .., evt ..): " << std::endl;}

			//Setting up the TLorentzVector 4momenta...
			//mcEles_HighestEt_p4.SetPxPyPzE(mcEles_HighestEt_OLDp4_ptr->Px(),mcEles_HighestEt_OLDp4_ptr->Py(),mcEles_HighestEt_OLDp4_ptr->Pz(),mcEles_HighestEt_OLDp4_ptr->E());
			mcEles_HighestEt_p4 = ConvertToTLorentzVector(mcEles_HighestEt_OLDp4_ptr);
			mcEles_2ndHighestEt_p4 = ConvertToTLorentzVector(mcEles_2ndHighestEt_OLDp4_ptr);
			mcZcandidate_p4 = ConvertToTLorentzVector(mcZcandidate_OLDp4_ptr);
	
			//Printing the properties of the 2 highest Et eles and the Z candidate, along with the trigger info, to screen...
			if(vFlg>1){
				if(hltPathA_decision)
					std::cout << "  ->This event passed HLT trigger path A (i.e. " << *hltPathA_name_ptr << ")" << std::endl;
				else
					std::cout << "  ->This event did *NOT* pass HLT trigger path A (i.e. " << *hltPathA_name_ptr << ")" << std::endl;

				std::cout << "  ->For the highest Et ele in the event..." << std::endl;
				std::cout << "       charge=" << mcEles_HighestEt_charge << "; PDGid=" << mcEles_HighestEt_PDGid << "; status=" << mcEles_HighestEt_status << std::endl;
				std::cout << "       P=" << mcEles_HighestEt_OLDp4_ptr->P() << "=" << mcEles_HighestEt_p4.P() << "?; pt=" << mcEles_HighestEt_pt << "=" << mcEles_HighestEt_OLDp4_ptr->Pt() << "=" << mcEles_HighestEt_p4.Pt() <<"?" << std::endl;
				std::cout << "       eta=" << mcEles_HighestEt_OLDp4_ptr->Eta() << "=" << mcEles_HighestEt_p4.Eta() << "?; phi=" << mcEles_HighestEt_OLDp4_ptr->Phi() << "=" << mcEles_HighestEt_p4.Phi() << "?" << std::endl;
	
				std::cout << "  ->For the 2nd highest Et ele in the event..." << std::endl;
				std::cout << "       charge=" << mcEles_2ndHighestEt_charge << "; PDGid=" << mcEles_2ndHighestEt_PDGid << "; status=" << mcEles_2ndHighestEt_status << std::endl;
				std::cout << "       P=" << mcEles_2ndHighestEt_OLDp4_ptr->P() << "=" << mcEles_2ndHighestEt_p4.P() << "?; pt=" << mcEles_2ndHighestEt_pt << "=" << mcEles_2ndHighestEt_OLDp4_ptr->Pt() << "=" << mcEles_2ndHighestEt_p4.Pt() <<"?" << std::endl;
				std::cout << "       eta=" << mcEles_2ndHighestEt_OLDp4_ptr->Eta() << "=" << mcEles_2ndHighestEt_p4.Eta() << "?; phi=" << mcEles_2ndHighestEt_OLDp4_ptr->Phi() << "=" << mcEles_2ndHighestEt_p4.Phi() << "?" << std::endl;
	
				std::cout << "  ->For the Z candidate formed from the highest two Et final state electrons in this event..." << std::endl;
				std::cout << "       P=" << mcZcandidate_OLDp4_ptr->P() << "=" << mcZcandidate_p4.P() << "?; Pt=" << mcZcandidate_pt << "=" << mcZcandidate_p4.Pt() << "; inv mass=" << mcZcandidate_p4.M() << std::endl;
				std::cout << "       eta=" << mcZcandidate_eta << "=" << mcZcandidate_p4.Eta() << "?; phi=" << mcZcandidate_phi << "=" << mcZcandidate_p4.Phi() << std::endl;
				std::cout << "       dEta=" << mcZcandidate_dEtaEles << "; dPhi=" << mcZcandidate_dPhiEles << "; dR=" << mcZcandidate_dREles << std::endl;
				std::cout << "       lab frame opening angle = " << mcZcandidate_openingAngle << "=" << mcEles_HighestEt_p4.Angle(mcEles_2ndHighestEt_p4.Vect()) << std::endl;
			}

			//Filling histograms for highest Et ele quantities...
			mcHighestEtEle_PtHist->Fill(mcEles_HighestEt_pt,evtWeight);
			mcHighestEtEle_LogPtHist->Fill(mcEles_HighestEt_pt,evtWeight);
			mcHighestEtEle_EtaHist->Fill(mcEles_HighestEt_p4.Eta(),evtWeight);
			mcHighestEtEle_PhiHist->Fill(mcEles_HighestEt_p4.Phi(),evtWeight);
			if(hltPathA_decision){
				mcHighestEtEle_PostTrigPtHist->Fill(mcEles_HighestEt_pt,evtWeight);
				mcHighestEtEle_PostTrigLogPtHist->Fill(mcEles_HighestEt_pt,evtWeight);
				mcHighestEtEle_PostTrigEtaHist->Fill(mcEles_HighestEt_p4.Eta(),evtWeight);
				mcHighestEtEle_PostTrigPhiHist->Fill(mcEles_HighestEt_p4.Phi(),evtWeight);
			}
			//Filling histograms for 2nd highest Et ele quantities...
			mc2ndHighestEtEle_PtHist->Fill(   mcEles_2ndHighestEt_pt,      evtWeight);
			mc2ndHighestEtEle_LogPtHist->Fill(mcEles_2ndHighestEt_pt,      evtWeight);
			mc2ndHighestEtEle_EtaHist->Fill(  mcEles_2ndHighestEt_p4.Eta(),evtWeight);
			mc2ndHighestEtEle_PhiHist->Fill(  mcEles_2ndHighestEt_p4.Phi(),evtWeight);
			if(hltPathA_decision){
				mc2ndHighestEtEle_PostTrigPtHist->Fill(   mcEles_2ndHighestEt_pt,     evtWeight);
				mc2ndHighestEtEle_PostTrigLogPtHist->Fill(mcEles_2ndHighestEt_pt,     evtWeight);
				mc2ndHighestEtEle_PostTrigEtaHist->Fill(  mcEles_2ndHighestEt_p4.Eta(),evtWeight);
				mc2ndHighestEtEle_PostTrigPhiHist->Fill(  mcEles_2ndHighestEt_p4.Phi(),evtWeight);
			}

			mcZcandidate_PtHist->Fill(mcZcandidate_pt,evtWeight);
			mcZcandidate_LogPtHist->Fill(mcZcandidate_pt,evtWeight);
			mcZcandidate_EtaHist->Fill(mcZcandidate_eta,evtWeight);
			mcZcandidate_PhiHist->Fill(mcZcandidate_phi,evtWeight);
			mcZcandidate_MassHist->Fill(mcZcandidate_mass,evtWeight);
			mcZcandidate_dEtaHist->Fill(mcZcandidate_dEtaEles,evtWeight);
			mcZcandidate_dPhiHist->Fill(mcZcandidate_dPhiEles,evtWeight);
			mcZcandidate_dRHist->Fill(  mcZcandidate_dREles,evtWeight);
			mcZcandidate_openingAngleHist->Fill(mcZcandidate_openingAngle,evtWeight);
			if(hltPathA_decision){
				mcZcandidate_PostTrigPtHist->Fill(mcZcandidate_pt,evtWeight);
				mcZcandidate_PostTrigLogPtHist->Fill(mcZcandidate_pt,evtWeight);
				mcZcandidate_PostTrigEtaHist->Fill(mcZcandidate_eta,evtWeight);
				mcZcandidate_PostTrigPhiHist->Fill(mcZcandidate_phi,evtWeight);
				mcZcandidate_PostTrigMassHist->Fill(mcZcandidate_mass,evtWeight);
				mcZcandidate_PostTrigdEtaHist->Fill(mcZcandidate_dEtaEles,evtWeight);
				mcZcandidate_PostTrigdPhiHist->Fill(mcZcandidate_dPhiEles,evtWeight);
				mcZcandidate_PostTrigdRHist->Fill(  mcZcandidate_dREles,evtWeight);
				mcZcandidate_PostTrigopeningAngleHist->Fill(mcZcandidate_openingAngle,evtWeight);
			}
			if(vFlg>1){std::cout << std::endl;}
		}

		std::cout << "   done." << std::endl;
	
		//Calculating errors on historgrams...
		SetPoissonErrorsOnHist(mcHighestEtEle_PtHist, evtWeight);
		SetPoissonErrorsOnHist(mcHighestEtEle_LogPtHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mcHighestEtEle_PostTrigPtHist, mcHighestEtEle_PtHist);
		//SetBinomialErrorsOnNumerator(mcHighestEtEle_PostTrigLogPtHist,mcHighestEtEle_LogPtHist);
		SetPoissonErrorsOnHist(mcHighestEtEle_EtaHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mcHighestEtEle_PostTrigEtaHist, mcHighestEtEle_EtaHist);
		SetPoissonErrorsOnHist(mcHighestEtEle_PhiHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mcHighestEtEle_PostTrigPhiHist, mcHighestEtEle_PhiHist);
	
		SetPoissonErrorsOnHist(mc2ndHighestEtEle_PtHist, evtWeight);
		SetPoissonErrorsOnHist(mc2ndHighestEtEle_LogPtHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mc2ndHighestEtEle_PostTrigPtHist, mc2ndHighestEtEle_PtHist);
		//SetBinomialErrorsOnNumerator(mc2ndHighestEtEle_PostTrigLogPtHist, mc2ndHighestEtEle_LogPtHist);
		SetPoissonErrorsOnHist(mc2ndHighestEtEle_EtaHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mc2ndHighestEtEle_PostTrigEtaHist, mc2ndHighestEtEle_EtaHist);
		SetPoissonErrorsOnHist(mc2ndHighestEtEle_PhiHist, evtWeight);
		//SetBinomialErrorsOnNumerator(mc2ndHighestEtEle_PostTrigPhiHist, mc2ndHighestEtEle_PhiHist);

		SetPoissonErrorsOnHist(mcZcandidate_PtHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_LogPtHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_EtaHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_PhiHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_MassHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_dEtaHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_dPhiHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_dRHist, evtWeight);
		SetPoissonErrorsOnHist(mcZcandidate_openingAngleHist, evtWeight);

		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigPtHist,    mcZcandidate_PtHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigLogPtHist, mcZcandidate_LogPtHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigEtaHist,   mcZcandidate_EtaHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigPhiHist,   mcZcandidate_PhiHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigMassHist,  mcZcandidate_MassHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigdEtaHist,  mcZcandidate_dEtaHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigdPhiHist,  mcZcandidate_dPhiHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigdRHist,    mcZcandidate_dRHist);
		//SetBinomialErrorsOnNumerator(mcZcandidate_PostTrigopeningAngleHist, mcZcandidate_openingAngleHist);
	

		//Calculating values and errors (binomial) for TrigEffi hists...
		CalcPreWeightEffiHist(mcHighestEtEle_PtTrigEffi,   mcHighestEtEle_PostTrigPtHist,    mcHighestEtEle_PtHist,    evtWeight);
		CalcPreWeightEffiHist(mcHighestEtEle_LogPtTrigEffi,mcHighestEtEle_PostTrigLogPtHist, mcHighestEtEle_LogPtHist, evtWeight);
		CalcPreWeightEffiHist(mcHighestEtEle_EtaTrigEffi,  mcHighestEtEle_PostTrigEtaHist,   mcHighestEtEle_EtaHist,   evtWeight);
		CalcPreWeightEffiHist(mcHighestEtEle_PhiTrigEffi,  mcHighestEtEle_PostTrigPhiHist,   mcHighestEtEle_PhiHist,   evtWeight);
	
		CalcPreWeightEffiHist(mc2ndHighestEtEle_PtTrigEffi,   mc2ndHighestEtEle_PostTrigPtHist,    mc2ndHighestEtEle_PtHist,    evtWeight);
		CalcPreWeightEffiHist(mc2ndHighestEtEle_LogPtTrigEffi,mc2ndHighestEtEle_PostTrigLogPtHist, mc2ndHighestEtEle_LogPtHist, evtWeight);
		CalcPreWeightEffiHist(mc2ndHighestEtEle_EtaTrigEffi,  mc2ndHighestEtEle_PostTrigEtaHist,   mc2ndHighestEtEle_EtaHist,   evtWeight);
		CalcPreWeightEffiHist(mc2ndHighestEtEle_PhiTrigEffi,  mc2ndHighestEtEle_PostTrigPhiHist,   mc2ndHighestEtEle_PhiHist,   evtWeight);

		CalcPreWeightEffiHist(mcZcandidate_PtTrigEffi,    mcZcandidate_PostTrigPtHist,    mcZcandidate_PtHist,    evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_LogPtTrigEffi, mcZcandidate_PostTrigLogPtHist, mcZcandidate_LogPtHist, evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_EtaTrigEffi,   mcZcandidate_PostTrigEtaHist,   mcZcandidate_EtaHist,   evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_PhiTrigEffi,   mcZcandidate_PostTrigPhiHist,   mcZcandidate_PhiHist,   evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_MassTrigEffi,  mcZcandidate_PostTrigMassHist,  mcZcandidate_MassHist,  evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_dEtaTrigEffi,  mcZcandidate_PostTrigdEtaHist,  mcZcandidate_dEtaHist,  evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_dPhiTrigEffi,  mcZcandidate_PostTrigdPhiHist,  mcZcandidate_dPhiHist,  evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_dRTrigEffi,    mcZcandidate_PostTrigdRHist,    mcZcandidate_dRHist,    evtWeight);
		CalcPreWeightEffiHist(mcZcandidate_openingAngleTrigEffi, mcZcandidate_PostTrigopeningAngleHist, mcZcandidate_openingAngleHist, evtWeight);


		//Output file creation...
		f->Close();

		TFile f_histos(outputFileName,"RECREATE");
	
		f_histos.Write();
	
		mcHighestEtEle_PtHist->Write();
		mcHighestEtEle_LogPtHist->Write();
		mcHighestEtEle_EtaHist->Write();
		mcHighestEtEle_PhiHist->Write();
		mcHighestEtEle_PostTrigPtHist->Write();
		mcHighestEtEle_PostTrigLogPtHist->Write();
		mcHighestEtEle_PostTrigEtaHist->Write();
		mcHighestEtEle_PostTrigPhiHist->Write();
		mcHighestEtEle_PtTrigEffi->Write();
		mcHighestEtEle_LogPtTrigEffi->Write();
		mcHighestEtEle_EtaTrigEffi->Write();
		mcHighestEtEle_PhiTrigEffi->Write();

		mc2ndHighestEtEle_PtHist->Write();
		mc2ndHighestEtEle_LogPtHist->Write();
		mc2ndHighestEtEle_EtaHist->Write();
		mc2ndHighestEtEle_PhiHist->Write();
		mc2ndHighestEtEle_PostTrigPtHist->Write();
		mc2ndHighestEtEle_PostTrigLogPtHist->Write();
		mc2ndHighestEtEle_PostTrigEtaHist->Write();
		mc2ndHighestEtEle_PostTrigPhiHist->Write();

		mcZcandidate_PtHist->Write();
		mcZcandidate_LogPtHist->Write();
		mcZcandidate_EtaHist->Write();
		mcZcandidate_PhiHist->Write();
		mcZcandidate_MassHist->Write();
		mcZcandidate_dEtaHist->Write();
		mcZcandidate_dPhiHist->Write();
		mcZcandidate_dRHist->Write();
		mcZcandidate_openingAngleHist->Write();

		mcZcandidate_PostTrigPtHist->Write();
		mcZcandidate_PostTrigLogPtHist->Write();
		mcZcandidate_PostTrigEtaHist->Write();
		mcZcandidate_PostTrigPhiHist->Write();
		mcZcandidate_PostTrigMassHist->Write();
		mcZcandidate_PostTrigdEtaHist->Write();
		mcZcandidate_PostTrigdPhiHist->Write();
		mcZcandidate_PostTrigdRHist->Write();
		mcZcandidate_PostTrigopeningAngleHist->Write();
	
		mcZcandidate_PtTrigEffi->Write();
		mcZcandidate_LogPtTrigEffi->Write();
		mcZcandidate_EtaTrigEffi->Write();
		mcZcandidate_PhiTrigEffi->Write();
		mcZcandidate_MassTrigEffi->Write();
		mcZcandidate_dEtaTrigEffi->Write();
		mcZcandidate_dPhiTrigEffi->Write();
		mcZcandidate_dRTrigEffi->Write();
		mcZcandidate_openingAngleTrigEffi->Write();
	
		f_histos.Close();*/
	}
	
	return 1;
}


