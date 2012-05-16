
// C++ includes
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

// ROOT includes

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
#include "TLorentzVector.h"
//Must include TROOT.h in order for the gROOT->ProcessLine("#inlcude <vector>") line to work.
#include "TROOT.h" 
//    gROOT->ProcessLine("#include <vector>");


// BstdZee includes
#include "BstdZeeFirst/Analyser/interface/tswEventHelper.h"
#include "BstdZeeFirst/Analyser/interface/tswHEEPEle.h"
#include "BstdZeeFirst/Analyser/interface/tswReconValidationHistos.h"
#include "BstdZeeFirst/Analyser/interface/tswUsefulFunctions.h"
#include "BstdZeeFirst/Analyser/interface/tswHEEPDiEle.h"
#include "BstdZeeFirst/Analyser/interface/tswDiEleDistns.h"
#include "BstdZeeFirst/Analyser/interface/tswDiEleDistnsByRgn.h"
#include "BstdZeeFirst/Analyser/interface/tswMCParticle.h"
#include "BstdZeeFirst/Analyser/interface/tswMCParticleDistns.h"
#include "BstdZeeFirst/Analyser/interface/tswMuonDistns.h"

#include "BstdZeeFirst/Analyser/interface/tswEleMuObject.h"
#include "BstdZeeFirst/Analyser/interface/tswEleMuDistns.h"


#ifdef __MAKECINT__
#pragma link C++ class vector<Int_t>+;
#pragma link C++ class vector< vector<float> >;
#endif

//Functions...


//=====================================================================================================
//-----------------------------------------------------------------------------------------------------
namespace tsw{
	class ABCDMethodTree{
	private:
		// PRIVATE MEMBERS
		TTree* abcdTree_;

		// For the event & di-ele branches ...
		Bool_t treeVar_trgDecision_;
		std::string* treeVar_trgNamePtr_; std::string treeVar_trgName_;
		Double_t treeVar_weight_;

		TLorentzVector* treeVar_p4Ptr_; TLorentzVector treeVar_p4_;
		Double_t treeVar_mass_;
		Double_t treeVar_pT_;
		Bool_t treeVar_passHEEPNoIso_;
		Bool_t treeVar_passModHEEPIso_;
		Bool_t treeVar_passModTrkIso_;
		Bool_t treeVar_passModEmHad1Iso_;

		// For the eleA & eleB branches ...
		TLorentzVector* treeVar_eleA_p4Ptr_; TLorentzVector treeVar_eleA_p4_;
		Int_t treeVar_eleA_charge_;
		Double_t treeVar_eleA_EOverP_;
		Bool_t treeVar_eleA_passHEEPNoIso_;

		TLorentzVector* treeVar_eleB_p4Ptr_; TLorentzVector treeVar_eleB_p4_;
		Int_t treeVar_eleB_charge_;
		Double_t treeVar_eleB_EOverP_;
		Bool_t treeVar_eleB_passHEEPNoIso_;

	public:
		ABCDMethodTree(){
			abcdTree_ = new TTree("abcdTree","ABCD method di-ele data");
			abcdTree_->SetDirectory(0); // This line is needed as a 'QUICK FIX' to stop the following error when running over very large nos. of events ...
			/* Error is as follows:
			 * Error in <TTree::Fill>: Failed filling branch:myTree.mass, nbytes=-1, entry=3990
			 *  This error is symptomatic of a Tree created as a memory-resident Tree
			 *  Instead of doing:
			 *     TTree *T = new TTree(...)
			 *     TFile *f = new TFile(...)
			 *  you should do:
			 *     TFile *f = new TFile(...)
			 *     TTree *T = new TTree(...)
			 */
			// Setting up the event / di-ele branches ...
			abcdTree_->Branch("trgDecision", &treeVar_trgDecision_, "trgDecision/O");
			abcdTree_->Branch("trgName",     &treeVar_trgNamePtr_);
			treeVar_trgNamePtr_ = &treeVar_trgName_;
			abcdTree_->Branch("weight",      &treeVar_weight_, "weight/D");

			abcdTree_->Branch("diEle_p4",   &treeVar_p4Ptr_);
			treeVar_p4Ptr_ = &treeVar_p4_;
			abcdTree_->Branch("diEle_mass", &treeVar_mass_,   "mass/D");
			abcdTree_->Branch("diEle_pT",   &treeVar_pT_,     "pT/D");
			abcdTree_->Branch("diEle_passHEEPNoIso",  &treeVar_passHEEPNoIso_,  "diEle_passHEEPNoIso/O");
			abcdTree_->Branch("diEle_passModHEEPIso", &treeVar_passModHEEPIso_, "diEle_passModHEEPIso/O");
			abcdTree_->Branch("diEle_passModTrkIso", &treeVar_passModTrkIso_, "diEle_passModTrkIso/O");
			abcdTree_->Branch("diEle_passModEmHad1Iso", &treeVar_passModEmHad1Iso_, "diEle_passModEmHad1Iso/O");

			// Setting up the branches for each electron ...
			abcdTree_->Branch("eleA_p4",     &treeVar_eleA_p4Ptr_);
			treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;
			abcdTree_->Branch("eleA_charge",   &treeVar_eleA_charge_,   "eleA_charge/I"); //Int_t
			abcdTree_->Branch("eleA_EOverP", &treeVar_eleA_EOverP_, "eleA_EOverP/D"); //Double_t
			abcdTree_->Branch("eleA_passHEEPNoIso", &treeVar_eleA_passHEEPNoIso_, "eleA_passHEEPNoIso/O");

			abcdTree_->Branch("eleB_p4",     &treeVar_eleB_p4Ptr_);
			treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;
			abcdTree_->Branch("eleB_charge",   &treeVar_eleB_charge_,   "eleB_charge/I"); //Int_t
			abcdTree_->Branch("eleB_EOverP", &treeVar_eleB_EOverP_, "eleB_EOverP/D"); //Double_t
			abcdTree_->Branch("eleB_passHEEPNoIso", &treeVar_eleB_passHEEPNoIso_, "eleB_passHEEPNoIso/O");
		}
		~ABCDMethodTree(){ delete abcdTree_; }

		void FillTree(tsw::HEEPDiEle* diEle, std::string trgName, bool trgDecision, double evtWeight){
			// Setting values of variables for eleA & eleB branches ...
			tsw::HEEPEle diEle_eleA = diEle->eleA();
			treeVar_eleA_p4_ = diEle_eleA.p4();
			treeVar_eleA_charge_ = diEle_eleA.charge();
			treeVar_eleA_EOverP_ = diEle_eleA.epIn();
			treeVar_eleA_passHEEPNoIso_ = diEle_eleA.ApplyHEEPCutsNoIso();

			tsw::HEEPEle diEle_eleB = diEle->eleB();
			treeVar_eleB_p4_ = diEle_eleB.p4();
			treeVar_eleB_charge_ = diEle_eleB.charge();
			treeVar_eleB_EOverP_ = diEle_eleB.epIn();
			treeVar_eleB_passHEEPNoIso_ = diEle_eleB.ApplyHEEPCutsNoIso();

			// Setting values of variables for event & di-ele branches ...
			treeVar_trgDecision_ = trgDecision;
			treeVar_trgName_ = trgName;
			treeVar_weight_ = evtWeight;

			treeVar_p4_ = diEle->totalP4();
			treeVar_mass_ = diEle->invMass();
			treeVar_pT_ = diEle->pT();
			treeVar_passHEEPNoIso_ = (treeVar_eleA_passHEEPNoIso_ && treeVar_eleB_passHEEPNoIso_);
			treeVar_passModTrkIso_ = diEle->ApplyDiEleTrkIsolCut();
			treeVar_passModEmHad1Iso_ = diEle->ApplyDiEleEmHad1IsolCut();
			treeVar_passModHEEPIso_ = (treeVar_passModTrkIso_ && treeVar_passModEmHad1Iso_);

			// And finally fill the tree ...
			abcdTree_->Fill();
		}
		void SaveToFile(TString outFileName){
			TFile f_tree(outFileName,"RECREATE");
			abcdTree_->Write();
			f_tree.Close();
		}
	};

	///////////////////////////////////////////////////////////////////////////////
	// EffiCalcTree: Simple class for generating trees containing small amount of selection cut & kinematic
	//               info about Z candidates (plus event information)
	//
	class DiEleTree{
	private:
		// PRIVATE MEMBERS
		TTree* diEleTree_;

		// For the event & kinematic branches ...
		Double_t treeVar_weight_;
		Double_t treeVar_genWeight_;
		Double_t treeVar_puWeight_;
		Double_t treeVar_xsecWeight_;

		UInt_t treeVar_runNum_;
		UInt_t treeVar_lumiNum_;
		UInt_t treeVar_evtNum_;

		Bool_t treeVar_trgDecision_;

		TLorentzVector* treeVar_Zp4Ptr_; TLorentzVector treeVar_Zp4_;
		Double_t treeVar_ZpT_;
		Double_t treeVar_Zmass_;
		Double_t treeVar_dR_;
		Double_t treeVar_dEta_;
		Double_t treeVar_dPhi_;

		TLorentzVector* treeVar_eleA_p4Ptr_; TLorentzVector treeVar_eleA_p4_;
		TLorentzVector* treeVar_eleB_p4Ptr_; TLorentzVector treeVar_eleB_p4_;
		Int_t treeVar_eleA_modHeepCutCode_;
		Int_t treeVar_eleB_modHeepCutCode_;

	public:
		DiEleTree(std::string treeName="zBosonTree"){
			diEleTree_ = new TTree(treeName.c_str(), "Tree of Z candidates");
			diEleTree_->SetDirectory(0); // This line is needed as a 'QUICK FIX' to stop the following error when running over very large nos. of events ...
			/* Error is as follows:
			 * Error in <TTree::Fill>: Failed filling branch:myTree.mass, nbytes=-1, entry=3990
			 *  This error is symptomatic of a Tree created as a memory-resident Tree
			 *  Instead of doing:
			 *     TTree *T = new TTree(...)
			 *     TFile *f = new TFile(...)
			 *  you should do:
			 *     TFile *f = new TFile(...)
			 *     TTree *T = new TTree(...)
			 */
			// Setting up the event / di-ele branches ...
			diEleTree_->Branch("weight",     &treeVar_weight_,     "weight/D");
			diEleTree_->Branch("genWeight",  &treeVar_genWeight_,  "genWeight/D");
			diEleTree_->Branch("puWeight",   &treeVar_puWeight_,   "puWeight/D");
			diEleTree_->Branch("xsecWeight", &treeVar_xsecWeight_, "xsecWeight/D");

			diEleTree_->Branch("run",    &treeVar_runNum_,   "run/i");
			diEleTree_->Branch("lumi",   &treeVar_lumiNum_,   "lumi/i");
			diEleTree_->Branch("evtNum", &treeVar_evtNum_,   "evtNum/i");

			diEleTree_->Branch("trgDecision", &treeVar_trgDecision_, "trgDecision/O");

			diEleTree_->Branch("Zp4",   &treeVar_Zp4Ptr_); treeVar_Zp4Ptr_ = &treeVar_Zp4_;
			diEleTree_->Branch("ZpT",   &treeVar_ZpT_,     "ZpT/D");
			diEleTree_->Branch("Zmass", &treeVar_Zmass_,   "Zmass/D");
			diEleTree_->Branch("dR",    &treeVar_dR_,      "dR/D");
			diEleTree_->Branch("dEta",  &treeVar_dEta_,    "dEta/D");
			diEleTree_->Branch("dPhi",  &treeVar_dPhi_,    "dPhi/D");

			diEleTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);  treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;
			diEleTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);	treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;

			diEleTree_->Branch("eleA_modHeepCutCode", &treeVar_eleA_modHeepCutCode_, "eleA_modHeepCutCode/I");
			diEleTree_->Branch("eleB_modHeepCutCode", &treeVar_eleB_modHeepCutCode_, "eleB_modHeepCutCode/I");


		}
		~DiEleTree(){ delete diEleTree_; }

		void FillTree(const tsw::HEEPDiEle& diEle, const tsw::EventHelper& evtHelper, bool trigDecision)
		{
			treeVar_weight_     = evtHelper.totWeight();
			treeVar_genWeight_  = evtHelper.genWeight();
			treeVar_puWeight_   = evtHelper.puWeight();
			treeVar_xsecWeight_ = evtHelper.xsecWeight();

			treeVar_runNum_ = evtHelper.runNum();
			treeVar_lumiNum_ = evtHelper.lumiSec();
			treeVar_evtNum_ = evtHelper.eventNum();

			treeVar_trgDecision_ = trigDecision;

			treeVar_Zp4_   = diEle.p4();
			treeVar_ZpT_   = diEle.pT();
			treeVar_Zmass_ = diEle.invMass();
			treeVar_dR_    = diEle.deltaR();
			treeVar_dEta_  = diEle.deltaEta();
			treeVar_dPhi_  = diEle.deltaPhi();

			treeVar_eleA_p4_ = diEle.eleA().p4();
			treeVar_eleB_p4_ = diEle.eleB().p4();

			treeVar_eleA_modHeepCutCode_ = diEle.eleA().heepIdModIsoCutCode(evtHelper);
			treeVar_eleB_modHeepCutCode_ = diEle.eleB().heepIdModIsoCutCode(evtHelper);

			// And finally fill the tree ...
			diEleTree_->Fill();
		}
		void SaveToFile(TString outFileName){
			TFile f_tree(outFileName,"RECREATE");
			diEleTree_->Write();
			f_tree.Close();
		}
	};

	///////////////////////////////////////////////////////////////////////////////
	// EffiCalcTree: Simple class used to generate trees that are used in calculating
	//               the efficiency of cuts on reconstruction of Z bosons
	//
	class EffiCalcTree{
	private:
		// PRIVATE MEMBERS
		TTree* effiCalcTree_;

		// Member variables for storing branch information (all are RECO-level unless otherwise specified in name) ...
		unsigned int treeVar_runNum_;
		unsigned int treeVar_lumiSec_;
		unsigned int treeVar_evtNum_;

		Double_t treeVar_weight_;
		float treeVar_mc_numVtx_;

		TLorentzVector* treeVar_mcZ_ele1_p4Ptr_; TLorentzVector treeVar_mcZ_ele1_p4_;
		TLorentzVector* treeVar_mcZ_ele2_p4Ptr_; TLorentzVector treeVar_mcZ_ele2_p4_;
		Bool_t treeVar_ptAcc_;
		Bool_t treeVar_ebebAcceptance_;
		Bool_t treeVar_ebeeAcceptance_;
		Bool_t treeVar_eeeeAcceptance_;
		Bool_t treeVar_bothRecod_;

		Double_t treeVar_ZpT_;
		Double_t treeVar_ZdEta_;
		Double_t treeVar_ZdPhi_;
		Double_t treeVar_ZdR_;
		TLorentzVector* treeVar_eleA_p4Ptr_; TLorentzVector treeVar_eleA_p4_;
		Double_t treeVar_eleA_dRmc_;
		TLorentzVector* treeVar_eleB_p4Ptr_;	TLorentzVector treeVar_eleB_p4_;
		Double_t treeVar_eleB_dRmc_;

		Bool_t treeVar_cut_both_fiducial_;
		Bool_t treeVar_cut_both_ecalDriven_;
		Bool_t treeVar_cut_both_dEta_;
		Bool_t treeVar_cut_both_dPhi_;
		Bool_t treeVar_cut_both_hOverE_;
		Bool_t treeVar_cut_both_showerShape_;
		Bool_t treeVar_cut_both_heepId_;

		Bool_t treeVar_cut_eleA_stdTrkIso_;
		Bool_t treeVar_cut_eleB_stdTrkIso_;
		Double_t treeVar_eleA_stdTrkIso_;
		Double_t treeVar_eleB_stdTrkIso_;
		Bool_t treeVar_cut_eleA_stdEmH1Iso_;
		Bool_t treeVar_cut_eleB_stdEmH1Iso_;
		Double_t treeVar_eleA_stdEmH1Iso_;
		Double_t treeVar_eleB_stdEmH1Iso_;

		unsigned int treeVar_eleA_nTrksInnerVeto_;
		unsigned int treeVar_eleB_nTrksInnerVeto_;

		Bool_t treeVar_cut_eleA_modTrkIso_;
		Bool_t treeVar_cut_eleB_modTrkIso_;
		Double_t treeVar_eleA_modTrkIso_;
		Double_t treeVar_eleB_modTrkIso_;
		Bool_t treeVar_cut_eleA_scModEmH1Iso_;
		Bool_t treeVar_cut_eleB_scModEmH1Iso_;
		Double_t treeVar_eleA_scModEmH1Iso_;
		Double_t treeVar_eleB_scModEmH1Iso_;

		Double_t treeVar_combThr_EmH1_;

		Double_t treeVar_eleA_isoDep_stdTrk_;
		Double_t treeVar_eleB_isoDep_stdTrk_;
		Double_t treeVar_eleA_isoDep_stdEmH1_;
		Double_t treeVar_eleB_isoDep_stdEmH1_;
		Double_t treeVar_cut_eleA_isoDep_stdEmH1_;
		Double_t treeVar_cut_eleB_isoDep_stdEmH1_;
		Double_t treeVar_eleA_isoDep_inrVetoModTrk_;
		Double_t treeVar_eleB_isoDep_inrVetoModTrk_;
		Double_t treeVar_eleA_isoDep_inrVetoModEmH1_;
		Double_t treeVar_eleB_isoDep_inrVetoModEmH1_;
		Bool_t treeVar_cut_eleA_isoDep_inrVetoModEmH1_;
		Bool_t treeVar_cut_eleB_isoDep_inrVetoModEmH1_;

		Double_t treeVar_eleA_inrXSVetoModTrk_;
		Double_t treeVar_eleA_inrSVetoModTrk_;
		Double_t treeVar_eleA_inrMVetoModTrk_;
		Double_t treeVar_eleA_inrLVetoModTrk_;
		Double_t treeVar_eleA_inrXLVetoModTrk_;

		Double_t treeVar_eleB_inrXSVetoModTrk_;
		Double_t treeVar_eleB_inrSVetoModTrk_;
		Double_t treeVar_eleB_inrMVetoModTrk_;
		Double_t treeVar_eleB_inrLVetoModTrk_;
		Double_t treeVar_eleB_inrXLVetoModTrk_;

		Double_t treeVar_eleA_inrXSVetoModEmH1_;
		Double_t treeVar_eleA_inrSVetoModEmH1_;
		Double_t treeVar_eleA_inrMVetoModEmH1_;
		Double_t treeVar_eleA_inrLVetoModEmH1_;
		Double_t treeVar_eleA_inrXLVetoModEmH1_;

		Double_t treeVar_eleB_inrXSVetoModEmH1_;
		Double_t treeVar_eleB_inrSVetoModEmH1_;
		Double_t treeVar_eleB_inrMVetoModEmH1_;
		Double_t treeVar_eleB_inrLVetoModEmH1_;
		Double_t treeVar_eleB_inrXLVetoModEmH1_;

		Bool_t treeVar_cut_eleA_inrMVetoModEmH1_;
		Bool_t treeVar_cut_eleB_inrMVetoModEmH1_;

		UInt_t treeVar_eleA_nGenHadronsDr04_;
		UInt_t treeVar_eleB_nGenHadronsDr04_;
		Double_t treeVar_eleA_ptSumGenHadronsDr04_;
		Double_t treeVar_eleB_ptSumGenHadronsDr04_;

		Double_t treeVar_eleA_EmH1RhoCorrn_;
		Double_t treeVar_eleB_EmH1RhoCorrn_;

		Int_t treeVar_eleA_heepIdModIsoCutCode_;
		Int_t treeVar_eleB_heepIdModIsoCutCode_;

	public:
		EffiCalcTree()
		{
			effiCalcTree_ = new TTree("zBosonEffiTree","Tree of Z candidate information for signal MC effi calc'ns");
			effiCalcTree_->SetDirectory(0); // This line is needed as a 'QUICK FIX' to stop the following error when running over very large nos. of events ...
			/* Error is as follows:
			 * Error in <TTree::Fill>: Failed filling branch:myTree.mass, nbytes=-1, entry=3990
			 *  This error is symptomatic of a Tree created as a memory-resident Tree
			 *  Instead of doing:
			 *     TTree *T = new TTree(...)
			 *     TFile *f = new TFile(...)
			 *  you should do:
			 *     TFile *f = new TFile(...)
			 *     TTree *T = new TTree(...)
			 */
			// Setting up the event / di-ele branches ...
			effiCalcTree_->Branch("weight",      &treeVar_weight_, "weight/D");
			effiCalcTree_->Branch("mc_numVtx", &treeVar_mc_numVtx_, "mc_numVtx/f");
			effiCalcTree_->Branch("run", &treeVar_runNum_, "run/i");
			effiCalcTree_->Branch("lumi", &treeVar_lumiSec_, "lumi/i");
			effiCalcTree_->Branch("evtNum", &treeVar_evtNum_, "evtNum/i");

			effiCalcTree_->Branch("mcZ_ele1_p4", &treeVar_mcZ_ele1_p4Ptr_);	treeVar_mcZ_ele1_p4Ptr_ = &treeVar_mcZ_ele1_p4_;
			effiCalcTree_->Branch("mcZ_ele2_p4", &treeVar_mcZ_ele2_p4Ptr_);	treeVar_mcZ_ele2_p4Ptr_ = &treeVar_mcZ_ele2_p4_;
			effiCalcTree_->Branch("mcAccept_pt",   &treeVar_ptAcc_, "mcAccept_pt/O");
			effiCalcTree_->Branch("mcAccept_ebeb", &treeVar_ebebAcceptance_, "mcAccept_ebeb/O");
			effiCalcTree_->Branch("mcAccept_ebee", &treeVar_ebeeAcceptance_, "mcAccept_ebee/O");
			effiCalcTree_->Branch("mcAccept_eeee", &treeVar_eeeeAcceptance_, "mcAccept_eeee/O");

			effiCalcTree_->Branch("bothRecod", &treeVar_bothRecod_, "bothRecod/O");

			effiCalcTree_->Branch("ZpT", &treeVar_ZpT_,     "ZpT/D");
			effiCalcTree_->Branch("ZdEta", &treeVar_ZdEta_,     "ZdEta/D");
			effiCalcTree_->Branch("ZdPhi", &treeVar_ZdPhi_,     "ZdPhi/D");
			effiCalcTree_->Branch("ZdR", &treeVar_ZdR_, "ZdR/D");
			effiCalcTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);	treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;
			effiCalcTree_->Branch("eleA_dRmc", &treeVar_eleA_dRmc_, "eleA_dRmc/D");
			effiCalcTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);	treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;
			effiCalcTree_->Branch("eleB_dRmc", &treeVar_eleB_dRmc_, "eleB_dRmc/D");


			effiCalcTree_->Branch("cut_both_fiducial", &treeVar_cut_both_fiducial_, "cut_both_fiducial/O");
			effiCalcTree_->Branch("cut_both_ecalDriven", &treeVar_cut_both_ecalDriven_, "cut_both_ecalDriven/O");
			effiCalcTree_->Branch("cut_both_dEta", &treeVar_cut_both_dEta_, "cut_both_dEta/O");
			effiCalcTree_->Branch("cut_both_dPhi", &treeVar_cut_both_dPhi_, "cut_both_dPhi/O");
			effiCalcTree_->Branch("cut_both_hOverE", &treeVar_cut_both_hOverE_, "cut_both_hOverE/O");
			effiCalcTree_->Branch("cut_both_showerShape", &treeVar_cut_both_showerShape_, "cut_both_showerShape/O");
			effiCalcTree_->Branch("cut_both_heepId", &treeVar_cut_both_heepId_, "cut_both_heepId/O");

			effiCalcTree_->Branch("cut_eleA_stdTrkIso", &treeVar_cut_eleA_stdTrkIso_, "cut_eleA_stdTrkIso/O");
			effiCalcTree_->Branch("cut_eleB_stdTrkIso", &treeVar_cut_eleB_stdTrkIso_, "cut_eleB_stdTrkIso/O");
			effiCalcTree_->Branch("eleA_stdTrkIso", &treeVar_eleA_stdTrkIso_, "eleA_stdTrkIso/D");
			effiCalcTree_->Branch("eleB_stdTrkIso", &treeVar_eleB_stdTrkIso_, "eleB_stdTrkIso/D");
			effiCalcTree_->Branch("cut_eleA_stdEmH1Iso", &treeVar_cut_eleA_stdEmH1Iso_, "cut_eleA_stdEmH1Iso/O");
			effiCalcTree_->Branch("cut_eleB_stdEmH1Iso", &treeVar_cut_eleB_stdEmH1Iso_, "cut_eleB_stdEmH1Iso/O");
			effiCalcTree_->Branch("eleA_stdEmH1Iso", &treeVar_eleA_stdEmH1Iso_, "eleA_stdEmH1Iso/D");
			effiCalcTree_->Branch("eleB_stdEmH1Iso", &treeVar_eleB_stdEmH1Iso_, "eleB_stdEmH1Iso/D");

			effiCalcTree_->Branch("eleA_nTrksInnerVeto", &treeVar_eleA_nTrksInnerVeto_, "eleA_nTrksInnerVeto/i");
			effiCalcTree_->Branch("eleB_nTrksInnerVeto", &treeVar_eleB_nTrksInnerVeto_, "eleB_nTrksInnerVeto/i");

			effiCalcTree_->Branch("cut_eleA_modTrkIso", &treeVar_cut_eleA_modTrkIso_, "cut_eleA_modTrkIso/O");
			effiCalcTree_->Branch("cut_eleB_modTrkIso", &treeVar_cut_eleB_modTrkIso_, "cut_eleB_modTrkIso/O");
			effiCalcTree_->Branch("eleA_modTrkIso", &treeVar_eleA_modTrkIso_, "eleA_modTrkIso/D");
			effiCalcTree_->Branch("eleB_modTrkIso", &treeVar_eleB_modTrkIso_, "eleB_modTrkIso/D");
			effiCalcTree_->Branch("cut_eleA_scModEmH1Iso", &treeVar_cut_eleA_scModEmH1Iso_, "cut_eleA_scModEmH1Iso/O");
			effiCalcTree_->Branch("cut_eleB_scModEmH1Iso", &treeVar_cut_eleB_scModEmH1Iso_, "cut_eleB_scModEmH1Iso/O");
			effiCalcTree_->Branch("eleA_scModEmH1Iso", &treeVar_eleA_scModEmH1Iso_, "eleA_scModEmH1Iso/D");
			effiCalcTree_->Branch("eleB_scModEmH1Iso", &treeVar_eleB_scModEmH1Iso_, "eleB_scModEmH1Iso/D");

			effiCalcTree_->Branch("combThr_EmH1", &treeVar_combThr_EmH1_, "combThr_EmH1/D");

			effiCalcTree_->Branch("eleA_isoDep_stdTrk",      &treeVar_eleA_isoDep_stdTrk_,  "eleA_isoDep_stdTrk/D");
			effiCalcTree_->Branch("eleB_isoDep_stdTrk",      &treeVar_eleB_isoDep_stdTrk_, "eleB_isoDep_stdTrk/D");
			effiCalcTree_->Branch("eleA_isoDep_stdEmH1",     &treeVar_eleA_isoDep_stdEmH1_, "eleA_isoDep_stdEmH1/D");
			effiCalcTree_->Branch("eleB_isoDep_stdEmH1",     &treeVar_eleB_isoDep_stdEmH1_, "eleB_isoDep_stdEmH1/D");
			effiCalcTree_->Branch("cut_eleA_isoDep_stdEmH1", &treeVar_cut_eleA_isoDep_stdEmH1_, "cut_eleA_isoDep_stdEmH1/D");
			effiCalcTree_->Branch("cut_eleB_isoDep_stdEmH1", &treeVar_cut_eleB_isoDep_stdEmH1_, "cut_eleB_isoDep_stdEmH1/D");
			effiCalcTree_->Branch("eleA_isoDep_inrVetoModTrk",     &treeVar_eleA_isoDep_inrVetoModTrk_, "eleA_isoDep_inrVetoModTrk/D");
			effiCalcTree_->Branch("eleB_isoDep_inrVetoModTrk",     &treeVar_eleB_isoDep_inrVetoModTrk_, "eleB_isoDep_inrVetoModTrk/D");
			effiCalcTree_->Branch("eleA_isoDep_inrVetoModEmH1",    &treeVar_eleA_isoDep_inrVetoModEmH1_, "eleA_isoDep_inrVetoModEmH1/D");
			effiCalcTree_->Branch("eleB_isoDep_inrVetoModEmH1",    &treeVar_eleB_isoDep_inrVetoModEmH1_, "eleB_isoDep_inrVetoModEmH1/D");
			effiCalcTree_->Branch("cut_eleA_isoDep_inrVetoModEmH1", &treeVar_cut_eleA_isoDep_inrVetoModEmH1_, "cut_eleA_isoDep_inrVetoModEmH1/O");
			effiCalcTree_->Branch("cut_eleB_isoDep_inrVetoModEmH1", &treeVar_cut_eleB_isoDep_inrVetoModEmH1_, "cut_eleB_isoDep_inrVetoModEmH1/O");

			effiCalcTree_->Branch("eleA_inrXSVetoModTrk", &treeVar_eleA_inrXSVetoModTrk_, "eleA_inrXSVetoModTrk/D");
			effiCalcTree_->Branch("eleA_inrSVetoModTrk",  &treeVar_eleA_inrSVetoModTrk_,  "eleA_inrSVetoModTrk/D");
			effiCalcTree_->Branch("eleA_inrMVetoModTrk",  &treeVar_eleA_inrMVetoModTrk_,  "eleA_inrMVetoModTrk/D");
			effiCalcTree_->Branch("eleA_inrLVetoModTrk",  &treeVar_eleA_inrLVetoModTrk_,  "eleA_inrLVetoModTrk/D");
			effiCalcTree_->Branch("eleA_inrXLVetoModTrk", &treeVar_eleA_inrXLVetoModTrk_, "eleA_inrXLVetoModTrk/D");

			effiCalcTree_->Branch("eleB_inrXSVetoModTrk", &treeVar_eleB_inrXSVetoModTrk_, "eleB_inrXSVetoModTrk/D");
			effiCalcTree_->Branch("eleB_inrSVetoModTrk",  &treeVar_eleB_inrSVetoModTrk_,  "eleB_inrSVetoModTrk/D");
			effiCalcTree_->Branch("eleB_inrMVetoModTrk",  &treeVar_eleB_inrMVetoModTrk_,  "eleB_inrMVetoModTrk/D");
			effiCalcTree_->Branch("eleB_inrLVetoModTrk",  &treeVar_eleB_inrLVetoModTrk_,  "eleB_inrLVetoModTrk/D");
			effiCalcTree_->Branch("eleB_inrXLVetoModTrk", &treeVar_eleB_inrXLVetoModTrk_, "eleB_inrXLVetoModTrk/D");

			effiCalcTree_->Branch("eleA_inrXSVetoModEmH1", &treeVar_eleA_inrXSVetoModEmH1_, "eleA_inrXSVetoModEmH1/D");
			effiCalcTree_->Branch("eleA_inrSVetoModEmH1",  &treeVar_eleA_inrSVetoModEmH1_,  "eleA_inrSVetoModEmH1/D");
			effiCalcTree_->Branch("eleA_inrMVetoModEmH1",  &treeVar_eleA_inrMVetoModEmH1_,  "eleA_inrMVetoModEmH1/D");
			effiCalcTree_->Branch("eleA_inrLVetoModEmH1",  &treeVar_eleA_inrLVetoModEmH1_,  "eleA_inrLVetoModEmH1/D");
			effiCalcTree_->Branch("eleA_inrXLVetoModEmH1", &treeVar_eleA_inrXLVetoModEmH1_, "eleA_inrXLVetoModEmH1/D");

			effiCalcTree_->Branch("eleB_inrXSVetoModEmH1", &treeVar_eleB_inrXSVetoModEmH1_, "eleB_inrXSVetoModEmH1/D");
			effiCalcTree_->Branch("eleB_inrSVetoModEmH1",  &treeVar_eleB_inrSVetoModEmH1_,  "eleB_inrSVetoModEmH1/D");
			effiCalcTree_->Branch("eleB_inrMVetoModEmH1",  &treeVar_eleB_inrMVetoModEmH1_,  "eleB_inrMVetoModEmH1/D");
			effiCalcTree_->Branch("eleB_inrLVetoModEmH1",  &treeVar_eleB_inrLVetoModEmH1_,  "eleB_inrLVetoModEmH1/D");
			effiCalcTree_->Branch("eleB_inrXLVetoModEmH1", &treeVar_eleB_inrXLVetoModEmH1_, "eleB_inrXLVetoModEmH1/D");

			effiCalcTree_->Branch("cut_eleA_inrMVetoModEmH1", &treeVar_cut_eleA_inrMVetoModEmH1_, "cut_eleA_inrMVetoModEmH1/O");
			effiCalcTree_->Branch("cut_eleB_inrMVetoModEmH1", &treeVar_cut_eleB_inrMVetoModEmH1_, "cut_eleB_inrMVetoModEmH1/O");

			effiCalcTree_->Branch("eleA_nGenHadronsDr04", &treeVar_eleA_nGenHadronsDr04_, "eleA_nGenHadronsDr04/i");
			effiCalcTree_->Branch("eleB_nGenHadronsDr04", &treeVar_eleB_nGenHadronsDr04_, "eleB_nGenHadronsDr04/i");
			effiCalcTree_->Branch("eleA_ptSumGenHadronsDr04", &treeVar_eleA_ptSumGenHadronsDr04_, "eleA_ptSumGenHadronsDr04/D");
			effiCalcTree_->Branch("eleB_ptSumGenHadronsDr04", &treeVar_eleB_ptSumGenHadronsDr04_, "eleB_ptSumGenHadronsDr04/D");

			effiCalcTree_->Branch("eleA_EmH1RhoCorrn", &treeVar_eleA_EmH1RhoCorrn_, "eleA_EmH1RhoCorrn_/D");
			effiCalcTree_->Branch("eleB_EmH1RhoCorrn", &treeVar_eleB_EmH1RhoCorrn_, "eleB_EmH1RhoCorrn_/D");

			effiCalcTree_->Branch("eleA_heepIdModIsoCutCode", &treeVar_eleA_heepIdModIsoCutCode_, "eleA_heepIdModIsoCutCode/I");
			effiCalcTree_->Branch("eleB_heepIdModIsoCutCode", &treeVar_eleB_heepIdModIsoCutCode_, "eleB_heepIdModIsoCutCode/I");
		}
		~EffiCalcTree(){ delete effiCalcTree_; }

		void FillTree(tsw::HEEPDiEle* diEle, const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2,
				const unsigned int runNum, const unsigned int lumiNum, const unsigned int evtNum, const tsw::EventHelper& eventHelper)
		{
			treeVar_weight_  = eventHelper.totWeight();
			treeVar_mc_numVtx_ = eventHelper.GetMCPU_nVtx();
			treeVar_runNum_  = runNum;
			treeVar_lumiSec_ = lumiNum;
			treeVar_evtNum_  = evtNum;

			treeVar_mcZ_ele1_p4_ = mcZboson_ele1;
			treeVar_mcZ_ele2_p4_ = mcZboson_ele2;

			bool mcZeles_zMass = fabs((mcZboson_ele1+mcZboson_ele2).M()-90.0)<30.0;
			bool mcZeles_pTacc = ( diEle->eleA().et()>35.0 && diEle->eleB().et()>35.0 );
			bool mcZele1_isEB = fabs(mcZboson_ele1.Eta())<1.4;
			bool mcZele2_isEB = fabs(mcZboson_ele2.Eta())<1.4;
			bool mcZele1_isEE = (fabs(mcZboson_ele1.Eta())>1.6 && fabs(mcZboson_ele1.Eta())<2.45);
			bool mcZele2_isEE = (fabs(mcZboson_ele2.Eta())>1.6 && fabs(mcZboson_ele2.Eta())<2.45);

			treeVar_ptAcc_ = mcZeles_pTacc;
			treeVar_ebebAcceptance_ = (mcZele1_isEB && mcZele2_isEB) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_ebeeAcceptance_ = ((mcZele1_isEB && mcZele2_isEE) || (mcZele1_isEE && mcZele2_isEB) ) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_eeeeAcceptance_ = (mcZele1_isEE && mcZele2_isEE) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_bothRecod_ = true;

			treeVar_ZpT_ = diEle->pT();
			treeVar_ZdEta_ = diEle->deltaEta();
			treeVar_ZdPhi_ = diEle->deltaPhi();
			treeVar_ZdR_   = diEle->deltaR();
			treeVar_eleA_p4_ = diEle->eleA().p4();
			treeVar_eleB_p4_ = diEle->eleB().p4();

			// Calculate dR between RECO & MC electrons ...
			const double dR_recoEleA_mcEle1 = treeVar_eleA_p4_.DeltaR(treeVar_mcZ_ele1_p4_);
			const double dR_recoEleB_mcEle2 = treeVar_eleB_p4_.DeltaR(treeVar_mcZ_ele2_p4_);
			const double dR_recoEleA_mcEle2 = treeVar_eleA_p4_.DeltaR(treeVar_mcZ_ele2_p4_);
			const double dR_recoEleB_mcEle1 = treeVar_eleB_p4_.DeltaR(treeVar_mcZ_ele1_p4_);
			if( std::min(dR_recoEleA_mcEle1,dR_recoEleB_mcEle2) <= std::min(dR_recoEleA_mcEle2,dR_recoEleB_mcEle1) ){
				treeVar_eleA_dRmc_ = dR_recoEleA_mcEle1;
				treeVar_eleB_dRmc_ = dR_recoEleB_mcEle2;
			}
			else{
				treeVar_eleA_dRmc_ = dR_recoEleA_mcEle2;
				treeVar_eleB_dRmc_ = dR_recoEleB_mcEle1;
			}

			const tsw::HEEPEle recoEleA = diEle->eleA();
			const tsw::HEEPEle recoEleB = diEle->eleB();
			treeVar_cut_both_fiducial_    = ( recoEleA.et()>35.0 && recoEleB.et()>35.0 )
															&& (fabs(recoEleA.scEta())<1.442 && fabs(recoEleB.scEta())<1.442 );
			treeVar_cut_both_ecalDriven_  = (recoEleA.isEcalDriven() && recoEleB.isEcalDriven() );
			treeVar_cut_both_dEta_        = ( fabs(recoEleA.dEtaIn())<0.005 && fabs(recoEleB.dEtaIn())<0.005 );
			treeVar_cut_both_dPhi_        = ( fabs(recoEleA.dPhiIn())<0.06  && fabs(recoEleB.dPhiIn())<0.06 );
			treeVar_cut_both_hOverE_      = ( recoEleA.hOverE()<0.05 && recoEleB.hOverE()<0.05 );
			treeVar_cut_both_showerShape_ = ( (recoEleA.e2x5MaxOver5x5()>0.94) || (recoEleA.e1x5Over5x5()>0.83) )
															&& ( (recoEleB.e2x5MaxOver5x5()>0.94) || (recoEleB.e1x5Over5x5()>0.83) );
			treeVar_cut_both_heepId_ = diEle->eleA().ApplyHEEPCutsNoIso() && diEle->eleB().ApplyHEEPCutsNoIso();

			treeVar_cut_eleA_stdTrkIso_ = diEle->eleA().ApplyHEEPIsoCut_Trk();
			treeVar_cut_eleB_stdTrkIso_ = diEle->eleB().ApplyHEEPIsoCut_Trk();
			treeVar_eleA_stdTrkIso_     = diEle->eleA().isolPtTrks();
			treeVar_eleB_stdTrkIso_     = diEle->eleB().isolPtTrks();

			treeVar_cut_eleA_modTrkIso_ = diEle->eleA_modTrkIsolCut();
			treeVar_cut_eleB_modTrkIso_ = diEle->eleB_modTrkIsolCut();
			treeVar_eleA_modTrkIso_     = diEle->eleA_modTrkIso();
			treeVar_eleB_modTrkIso_     = diEle->eleB_modTrkIso();

			treeVar_cut_eleA_stdEmH1Iso_ = diEle->eleA().ApplyHEEPIsoCut_EmHad1();
			treeVar_cut_eleB_stdEmH1Iso_ = diEle->eleB().ApplyHEEPIsoCut_EmHad1();
			treeVar_eleA_stdEmH1Iso_     = diEle->eleA().isolEmHadDepth1();
			treeVar_eleB_stdEmH1Iso_     = diEle->eleB().isolEmHadDepth1();

			treeVar_eleA_nTrksInnerVeto_ = diEle->eleA().numInnerIsoConeTrks();
			treeVar_eleB_nTrksInnerVeto_ = diEle->eleB().numInnerIsoConeTrks();

			treeVar_cut_eleA_scModEmH1Iso_ = diEle->eleA_modEmHad1IsoCut();
			treeVar_cut_eleB_scModEmH1Iso_ = diEle->eleB_modEmHad1IsoCut();
			treeVar_eleA_scModEmH1Iso_     = diEle->eleA_modEmHad1Iso();
			treeVar_eleB_scModEmH1Iso_     = diEle->eleB_modEmHad1Iso();

			treeVar_combThr_EmH1_ = 4.0 + 0.03*treeVar_eleA_p4_.Et() + 0.03*treeVar_eleB_p4_.Et();

			// Alternative isolation values, from the isoDeps modules ...
			treeVar_eleA_isoDep_stdTrk_  = diEle->eleA().isol_isoDep_stdTrk();
			treeVar_eleB_isoDep_stdTrk_  = diEle->eleB().isol_isoDep_stdTrk();
			treeVar_eleA_isoDep_stdEmH1_ = diEle->eleA().isol_isoDep_stdEmHadD1();
			treeVar_eleB_isoDep_stdEmH1_ = diEle->eleB().isol_isoDep_stdEmHadD1();
			treeVar_cut_eleA_isoDep_stdEmH1_ = diEle->eleA().isolCut_isoDep_stdEmHadD1();
			treeVar_cut_eleB_isoDep_stdEmH1_ = diEle->eleB().isolCut_isoDep_stdEmHadD1();
			treeVar_eleA_isoDep_inrVetoModTrk_ = diEle->eleA().isol_isoDep_inrVetoModTrk();
			treeVar_eleB_isoDep_inrVetoModTrk_ = diEle->eleB().isol_isoDep_inrVetoModTrk();
			treeVar_eleA_isoDep_inrVetoModEmH1_ = diEle->eleA().isol_isoDep_inrVetoModEmHadD1();
			treeVar_eleB_isoDep_inrVetoModEmH1_ = diEle->eleB().isol_isoDep_inrVetoModEmHadD1();
			treeVar_cut_eleA_isoDep_inrVetoModEmH1_ = diEle->eleA().isolCut_isoDep_inrVetoModEmHadD1();
			treeVar_cut_eleB_isoDep_inrVetoModEmH1_ = diEle->eleB().isolCut_isoDep_inrVetoModEmHadD1();

			// Alternative isolation values, from the BstdZee EDProducer ...
			treeVar_eleA_inrXSVetoModTrk_  = diEle->eleA().isol_inrVetoModTrk(tsw::Event::xSmallVeto);
			treeVar_eleA_inrSVetoModTrk_  = diEle->eleA().isol_inrVetoModTrk(tsw::Event::smallVeto);
			treeVar_eleA_inrMVetoModTrk_  = diEle->eleA().isol_inrVetoModTrk(tsw::Event::mediumVeto);
			treeVar_eleA_inrLVetoModTrk_  = diEle->eleA().isol_inrVetoModTrk(tsw::Event::largeVeto);
			treeVar_eleA_inrXLVetoModTrk_  = diEle->eleA().isol_inrVetoModTrk(tsw::Event::xLargeVeto);

			treeVar_eleB_inrXSVetoModTrk_  = diEle->eleB().isol_inrVetoModTrk(tsw::Event::xSmallVeto);
			treeVar_eleB_inrSVetoModTrk_  = diEle->eleB().isol_inrVetoModTrk(tsw::Event::smallVeto);
			treeVar_eleB_inrMVetoModTrk_  = diEle->eleB().isol_inrVetoModTrk(tsw::Event::mediumVeto);
			treeVar_eleB_inrLVetoModTrk_  = diEle->eleB().isol_inrVetoModTrk(tsw::Event::largeVeto);
			treeVar_eleB_inrXLVetoModTrk_  = diEle->eleB().isol_inrVetoModTrk(tsw::Event::xLargeVeto);

			treeVar_eleA_inrXSVetoModEmH1_ = diEle->eleA().isol_inrVetoModEmHadD1(tsw::Event::xSmallVeto);
			treeVar_eleA_inrSVetoModEmH1_  = diEle->eleA().isol_inrVetoModEmHadD1(tsw::Event::smallVeto);
			treeVar_eleA_inrMVetoModEmH1_  = diEle->eleA().isol_inrVetoModEmHadD1(tsw::Event::mediumVeto);
			treeVar_eleA_inrLVetoModEmH1_  = diEle->eleA().isol_inrVetoModEmHadD1(tsw::Event::largeVeto);
			treeVar_eleA_inrXLVetoModEmH1_ = diEle->eleA().isol_inrVetoModEmHadD1(tsw::Event::xLargeVeto);

			treeVar_eleB_inrXSVetoModEmH1_ = diEle->eleB().isol_inrVetoModEmHadD1(tsw::Event::xSmallVeto);
			treeVar_eleB_inrSVetoModEmH1_  = diEle->eleB().isol_inrVetoModEmHadD1(tsw::Event::smallVeto);
			treeVar_eleB_inrMVetoModEmH1_  = diEle->eleB().isol_inrVetoModEmHadD1(tsw::Event::mediumVeto);
			treeVar_eleB_inrLVetoModEmH1_  = diEle->eleB().isol_inrVetoModEmHadD1(tsw::Event::largeVeto);
			treeVar_eleB_inrXLVetoModEmH1_ = diEle->eleB().isol_inrVetoModEmHadD1(tsw::Event::xLargeVeto);

			treeVar_cut_eleA_inrMVetoModEmH1_ = diEle->eleA().isolCut_inrVetoModEmHadD1(tsw::Event::mediumVeto);
			treeVar_cut_eleB_inrMVetoModEmH1_ = diEle->eleB().isolCut_inrVetoModEmHadD1(tsw::Event::mediumVeto);

			treeVar_eleA_nGenHadronsDr04_ = diEle->eleA().isol_nGenHadronsDr04();
			treeVar_eleB_nGenHadronsDr04_ = diEle->eleB().isol_nGenHadronsDr04();
			treeVar_eleA_ptSumGenHadronsDr04_ = diEle->eleA().isol_ptSumGenHadronsDr04();
			treeVar_eleB_ptSumGenHadronsDr04_ = diEle->eleB().isol_ptSumGenHadronsDr04();

			// Rho PU correction to EmH1 isol'n values
			treeVar_eleA_EmH1RhoCorrn_ = diEle->eleA().isol_rhoCorrnEmH1(eventHelper);
			treeVar_eleB_EmH1RhoCorrn_ = diEle->eleB().isol_rhoCorrnEmH1(eventHelper);

			treeVar_eleA_heepIdModIsoCutCode_ = diEle->eleA().heepIdModIsoCutCode_v40_v00(eventHelper);
			treeVar_eleB_heepIdModIsoCutCode_ = diEle->eleB().heepIdModIsoCutCode_v40_v00(eventHelper);

			// And finally fill the tree ...
			effiCalcTree_->Fill();
		}

		void FillTree_NotReconstructed(const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2,
				const unsigned int runNum, const unsigned int lumiNum, const unsigned int evtNum, const tsw::EventHelper& eventHelper)
		{
			treeVar_weight_  = eventHelper.totWeight();
			treeVar_mc_numVtx_ = eventHelper.GetMCPU_nVtx();
			treeVar_runNum_  = runNum;
			treeVar_lumiSec_ = lumiNum;
			treeVar_evtNum_  = evtNum;

			treeVar_mcZ_ele1_p4_ = mcZboson_ele1;
			treeVar_mcZ_ele2_p4_ = mcZboson_ele2;

			bool mcZeles_zMass = fabs((mcZboson_ele1+mcZboson_ele2).M()-90.0)<30.0;
			bool mcZeles_pTacc = ( mcZboson_ele1.Pt()>40.0 && mcZboson_ele2.Pt()>40.0 );
			bool mcZele1_isEB = fabs(mcZboson_ele1.Eta())<1.4;
			bool mcZele2_isEB = fabs(mcZboson_ele2.Eta())<1.4;
			bool mcZele1_isEE = (fabs(mcZboson_ele1.Eta())>1.6 && fabs(mcZboson_ele1.Eta())<2.45);
			bool mcZele2_isEE = (fabs(mcZboson_ele2.Eta())>1.6 && fabs(mcZboson_ele2.Eta())<2.45);

			treeVar_ptAcc_ = mcZeles_pTacc;
			treeVar_ebebAcceptance_ = (mcZele1_isEB && mcZele2_isEB) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_ebeeAcceptance_ = ((mcZele1_isEB && mcZele2_isEE) || (mcZele1_isEE && mcZele2_isEB) ) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_eeeeAcceptance_ = (mcZele1_isEE && mcZele2_isEE) && (mcZeles_pTacc && mcZeles_zMass);
			treeVar_bothRecod_ = false;

			// Setting various variables to default values ...
			treeVar_ZpT_ = (mcZboson_ele1+mcZboson_ele2).Pt();
			treeVar_ZdEta_ = mcZboson_ele1.Eta() - mcZboson_ele2.Eta();
			treeVar_ZdPhi_ = mcZboson_ele1.DeltaPhi(mcZboson_ele2);
			treeVar_ZdR_   = mcZboson_ele1.DeltaR(mcZboson_ele2);
			if(mcZboson_ele1.Pt()>mcZboson_ele2.Pt()){
				treeVar_eleA_p4_ = mcZboson_ele1;
				treeVar_eleB_p4_ = mcZboson_ele2;
			}
			else{
				treeVar_eleA_p4_ = mcZboson_ele2;
				treeVar_eleB_p4_ = mcZboson_ele1;
			}

			treeVar_cut_both_fiducial_    = false;
			treeVar_cut_both_ecalDriven_  = false;
			treeVar_cut_both_dEta_        = false;
			treeVar_cut_both_dPhi_        = false;
			treeVar_cut_both_hOverE_      = false;
			treeVar_cut_both_showerShape_ = false;
			treeVar_cut_both_heepId_ = false;

			treeVar_cut_eleA_stdTrkIso_ = 9999.9;
			treeVar_cut_eleB_stdTrkIso_ = 9999.9;
			treeVar_eleA_stdTrkIso_     = 9999.9;
			treeVar_eleB_stdTrkIso_     = 9999.9;

			treeVar_cut_eleA_modTrkIso_ = 9999.9;
			treeVar_cut_eleB_modTrkIso_ = 9999.9;
			treeVar_eleA_modTrkIso_     = 9999.9;
			treeVar_eleB_modTrkIso_     = 9999.9;

			treeVar_cut_eleA_stdEmH1Iso_ = 9999.9;
			treeVar_cut_eleB_stdEmH1Iso_ = 9999.9;
			treeVar_eleA_stdEmH1Iso_     = 9999.9;
			treeVar_eleB_stdEmH1Iso_     = 9999.9;

			treeVar_eleA_nTrksInnerVeto_ = 9999.9;
			treeVar_eleB_nTrksInnerVeto_ = 9999.9;

			treeVar_cut_eleA_scModEmH1Iso_ = 9999.9;
			treeVar_cut_eleB_scModEmH1Iso_ = 9999.9;
			treeVar_eleA_scModEmH1Iso_     = 9999.9;
			treeVar_eleB_scModEmH1Iso_     = 9999.9;

			treeVar_combThr_EmH1_ = 9999.9;

			// Alternative isolation values, from the isoDeps modules ...
			treeVar_eleA_isoDep_stdTrk_  = 9999.9;
			treeVar_eleB_isoDep_stdTrk_  = 9999.9;
			treeVar_eleA_isoDep_stdEmH1_ = 9999.9;
			treeVar_eleB_isoDep_stdEmH1_ = 9999.9;
			treeVar_cut_eleA_isoDep_stdEmH1_ = 9999.9;
			treeVar_cut_eleB_isoDep_stdEmH1_ = 9999.9;
			treeVar_eleA_isoDep_inrVetoModTrk_ = 9999.9;
			treeVar_eleB_isoDep_inrVetoModTrk_ = 9999.9;
			treeVar_eleA_isoDep_inrVetoModEmH1_ = 9999.9;
			treeVar_eleB_isoDep_inrVetoModEmH1_ = 9999.9;
			treeVar_cut_eleA_isoDep_inrVetoModEmH1_ = 9999.9;
			treeVar_cut_eleB_isoDep_inrVetoModEmH1_ = 9999.9;

			// Alternative isolation values, from the BstdZee EDProducer ...
			treeVar_eleA_inrXSVetoModTrk_  = 9999.9;
			treeVar_eleA_inrSVetoModTrk_  = 9999.9;
			treeVar_eleA_inrMVetoModTrk_  = 9999.9;
			treeVar_eleA_inrLVetoModTrk_  = 9999.9;
			treeVar_eleA_inrXLVetoModTrk_  = 9999.9;

			treeVar_eleB_inrXSVetoModTrk_  = 9999.9;
			treeVar_eleB_inrSVetoModTrk_  = 9999.9;
			treeVar_eleB_inrMVetoModTrk_  = 9999.9;
			treeVar_eleB_inrLVetoModTrk_  = 9999.9;
			treeVar_eleB_inrXLVetoModTrk_  = 9999.9;

			treeVar_eleA_inrXSVetoModEmH1_ = 9999.9;
			treeVar_eleA_inrSVetoModEmH1_  = 9999.9;
			treeVar_eleA_inrMVetoModEmH1_  = 9999.9;
			treeVar_eleA_inrLVetoModEmH1_  = 9999.9;
			treeVar_eleA_inrXLVetoModEmH1_ = 9999.9;

			treeVar_eleB_inrXSVetoModEmH1_ = 9999.9;
			treeVar_eleB_inrSVetoModEmH1_  = 9999.9;
			treeVar_eleB_inrMVetoModEmH1_  = 9999.9;
			treeVar_eleB_inrLVetoModEmH1_  = 9999.9;
			treeVar_eleB_inrXLVetoModEmH1_ = 9999.9;

			treeVar_cut_eleA_inrMVetoModEmH1_ = 9999.9;
			treeVar_cut_eleB_inrMVetoModEmH1_ = 9999.9;

			treeVar_eleA_nGenHadronsDr04_ = 9999;
			treeVar_eleB_nGenHadronsDr04_ = 9999;
			treeVar_eleA_ptSumGenHadronsDr04_ = 9999.9;
			treeVar_eleB_ptSumGenHadronsDr04_ = 9999.9;

			// Rho PU correction to EmH1 isol'n values
			treeVar_eleA_EmH1RhoCorrn_ = 0.0;
			treeVar_eleB_EmH1RhoCorrn_ = 0.0;

			treeVar_eleA_heepIdModIsoCutCode_ = 0x01000-1;
			treeVar_eleB_heepIdModIsoCutCode_ = 0x01000-1;

			// And finally fill the tree ...
			effiCalcTree_->Fill();
		}

		void SaveToFile(TString outFileName){
			TFile f_tree(outFileName,"RECREATE");
			effiCalcTree_->Write();
			f_tree.Close();
		}
	};

}

//-----------------------------------------------------------------------------------------------------
//=========================================== Analyser class ==========================================

class BstdZeeFirstAnalyser{
	public:
		BstdZeeFirstAnalyser(int runMode, int numEvts, bool isMC, const TString& inFileName, const TString& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax);
		~BstdZeeFirstAnalyser();
		void DoAnalysis(const Double_t evtWeight);
		
	private:
		void AnalyseEvent(const Double_t );
		void SetupEleClassVectors();
		void SetupMuonCollection();
		void PrintOutBranchVariables();
		void FinishOffAnalysis(); //Method to be called after all events analysed - i.e. for normalisation of histograms, calculating errors on each bins, etc

		//Methods for application of HEEP cuts ...
		std::vector<bool> HEEPCutsWithoutIso(int eleType);

		// Method for obtaining MC-matched reconstructed di-electron
		tsw::HEEPDiEle* getMcMatchedDiEle(const TLorentzVector& , const TLorentzVector& , const std::vector<tsw::HEEPEle>& );

		//Methods for filling sets of histograms...
		void FillHistograms();
		void DoEMuAnalysis();

	public:
		// Method for retrieving the actual number of events run over ...
		unsigned int GetNumEvtsRunOver(){
			return numEvtsRunOver_;
		}
		// Methods for skipping events (for purposes of merging multiple MC samples of single process w/o overlap)
		bool SkipEvent(){
			if(skipFlg_highMCZpTEvts_ && mcZ_p4_ptr_->Pt()>skipThr_highMCZpTEvts_)
				return true;
			else if(skipFlg_lowMCZpTEvts_ && mcZ_p4_ptr_->Pt()<=skipThr_lowMCZpTEvts_)
				return true;
			else
				return false;
		}
		void SkipEvtIfMCZpTAbove(Float_t ptThreshold){
			skipFlg_highMCZpTEvts_ = true;
			skipThr_highMCZpTEvts_ = ptThreshold;
		}
		void SkipEvtIfMCZpTBelow(Float_t ptThreshold){
			skipFlg_lowMCZpTEvts_ = true;
			skipThr_lowMCZpTEvts_ = ptThreshold;
		}

	private:
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
		const int numEvtsRequested_;
		unsigned int numEvtsRunOver_;
		bool isMC_;
		const bool readInBstdEles_;

		TStopwatch timer_DoAnalysis_readIn_;
		TStopwatch timer_DoAnalysis_setup_;
		TStopwatch timer_DoAnalysis_FillHistograms_;
		TStopwatch timer_DoAnalysis_DoEMuMethod_;


		//Input and output files...
		TFile* inputFile_;
		TString inputFile_name_;
		TTree* inputFile_tree_;
		TString outputFile_name_;

		// Variables associated with skipping of events ...
		bool skipFlg_highMCZpTEvts_;
		Float_t skipThr_highMCZpTEvts_;
		bool skipFlg_lowMCZpTEvts_;
		Float_t skipThr_lowMCZpTEvts_;

		//------------//
		//Variables for storing contents of event information branches...
		tsw::Event dummyEvent_;
		tsw::Event* event_; //This has to be initialised to 0 in CTOR - otherwise weird runtime errors sometimes crop-up upon read out of branches ..
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
		std::vector<float>* normHEEPEles_ecalEnergyError_;
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
		std::vector<float>* normHEEPEles_closestCtfTrk_pt_;
		std::vector<float>* normHEEPEles_closestCtfTrk_eta_;
		std::vector<float>* normHEEPEles_closestCtfTrk_phi_;
		std::vector<float>* normHEEPEles_closestCtfTrk_innerPt_;
		std::vector<float>* normHEEPEles_closestCtfTrk_innerEta_;
		std::vector<float>* normHEEPEles_closestCtfTrk_innerPhi_;
		std::vector<float>* normHEEPEles_closestCtfTrk_outerPt_;
		std::vector<float>* normHEEPEles_closestCtfTrk_outerEta_;
		std::vector<float>* normHEEPEles_closestCtfTrk_outerPhi_;

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

		// SC information ...
		std::vector<float>* normHEEPEles_SCposn_eta_;
		std::vector<float>* normHEEPEles_SCposn_phi_;
		std::vector<float>* normHEEPEles_SC_rawEnergy_;
		std::vector< std::vector<float> >* normHEEPEles_SC_recHits_Et_;
		std::vector< std::vector<float> >* normHEEPEles_SC_recHits_eta_;
		std::vector< std::vector<float> >* normHEEPEles_SC_recHits_phi_;
		std::vector< std::vector<bool> >*  normHEEPEles_SC_recHits_isFromEB_;
		std::vector<float>* normHEEPEles_SC_totEnergyRecHits_;
		std::vector<unsigned int>* normHEEPEles_SC_totNumRecHits_;

		std::vector<float>* normHEEPEles_gsfTrk_eta_;
		std::vector<float>* normHEEPEles_gsfTrk_phi_;
		std::vector<float>* normHEEPEles_gsfTrk_vz_;
		std::vector< std::vector<float> >* normHEEPEles_innerIsoConeTrks_pt_;
		std::vector< std::vector<float> >* normHEEPEles_innerIsoConeTrks_eta_;
		std::vector< std::vector<float> >* normHEEPEles_innerIsoConeTrks_phi_;
		std::vector< std::vector<float> >* normHEEPEles_innerIsoConeTrks_vz_;

		std::vector<unsigned int>* normHEEPEles_numMissInnerHits_;

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
		tsw::MuonCollection normMuons_barrel_tight_;

		// Muon histograms ...
		tsw::MuonDistns normMuons_1stpT_Histos_;
		tsw::MuonDistns normMuons_tight_1stpT_Histos_;
		tsw::MuonDistns normMuons_barrel_tight_1stpT_Histos_;

		// emu method: emu object & output tree stuff
		tsw::HEEPEle normEles_EB_HEEPNoIso_1stpT_;
		bool normEles_EB_HEEPNoIso_1stpT_exists_;
		tsw::HEEPDiEle normDiEle_HEEPNoIso_;
		bool normDiEle_HEEPNoIso_exists_;
		tsw::EleMuObject eleMu_EB_HEEPNoIso_muB_tight_;

//		tsw::ABCDMethodTree frPreDiEleTree_;
		tsw::DiEleTree noIsoZCandDiEleTree_;
		tsw::DiEleTree modIsoZCandDiEleTree_;
		tsw::EffiCalcTree zCandEffiTree_;

		TTree* eleMuTree_;
		Double_t eleMuTreeVar_mass_;
		Double_t eleMuTreeVar_pT_;
		Double_t eleMuTreeVar_weight_;

		TTree* diEleTree_;
		Double_t diEleTreeVar_mass_;
		Double_t diEleTreeVar_pT_;
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

		//Di-electron pair histograms ...
		tsw::DiEleDistns normDiEle_AllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_AllHEEP_trgA_Histos_;
		tsw::DiEleDistns normDiEle_AllHEEP_MZ_Histos_;
		tsw::DiEleDistns normDiEle_AllHEEP_MZ_trgA_Histos_;
		tsw::DiEleDistnsByRgn normDiEle_AllHEEP_MZ_trgA_Rgnl_Histos_;
		tsw::DiEleDistns normDiEle_HEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_HEEPNoIso_trgA_Histos_;
		tsw::DiEleDistns normDiEle_HEEPNoIso_MZ_Histos_;
		tsw::DiEleDistns normDiEle_HEEPNoIso_MZ_trgA_Histos_;
		tsw::DiEleDistnsByRgn normDiEle_HEEPNoIso_MZ_trgA_Rgnl_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_trgA_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_trkIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_trgA_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_notEmHad1Iso_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_trgA_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIso_MZ_stdHEEPIso_Histos_;

		// ---------------------------------------------------- //
		// Event counters for checking effect of order of cuts in QCD estimation histos ...
		unsigned int numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_;
		unsigned int numEvts_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_;
		unsigned int numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_;
		//
		// Di-ele histograms for QCD estimation ...
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothHEEPNoIso_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothModAllHEEP_Histos_;
		tsw::DiEleDistns normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothModAllHEEP_Histos_;

		tsw::DiEleDistns normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_Histos_;
		tsw::DiEleDistns normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_modHEEPIsoOnHEEPNoIso_Histos_;
		// ---------------------------------------------------- //

		// E-mu object histograms ...
		tsw::EleMuDistns h_eleMu_EB_HEEPNoIso_muB_tight_MZ_;
		tsw::EleMuDistns h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_;
		tsw::EleMuDistns h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_HEEPIso_;

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

BstdZeeFirstAnalyser::BstdZeeFirstAnalyser(int runMode, int numEvts, bool isMC, const TString& inFileName, const TString& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax):
	numEvtsRequested_(numEvts),
	readInBstdEles_(false),
	skipFlg_highMCZpTEvts_(false),
	skipThr_highMCZpTEvts_(99999.9),
	skipFlg_lowMCZpTEvts_(false),
	skipThr_lowMCZpTEvts_(-99999.9),
	dummyEvent_(),
	event_(0),
	normMuons_1stpT_Histos_(             "h_normMuons_1stpT_",       "normal", "",           1, 1, 50, 1000.0, 50, 1000.0),
	normMuons_tight_1stpT_Histos_(       "h_normMuons_tight_1stpT_", "normal", "tight cuts", 2, 1, 50, 1000.0, 50, 1000.0),
	normMuons_barrel_tight_1stpT_Histos_("h_normMuons_barrel_tight_1stpT_", "normal", "barrel, tight cuts", 2, 1, 50, 1000.0, 50, 1000.0),
	//
	noIsoZCandDiEleTree_("noIsoZBosonTree"),
	modIsoZCandDiEleTree_("modIsoZBosonTree"),
	//
	normEles_reconValidationHistos_(    "h_normEles_",      "standard",  "",              1,  false),
	normEles_simpleCuts_reconValHistos_("h_normEles_sCuts_", "standard", ", simple cuts", 1,  false), //was 2
	normEles_HEEPCuts_reconValHistos_(  "h_normEles_HEEP_",  "standard", ", HEEP cuts",   1,  true),  //was 8
	bstdEles_reconValidationHistos_(    "h_bstdEles_",       "special",   "",             4,  false),//was 4
	bstdEles_simpleCuts_reconValHistos_("h_bstdEles_sCuts_", "special",  ", simple cuts", 4, false), //was 38
	bstdEles_HEEPCuts_reconValHistos_(  "h_bstdEles_HEEP_",  "special",  ", HEEP cuts",   4,  true),  //was 6
	//
	normDiEle_AllHEEP_Histos_(          "h_normDiEle_AllHEEP_",         "standard", "all HEEP cuts",              1, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_AllHEEP_trgA_Histos_(     "h_normDiEle_AllHEEP_trgA_",    "standard", "all HEEP cuts, trgA",        1, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_AllHEEP_MZ_Histos_(       "h_normDiEle_AllHEEP_MZ_",      "standard", "all HEEP cuts, Z mass",      1, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_AllHEEP_MZ_trgA_Histos_(  "h_normDiEle_AllHEEP_MZ_trgA_", "standard", "all HEEP cuts, Z mass, trgA",1, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_AllHEEP_MZ_trgA_Rgnl_Histos_("h_normDiEle_AllHEEP_MZ_trgA_", "standard", "all HEEP cuts, Z mass, trgA",1, 1, nbins_mass, nbins_pt, ptmax),//
	//
	normDiEle_HEEPNoIso_Histos_(        "h_normDiEle_HEEPNoIso_",         "standard", "HEEP cuts w/o isol",              2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_HEEPNoIso_trgA_Histos_(   "h_normDiEle_HEEPNoIso_trgA_",    "standard", "HEEP cuts w/o isol, trgA",        2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_HEEPNoIso_MZ_Histos_(     "h_normDiEle_HEEPNoIso_MZ_",      "standard", "HEEP cuts w/o isol, Z mass",      2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_HEEPNoIso_MZ_trgA_Histos_("h_normDiEle_HEEPNoIso_MZ_trgA_", "standard", "HEEP cuts w/o isol, Z mass, trgA",2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_HEEPNoIso_MZ_trgA_Rgnl_Histos_("h_normDiEle_HEEPNoIso_MZ_trgA_Rgnl_", "standard", "HEEP cuts w/o isol, Z mass, trgA",2, 1, nbins_mass, nbins_pt, ptmax),//
	//
	normDiEle_EB_HEEPNoIso_Histos_(        "h_normDiEle_EB_HEEPNoIso_",         "standard", "EB and HEEP cuts w/o isol",              2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_trgA_Histos_(   "h_normDiEle_EB_HEEPNoIso_trgA_",    "standard", "EB and HEEP cuts w/o isol, trgA",        2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_Histos_(     "h_normDiEle_EB_HEEPNoIso_MZ_",      "standard", "EB and HEEP cuts w/o isol, Z mass",      2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_trkIso_Histos_(     "h_normDiEle_EB_HEEPNoIso_MZ_trkIso_",      "standard", "EB and HEEP cuts w/o isol, Z mass, trkIso",      2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_Histos_(  "h_normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_",   "standard", "EB and HEEP cuts w/o isol, Z mass, EmHad1Iso",   2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_Histos_("h_normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_",   "standard", "EB and HEEP cuts w/o isol, Z mass, EmHad1Iso AND trkIso",   2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_trgA_Histos_("h_normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_trgA_",   "standard", "EB and HEEP cuts w/o isol, Z mass, EmHad1Iso AND trkIso, trgA passed",   2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_notEmHad1Iso_Histos_(  "h_normDiEle_EB_HEEPNoIso_MZ_notEmHad1Iso_",   "standard", "EB and HEEP cuts w/o isol, Z mass, notEmHad1Iso",   2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIso_MZ_trgA_Histos_("h_normDiEle_EB_HEEPNoIso_MZ_trgA_", "standard", "EB and HEEP cuts w/o isol, Z mass, trgA",2, 1, nbins_mass, nbins_pt, ptmax),//
	normDiEle_EB_HEEPNoIso_MZ_stdHEEPIso_Histos_("h_normDiEle_EB_HEEPNoIso_MZ_stdHEEPIso_", "standard", "EB and HEEP cuts w/o isol, Z mass, stdHEEPIso",2, 1, nbins_mass, nbins_pt, ptmax),//
	// QCD estimation member vars ...
	numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_(0),
	numEvts_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_(0),
	numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_(0),
	//
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_Histos_(                  "normDiEle_EB_fidECALDrFRPre_trgA_MZ_",                 "standard", "EB and fidECAL and FRPre, trgA, Z mass",                 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothHEEPNoIso_Histos_(    "normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothHEEPNoIso_",   "standard", "EB and fidECAL and FRPre, trgA, Z mass, bothHEEPNoIso",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothHEEPNoIso_Histos_( "normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothHEEPNoIso_",   "standard", "EB and fidECAL and FRPre, trgA, Z mass, <=1 pass HEEPNoIso",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothModAllHEEP_Histos_(   "normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothModAllHEEP_",  "standard", "EB and fidECAL and FRPre, trgA, Z mass, bothModAllHEEP", 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothModAllHEEP_Histos_("normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothModAllHEEP_",  "standard", "EB and fidECAL and FRPre, trgA, Z mass, <=1 pass ModAllHEEP", 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneHEEPNoIso_Histos_(    "normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneHEEPNoIso_",    "standard", "EB and fidECAL and FRPre, trgA, Z mass, oneHEEPNoIso",   2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneModAllHEEP_Histos_(   "normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneModAllHEEP_",   "standard", "EB and fidECAL and FRPre, trgA, Z mass, oneModAllHEEP",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroHEEPNoIso_Histos_(   "normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroHEEPNoIso_",   "standard", "EB and fidECAL and FRPre, trgA, Z mass, zeroHEEPNoIso",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroModAllHEEP_Histos_(  "normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroModAllHEEP_",  "standard", "EB and fidECAL and FRPre, trgA, Z mass, zeroModAllHEEP", 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneHEEPNoIso_Histos_( "normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneHEEPNoIso_", "standard", "EB and fidECAL and FRPre, trgA, Z mass, geqOneHEEPNoIso",2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneModAllHEEP_Histos_("normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneModAllHEEP_","standard", "EB and fidECAL and FRPre, trgA, Z mass, geqOneModAllHEEP",2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_Histos_(                  "normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_",                 "standard", "EB and fidECAL and FRPre, trgA, > Z mass",                 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothHEEPNoIso_Histos_(    "normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothHEEPNoIso_",   "standard", "EB and fidECAL and FRPre, trgA, > Z mass, bothHEEPNoIso",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothHEEPNoIso_Histos_( "normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothHEEPNoIso_",   "standard", "EB and fidECAL and FRPre, trgA, > Z mass, <=1 pass HEEPNoIso",  2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothModAllHEEP_Histos_(   "normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothModAllHEEP_",  "standard", "EB and fidECAL and FRPre, trgA, > Z mass, bothModAllHEEP", 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothModAllHEEP_Histos_("normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothModAllHEEP_",  "standard", "EB and fidECAL and FRPre, trgA, > Z mass, <=1 pass ModAllHEEP", 2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_Histos_(                      "normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_",                      "standard", "EB HEEPNoIso + inclGSF, trgA, Z mass",2, 1, nbins_mass, nbins_pt, ptmax),
	normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_modHEEPIsoOnHEEPNoIso_Histos_("normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_modHEEPIsoOnHEEPNoIso_","standard", "EB HEEPNoIso + inclGSF, trgA, Z mass, modHEEPIsoOnHEEPNoIso",2, 1, nbins_mass, nbins_pt, ptmax),
	//
	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_( "h_eleMu_EB_HEEPNoIso_muB_tight_MZ_", "standard", "EB+HEEPNoIso cuts, barrel+tight cuts, Z mass", 2, 1, nbins_mass, nbins_pt, ptmax ),
	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_( "h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_", "standard", "EB+HEEPNoIso cuts, barrel+tight cuts, Z mass, e-mu trigger", 2, 1, nbins_mass, nbins_pt, ptmax ),
	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_HEEPIso_( "h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_HEEPIso_", "standard", "EB+HEEPNoIso cuts, barrel+tight cuts, Z mass, e-mu trigger, HEEPIso", 2, 1, nbins_mass, nbins_pt, ptmax ),
	//
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

	timer_DoAnalysis_readIn_.Stop(); timer_DoAnalysis_readIn_.Reset();
	timer_DoAnalysis_setup_.Stop(); timer_DoAnalysis_setup_.Reset();
	timer_DoAnalysis_FillHistograms_.Stop(); timer_DoAnalysis_FillHistograms_.Reset();
	timer_DoAnalysis_DoEMuMethod_.Stop(); timer_DoAnalysis_DoEMuMethod_.Reset();

	//Initialise the member variables that are arguments of the CTOR...
	runMode_ = runMode;
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

	// ----- OPENING INPUT NTUPLE AND SETTING UP ASSOCIATED BRANCH LINKS ----- //
	//Opening the datafile(s)
	std::cout << " Opening the ntuple datafile..." << std::endl;
	inputFile_ = TFile::Open(inputFile_name_,"READ");
	/*if(!inputFile_ptr)
		std::cout << "ERROR in opening file." << std::endl;
	else
		std::cout << "   file opened successfully." << std::endl;*/

	//Changing into the correct directory within the input file ..
	inputFile_->cd("demo");
	//Getting a pointer to the ntuple tree in the file
	inputFile_tree_ = (TTree*)gDirectory->Get("EventDataTree");

	//Setup the branch links...
	SetupBranchLinks(inputFile_);

	inputFile_tree_->SetCacheSize(10000000);
	inputFile_tree_->AddBranchToCache("*");

	// ---------------------------------------------------- //
	// Determine the number of events being run over ...
	unsigned int inputTree_numEvts = inputFile_tree_->GetEntries();
	if( numEvtsRequested_<0 || (numEvtsRequested_>static_cast<int>(inputTree_numEvts)) )
		numEvtsRunOver_ = inputTree_numEvts;
	else
		numEvtsRunOver_ = numEvtsRequested_;
}

BstdZeeFirstAnalyser::~BstdZeeFirstAnalyser(){
	//Clear any heap variables here...
	DeleteReconValidationHistos();
}

void BstdZeeFirstAnalyser::DoAnalysis(const Double_t weightFromXsec)
{

	//Call AnalyseEvent method for each event...
	TStopwatch timer_AnalysingEvents; timer_AnalysingEvents.Start();
	for(unsigned int evtIdx=0; evtIdx<numEvtsRunOver_; evtIdx++){
		if(vFlg_>0){std::cout << std::endl << " Analysing event no. " << evtIdx << std::endl;}
		//Load in data for the evtIdx'th event...
		timer_DoAnalysis_readIn_.Start(false);
		Long64_t dataTreeEntry = inputFile_tree_->LoadTree(evtIdx);
		inputFile_tree_->GetEntry(dataTreeEntry);
		timer_DoAnalysis_readIn_.Stop();

		if( (vFlg_>-2) && (evtIdx%500000==0 || evtIdx==(numEvtsRunOver_-1)) ){std::cout << " *** Event no. " << evtIdx << " reached." << std::endl;}
		//TODO - Put setting up of TLorentzVector 4 momenta here ...

		// Skip to next pass through for loop IFF event should be skipped ...
		if(skipFlg_highMCZpTEvts_ && skipFlg_lowMCZpTEvts_){
			std::cout << "     ************************************************************************************ " << std::endl;
			std::cout << "     * ERROR: Both skipFlg_highMCZpTEvts_ and skipFlg_highMCZpTEvts_ are set to 'true'  * " << std::endl;
			std::cout << "     *     This may result in weird effects, and so am now breaking from event for loop * " << std::endl;
			std::cout << "     *     Please work out source of error and rectify.                                 * " << std::endl;
			std::cout << "     ************************************************************************************ " << std::endl;
			break;
		}
		if(SkipEvent())
			continue;

		//Analyse this data...
		AnalyseEvent(weightFromXsec);
	}
	timer_AnalysingEvents.Stop();
	timer_AnalysingEvents.Print();
	std::cout << "*** For timer_DoAnalysis_readIn_:" << std::endl;
	timer_DoAnalysis_readIn_.Print();
	std::cout << "*** For timer_DoAnalysis_setup_:" << std::endl;
	timer_DoAnalysis_setup_.Print();
	std::cout << "*** For timer_DoAnalysis_FillHistograms_:" << std::endl;
	timer_DoAnalysis_FillHistograms_.Print();
	std::cout << "*** For timer_DoAnalysis_DoEMuMethod_:" << std::endl;
	timer_DoAnalysis_DoEMuMethod_.Print();

	//Close the input file ...
	inputFile_->Close();

	FinishOffAnalysis(); //Output file is opened in here ...
}

//-------------------------------------------------------------------//
//---------------------- Private methods ----------------------------//
//-------------------------------------------------------------------//

void BstdZeeFirstAnalyser::AnalyseEvent(const Double_t weightFromXsec)
{
	eventHelper_ = tsw::EventHelper(event_);
	eventHelper_.setXsecWeight(weightFromXsec);
	//Setting up the TLorentzVector 4momenta...
	timer_DoAnalysis_setup_.Start(false);

	//mcEles_HighestEt_p4.SetPxPyPzE(mcEles_HighestEt_OLDp4_ptr->Px(),mcEles_HighestEt_OLDp4_ptr->Py(),mcEles_HighestEt_OLDp4_ptr->Pz(),mcEles_HighestEt_OLDp4_ptr->E());
	if(isMC_){
		mcEles_HighestEt_p4_ = ConvertToTLorentzVector(mcEles_HighestEt_OLDp4_ptr_);
		mcEles_2ndHighestEt_p4_ = ConvertToTLorentzVector(mcEles_2ndHighestEt_OLDp4_ptr_);
		mcZcandidate_p4_ = ConvertToTLorentzVector(mcZcandidate_OLDp4_ptr_);
		mcZ_daughterA_p4_ = ConvertToTLorentzVector(mcZ_daughterA_OLDp4_ptr_);
		mcZ_daughterB_p4_ = ConvertToTLorentzVector(mcZ_daughterB_OLDp4_ptr_);
	}

	normGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	/*ROOT::Math::XYZTVector tmpXYZTVector;
	for(unsigned int iEle = 0; iEle<normGsfEles_number_; iEle++){
		tmpXYZTVector = normGsfEles_OLDp4_->at(iEle);
		normGsfEles_p4_.push_back(  ConvertToTLorentzVector(&tmpXYZTVector)  );
	}
	bstdGsfEles_p4_.clear(); //std::vector<TLorentzVector>
	for(unsigned int iEle = 0; iEle<bstdGsfEles_number_; iEle++){
		tmpXYZTVector = bstdGsfEles_OLDp4_->at(iEle);
		bstdGsfEles_p4_.push_back(  ConvertToTLorentzVector(&tmpXYZTVector)  );
	}*/

	SetupEleClassVectors();

	if(isMC_)
		mcZboson_ = tsw::MCParticle(mcZ_pdgId_, mcZ_status_, -999, *mcZ_p4_ptr_);
	else{
		mcZboson_ = tsw::MCParticle();
		mcZ_daughterA_p4_.SetPxPyPzE(99999.9, 0.0, 0.0, 99999.9);
		mcZ_daughterB_p4_.SetPxPyPzE(99999.9, 0.0, 0.0, 99999.9);}


	if(vFlg_>0){PrintOutBranchVariables();}
	normEles_ = OrderHEEPElesByEt(normEles_);
	bstdEles_ = OrderHEEPElesByEt(bstdEles_);

	SetupMuonCollection();

	timer_DoAnalysis_setup_.Stop();
	///////////////////////////////////////
	// Now, can actually do some physics!
	timer_DoAnalysis_FillHistograms_.Start(false);
	FillHistograms();
	timer_DoAnalysis_FillHistograms_.Stop();

	timer_DoAnalysis_DoEMuMethod_.Start(false);
	DoEMuAnalysis();
	timer_DoAnalysis_DoEMuMethod_.Stop();
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
//		ithEleStruct.gsfEnergy_  = normHEEPEles_gsfEnergy_->at(iEle);
		ithEleStruct.caloEnergy_ = normHEEPEles_caloEnergy_->at(iEle);
//		ithEleStruct.ecalEnergyError_ = normHEEPEles_ecalEnergyError_->at(iEle);			/* TEMP v1f/g FIX */
		ithEleStruct.eta_        = normHEEPEles_eta_->at(iEle);
		ithEleStruct.scEta_      = normHEEPEles_scEta_->at(iEle);
//		ithEleStruct.detEta_     = normHEEPEles_detEta_->at(iEle);
//		ithEleStruct.detEtaAbs_  = normHEEPEles_detEtaAbs_->at(iEle);
		ithEleStruct.phi_        = normHEEPEles_phi_->at(iEle);
		ithEleStruct.scPhi_      = normHEEPEles_scPhi_->at(iEle);
//		ithEleStruct.detPhi_     = normHEEPEles_detPhi_->at(iEle);
//		ithEleStruct.zVtx_       = normHEEPEles_zVtx_->at(iEle);
		ithEleStruct.p4_         = normHEEPEles_p4ptr_->at(iEle);
		ithEleStruct.gsfP4_      = normHEEPEles_gsfP4ptr_->at(iEle);

		//classification (couldnt they have just named it 'type')
//		ithEleStruct.classification_ =  normHEEPEles_classification_->at(iEle);
		ithEleStruct.isEcalDriven_ = normHEEPEles_isEcalDriven_->at(iEle);
//		ithEleStruct.isTrackerDriven_ = normHEEPEles_isTrackerDriven_->at(iEle);
		ithEleStruct.isEB_ = normHEEPEles_isEB_->at(iEle);
		ithEleStruct.isEE_ = normHEEPEles_isEE_->at(iEle);

		//track methods
		ithEleStruct.charge_ = normHEEPEles_charge_->at(iEle);
//		ithEleStruct.trkCharge_ = normHEEPEles_trkCharge_->at(iEle);
//		ithEleStruct.pVtx_ = normHEEPEles_pVtx_->at(iEle);
//		ithEleStruct.pCalo_ = normHEEPEles_pCalo_->at(iEle);
		ithEleStruct.ptVtx_ = normHEEPEles_ptVtx_->at(iEle);
		ithEleStruct.ptCalo_ = normHEEPEles_ptCalo_->at(iEle);
		ithEleStruct.closestCtfTrk_pt_  = normHEEPEles_closestCtfTrk_pt_->at(iEle);
		ithEleStruct.closestCtfTrk_eta_ = normHEEPEles_closestCtfTrk_eta_->at(iEle);
		ithEleStruct.closestCtfTrk_phi_ = normHEEPEles_closestCtfTrk_phi_->at(iEle);
//		ithEleStruct.closestCtfTrk_innerPt_  = normHEEPEles_closestCtfTrk_innerPt_->at(iEle);
//		ithEleStruct.closestCtfTrk_innerEta_ = normHEEPEles_closestCtfTrk_innerEta_->at(iEle);
//		ithEleStruct.closestCtfTrk_innerPhi_ = normHEEPEles_closestCtfTrk_innerPhi_->at(iEle);
//		ithEleStruct.closestCtfTrk_outerPt_  = normHEEPEles_closestCtfTrk_outerPt_->at(iEle);
//		ithEleStruct.closestCtfTrk_outerEta_ = normHEEPEles_closestCtfTrk_outerEta_->at(iEle);
//		ithEleStruct.closestCtfTrk_outerPhi_ = normHEEPEles_closestCtfTrk_outerPhi_->at(iEle);

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		ithEleStruct.hOverE_ = normHEEPEles_hOverE_->at(iEle);
		ithEleStruct.dEtaIn_ = normHEEPEles_dEtaIn_->at(iEle);
		ithEleStruct.dPhiIn_ = normHEEPEles_dPhiIn_->at(iEle);
//		ithEleStruct.dPhiOut_ = normHEEPEles_dPhiOut_->at(iEle);
		ithEleStruct.epIn_ = normHEEPEles_epIn_->at(iEle);
		ithEleStruct.epOut_ = normHEEPEles_epOut_->at(iEle);
//		ithEleStruct.fbrem_ = normHEEPEles_fbrem_->at(iEle);
//		ithEleStruct.bremFrac_ = normHEEPEles_bremFrac_->at(iEle);
		ithEleStruct.invEOverInvP_ = normHEEPEles_invEOverInvP_->at(iEle);

		//shower shape variables
//		ithEleStruct.sigmaEtaEta_ = normHEEPEles_sigmaEtaEta_->at(iEle);
//		ithEleStruct.sigmaEtaEtaUnCorr_ = normHEEPEles_sigmaEtaEtaUnCorr_->at(iEle);
		ithEleStruct.sigmaIEtaIEta_ = normHEEPEles_sigmaIEtaIEta_->at(iEle);
//		ithEleStruct.e1x5_ = normHEEPEles_e1x5_->at(iEle);
//		ithEleStruct.e2x5Max_ = normHEEPEles_e2x5Max_->at(iEle);
//		ithEleStruct.e5x5_ = normHEEPEles_e5x5_->at(iEle);
		ithEleStruct.e1x5Over5x5_ = normHEEPEles_e1x5Over5x5_->at(iEle);
		ithEleStruct.e2x5MaxOver5x5_ = normHEEPEles_e2x5MaxOver5x5_->at(iEle);

		//isolation, we use cone of 0.3
		ithEleStruct.isolEm_ = normHEEPEles_isolEm_->at(iEle);
//		ithEleStruct.isolHad_ = normHEEPEles_isolHad_->at(iEle);
		ithEleStruct.isolHadDepth1_ = normHEEPEles_isolHadDepth1_->at(iEle);
		ithEleStruct.isolHadDepth2_ = normHEEPEles_isolHadDepth2_->at(iEle);
		ithEleStruct.isolPtTrks_ = normHEEPEles_isolPtTrks_->at(iEle);
		ithEleStruct.isolEmHadDepth1_ = normHEEPEles_isolEmHadDepth1_->at(iEle);

	  	ithEleStruct.isol_isoDep_stdTrk_   = eventHelper_.GetNormEle_IsoDep_stdTrkIso(iEle);
	  	ithEleStruct.isol_isoDep_stdEm_    = eventHelper_.GetNormEle_IsoDep_stdEcalIso(iEle);
	  	ithEleStruct.isol_isoDep_stdHadD1_ = eventHelper_.GetNormEle_IsoDep_stdHcalD1Iso(iEle);
	  	ithEleStruct.isol_isoDep_inrVetoModTrk_   = eventHelper_.GetNormEle_IsoDep_inrVetoModTrkIso(iEle);
	  	ithEleStruct.isol_isoDep_inrVetoModEm_    = eventHelper_.GetNormEle_IsoDep_inrVetoModEcalIso(iEle);
	  	ithEleStruct.isol_isoDep_inrVetoModHadD1_ = eventHelper_.GetNormEle_IsoDep_inrVetoModHcalD1Iso(iEle);

	  	ithEleStruct.isol_inrXSVetoModTrk_ = eventHelper_.GetNormEle_inrVetoModTrkIso(iEle, tsw::Event::xSmallVeto);
	  	ithEleStruct.isol_inrSVetoModTrk_ = eventHelper_.GetNormEle_inrVetoModTrkIso(iEle, tsw::Event::smallVeto);
	  	ithEleStruct.isol_inrMVetoModTrk_ = eventHelper_.GetNormEle_inrVetoModTrkIso(iEle, tsw::Event::mediumVeto);
	  	ithEleStruct.isol_inrLVetoModTrk_ = eventHelper_.GetNormEle_inrVetoModTrkIso(iEle, tsw::Event::largeVeto);
	  	ithEleStruct.isol_inrXLVetoModTrk_ = eventHelper_.GetNormEle_inrVetoModTrkIso(iEle, tsw::Event::xLargeVeto);

	  	ithEleStruct.isol_inrXSVetoModEm_ = eventHelper_.GetNormEle_inrVetoModEmIso(iEle, tsw::Event::xSmallVeto);
	  	ithEleStruct.isol_inrSVetoModEm_ = eventHelper_.GetNormEle_inrVetoModEmIso(iEle, tsw::Event::smallVeto);
	  	ithEleStruct.isol_inrMVetoModEm_ = eventHelper_.GetNormEle_inrVetoModEmIso(iEle, tsw::Event::mediumVeto);
	  	ithEleStruct.isol_inrLVetoModEm_ = eventHelper_.GetNormEle_inrVetoModEmIso(iEle, tsw::Event::largeVeto);
	  	ithEleStruct.isol_inrXLVetoModEm_ = eventHelper_.GetNormEle_inrVetoModEmIso(iEle, tsw::Event::xLargeVeto);

	  	ithEleStruct.isol_inrXSVetoModHadD1_ = eventHelper_.GetNormEle_inrVetoModHadD1Iso(iEle, tsw::Event::xSmallVeto);
	  	ithEleStruct.isol_inrSVetoModHadD1_ = eventHelper_.GetNormEle_inrVetoModHadD1Iso(iEle, tsw::Event::smallVeto);
	  	ithEleStruct.isol_inrMVetoModHadD1_ = eventHelper_.GetNormEle_inrVetoModHadD1Iso(iEle, tsw::Event::mediumVeto);
	  	ithEleStruct.isol_inrLVetoModHadD1_ = eventHelper_.GetNormEle_inrVetoModHadD1Iso(iEle, tsw::Event::largeVeto);
	  	ithEleStruct.isol_inrXLVetoModHadD1_ = eventHelper_.GetNormEle_inrVetoModHadD1Iso(iEle, tsw::Event::xLargeVeto);

		ithEleStruct.isol_nGenHadronsDr04_     = eventHelper_.GetNormEle_nGenHadronsDr04(iEle);
		ithEleStruct.isol_ptSumGenHadronsDr04_ = eventHelper_.GetNormEle_ptSumGenHadronsDr04(iEle);

//		ithEleStruct.SC_posn_eta_ = normHEEPEles_SCposn_eta_->at(iEle);							/* TEMP v1f/g FIX */
//		ithEleStruct.SC_posn_phi_ = normHEEPEles_SCposn_phi_->at(iEle);							/* TEMP v1f/g FIX */
//		ithEleStruct.SC_rawEnergy_ = normHEEPEles_SC_rawEnergy_->at(iEle);						/* TEMP v1f/g FIX */
		ithEleStruct.SC_recHits_Et_  = normHEEPEles_SC_recHits_Et_->at(iEle);
		ithEleStruct.SC_recHits_eta_ = normHEEPEles_SC_recHits_eta_->at(iEle);
		ithEleStruct.SC_recHits_phi_ = normHEEPEles_SC_recHits_phi_->at(iEle);
		ithEleStruct.SC_recHits_isFromEB_ = normHEEPEles_SC_recHits_isFromEB_->at(iEle);
//		ithEleStruct.SC_totEnergyRecHits_ = normHEEPEles_SC_totEnergyRecHits_->at(iEle);		/* TEMP v1f/g FIX */
//		ithEleStruct.SC_totNumRecHits_ = normHEEPEles_SC_totNumRecHits_->at(iEle);				/* TEMP v1f/g FIX */

		ithEleStruct.gsfTrk_eta_ = normHEEPEles_gsfTrk_eta_->at(iEle);
		ithEleStruct.gsfTrk_phi_ = normHEEPEles_gsfTrk_phi_->at(iEle);
		ithEleStruct.gsfTrk_vz_  = normHEEPEles_gsfTrk_vz_->at(iEle);

		ithEleStruct.innerIsoConeTrks_pt_  = normHEEPEles_innerIsoConeTrks_pt_->at(iEle);
		ithEleStruct.innerIsoConeTrks_eta_ = normHEEPEles_innerIsoConeTrks_eta_->at(iEle);
		ithEleStruct.innerIsoConeTrks_phi_ = normHEEPEles_innerIsoConeTrks_phi_->at(iEle);
		ithEleStruct.innerIsoConeTrks_vz_  = normHEEPEles_innerIsoConeTrks_vz_->at(iEle);

		ithEleStruct.numMissInnerHits_ = normHEEPEles_numMissInnerHits_->at(iEle);

		ithHEEPEle = tsw::HEEPEle::HEEPEle(ithEleStruct);
		normEles_.push_back(ithHEEPEle);
	}

	if(readInBstdEles_){
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
	} // End of if(readInBstdEles_)

	if(vFlg_>0){std::cout << "   done." << std::endl;}
}

void BstdZeeFirstAnalyser::SetupMuonCollection()
{

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

	normMuons_barrel_tight_ = normMuons_tight_.GetBarrelMuons();
	normMuons_barrel_tight_.OrderBypT();

	if(vFlg_>0)
		normMuons_tight_.Print();
}

void BstdZeeFirstAnalyser::PrintOutBranchVariables(){
	//Printing the event information to screen ...
	std::cout << "  ->Evt info (from event branch):"; event_->PrintBasicEventInformation();
	std::cout << "  ->Evt info: run no. " << evt_runNum_ << ", lumi sec. " << evt_lumiSec_<< ", event num. " << evt_evtNum_ << std::endl;

	event_->PrintPUVtxInfo();
	event_->PrintMCGenWeight();
	event_->PrintLHEZbosonInfo();

	//Printing the trigger decision to screen ...
	if(trg_PathA_decision_)
		std::cout << "  ->This event passed HLT trigger path A (i.e. " << *trg_PathA_name_ptr_ << ")" << std::endl;
	else
		std::cout << "  ->This event did *NOT* pass HLT trigger path A (i.e. " << *trg_PathA_name_ptr_ << ")" << std::endl;
	std::cout << "      Et of highest Et obj passing this path was: " << trg_PathA_highestTrigObjEt_ << std::endl;
	std::cout << "      (Last filter for this path was: " << *trg_PathA_nameOfLastFilter_ptr_ << ")" << std::endl;

	// Printing out e-mu trigger information to screen ...
	std::string eMuTrigger_name = eventHelper_.GetTrigInfo_eMuPath_name();
	if(eventHelper_.GetTrigInfo_eMuPath_decision())
		std::cout << "  ->This event passed the e-mu trigger (i.e. " << eMuTrigger_name << ")" << std::endl;
	else
		std::cout << "  ->This event did *NOT* pass the e-mu trigger (i.e. " << eMuTrigger_name << ")" << std::endl;


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
//	std::cout << "  ->There are " << normGsfEles_number_ << "=" << normGsfEles_charge_->size() << "=" << normGsfEles_OLDp4_->size() << "=" << normGsfEles_HEEP_Et_->size() << "=" << normGsfEles_ecalDriven_->size() << "?? 'normal' GSF electrons in this event." << std::endl;
	if(normGsfEles_number_!=normEles_.size()){
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
		std::cout << "         ***  ERROR: normGsfEles_number_ does NOT equal normEles_.size()" << std::endl;
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
	}
	else if(vFlg_>1)
		std::cout << "     normEles_.size()=" << normEles_.size() << std::endl;

	for(unsigned int iEle = 0; iEle < normGsfEles_number_; iEle++){
		std::cout << "    ->Normal GSF ele no. " << iEle << ":" << std::endl;
//		std::cout << " Charge = " << normGsfEles_charge_->at(iEle) << std::endl;
//
//		std::cout << "       p4: Px=" << normGsfEles_OLDp4_->at(iEle).Px() << "; Py=" << normGsfEles_OLDp4_->at(iEle).Py() << "; Pz=" << normGsfEles_OLDp4_->at(iEle).Pz() << std::endl;
//		std::cout << "           Px=" << normGsfEles_p4_.at(iEle).Px() << "; Py="  << normGsfEles_p4_.at(iEle).Py()  << "; Pz="  << normGsfEles_p4_.at(iEle).Pz() << std::endl;
//		std::cout << "           Pt=" << normGsfEles_p4_.at(iEle).Pt() << "; Eta=" << normGsfEles_p4_.at(iEle).Eta() << "; Phi=" << normGsfEles_p4_.at(iEle).Phi() << std::endl;
//
//		std::cout << "       Et=" << normGsfEles_Et_->at(iEle);
//		std::cout << "; HEEP_Et=" << normGsfEles_HEEP_Et_->at(iEle) << std::endl;
//		std::cout << "       Eta=" << normGsfEles_Eta_->at(iEle);
//		std::cout << "; scEta=" << normGsfEles_scEta_->at(iEle) << std::endl;
//		std::cout << "       ecalDriven=" << normGsfEles_ecalDriven_->at(iEle);
//		std::cout << "; ecalDrivenSeed=" << normGsfEles_ecalDrivenSeed_->at(iEle) << std::endl;
//
//		std::cout << "       HEEP_dEtaIn=" << normGsfEles_HEEP_dEtaIn_->at(iEle);
//		std::cout << "; HEEP_dPhiIn=" << normGsfEles_HEEP_dPhiIn_->at(iEle) << std::endl;
//		std::cout << "       HoverE=" << normGsfEles_HoverE_->at(iEle);
//		std::cout << "; sigmaIetaIeta=" << normGsfEles_sigmaIetaIeta_->at(iEle);
//		std::cout << "; scSigmaIetaIeta=" << normGsfEles_scSigmaIetaIeta_->at(iEle) << std::endl;
//
//		std::cout << "       dr03EmIsoEt=" << normGsfEles_dr03EmIsoEt_->at(iEle);
//		std::cout << "; dr03HadDepth1IsoEt=" << normGsfEles_dr03HadDepth1IsoEt_->at(iEle) << std::endl;
//		std::cout << "       dr03HadDepth2IsoEt=" << normGsfEles_dr03HadDepth2IsoEt_->at(iEle);
//		std::cout << "; dr03TkIsoPt=" << normGsfEles_dr03TkIsoPt_->at(iEle) << std::endl;
//
//		std::cout << "       e2x5Max=" << normGsfEles_e2x5Max_->at(iEle);
//		std::cout << "; e5x5=" << normGsfEles_e5x5_->at(iEle) << std::endl;

		std::cout << "     ->HEEP variables ..." << std::endl;
		normEles_.at(iEle).PrintOutVariables();
	}

	//Printing the information from the special boosted Z(ee) reconstruction GSF electron branches to screen ...
//	std::cout << "  ->There are " << bstdGsfEles_number_ << "=" << bstdGsfEles_charge_->size() << "=" << bstdGsfEles_OLDp4_->size() << "=" << bstdGsfEles_HEEP_Et_->size() << "=" << bstdGsfEles_ecalDriven_->size() << " 'special reconstruction' boosted electrons in this event." << std::endl;

	if(bstdGsfEles_number_!=bstdEles_.size()){
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
		std::cout << "         ***  ERROR: bstdGsfEles_number_ does NOT equal bstdEles_.size()" << std::endl;
		std::cout << "         ***--------------------------------------------------------***" << std::endl;
	}
	else if(vFlg_>1)
		std::cout << "     bstdEles_.size()=" << bstdEles_.size() << std::endl;

	for(unsigned int iEle = 0; iEle < bstdGsfEles_number_; iEle++){
		std::cout << "    ->Special reco'n GSF ele no. " << iEle << ":" << std::endl;
//		std::cout << " Charge = " << bstdGsfEles_charge_->at(iEle) << std::endl;
//
//		std::cout << "       p4: Px=" << bstdGsfEles_OLDp4_->at(iEle).Px() << "; Py=" << bstdGsfEles_OLDp4_->at(iEle).Py() << "; Pz=" << bstdGsfEles_OLDp4_->at(iEle).Pz() << std::endl;
//		std::cout << "           Px=" << bstdGsfEles_p4_.at(iEle).Px() << "; Py="  << bstdGsfEles_p4_.at(iEle).Py()  << "; Pz="  << bstdGsfEles_p4_.at(iEle).Pz() << std::endl;
//		std::cout << "           Pt=" << bstdGsfEles_p4_.at(iEle).Pt() << "; Eta=" << bstdGsfEles_p4_.at(iEle).Eta() << "; Phi=" << bstdGsfEles_p4_.at(iEle).Phi() << std::endl;
//
//		std::cout << "       Et=" << bstdGsfEles_Et_->at(iEle);
//		std::cout << "; HEEP_Et=" << bstdGsfEles_HEEP_Et_->at(iEle) << std::endl;
//		std::cout << "       Eta=" << bstdGsfEles_Eta_->at(iEle);
//		std::cout << "; scEta=" << bstdGsfEles_scEta_->at(iEle) << std::endl;
//		std::cout << "       ecalDriven=" << bstdGsfEles_ecalDriven_->at(iEle);
//		std::cout << "; ecalDrivenSeed=" << bstdGsfEles_ecalDrivenSeed_->at(iEle) << std::endl;
//
//		std::cout << "       HEEP_dEtaIn=" << bstdGsfEles_HEEP_dEtaIn_->at(iEle);
//		std::cout << "; HEEP_dPhiIn=" << bstdGsfEles_HEEP_dPhiIn_->at(iEle) << std::endl;
//		std::cout << "       HoverE=" << bstdGsfEles_HoverE_->at(iEle);
//		std::cout << "; sigmaIetaIeta=" << bstdGsfEles_sigmaIetaIeta_->at(iEle);
//		std::cout << "; scSigmaIetaIeta=" << bstdGsfEles_scSigmaIetaIeta_->at(iEle) << std::endl;
//
//		std::cout << "       dr03EmIsoEt=" << bstdGsfEles_dr03EmIsoEt_->at(iEle);
//		std::cout << "; dr03HadDepth1IsoEt=" << bstdGsfEles_dr03HadDepth1IsoEt_->at(iEle) << std::endl;
//		std::cout << "       dr03HadDepth2IsoEt=" << bstdGsfEles_dr03HadDepth2IsoEt_->at(iEle);
//		std::cout << "; dr03TkIsoPt=" << bstdGsfEles_dr03TkIsoPt_->at(iEle) << std::endl;
//
//		std::cout << "       e2x5Max=" << bstdGsfEles_e2x5Max_->at(iEle);
//		std::cout << "; e5x5=" << bstdGsfEles_e5x5_->at(iEle) << std::endl;

		std::cout << "     ->HEEP variables ..." << std::endl;
		bstdEles_.at(iEle).PrintOutVariables();
	}
}

void BstdZeeFirstAnalyser::FinishOffAnalysis(){

	//Calculate errors on histogram bins ...



	//Open the output file ...
	TFile f_histos(outputFile_name_+"_histos.root","RECREATE");
	f_histos.Write();

	//Save all of the histograms from analysis ...
	SaveReconValidationHistos(&f_histos);
	//Close the output file ...
	f_histos.Close();

	// Open up eMu method files, save the appropriate tree, and then
	TFile f_eleMuTree("eleMuTree/" + outputFile_name_ + "_eMuTree.root","RECREATE");
	eleMuTree_->Write();
	f_eleMuTree.Close();
	delete eleMuTree_;

	TFile f_diEleTree("diEleTree/" + outputFile_name_ + "_diEleTree.root","RECREATE");
	diEleTree_->Write();
	f_diEleTree.Close();
	delete diEleTree_;

	// Save the ABCD QCD estimation tree ...
//	frPreDiEleTree_.SaveToFile("/opt/ppd/newscratch/williams/Datafiles/abcdDiEleTrees/" + outputFile_name_ + "_abcdTree.root");

	// Save the Z candidate di-ele tree ...
	noIsoZCandDiEleTree_.SaveToFile(outputFile_name_ + "_noIsoZCandTree.root");
	modIsoZCandDiEleTree_.SaveToFile(outputFile_name_ + "_modIsoZCandTree.root");
	zCandEffiTree_.SaveToFile(outputFile_name_ + "_zEffiTree.root");

	// Output information to screen about diff ordering of cuts in QCD estimation
	//numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_
	std::cout << std::endl;
	std::cout << " ***** ---------------- QCD nos ----------------- *****" << std::endl;
	std::cout << " * Number of events with normDiEle_EB_HEEPNoIso_MZ_trgA:" << std::endl;
	std::cout << " *             " << numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_ << std::endl;
	//numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_
	std::cout << " * Number of these events in which normDiEle_EB_HEEPNoIso_MZ_trgA != normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso:" << std::endl;
	std::cout << " *             " << numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_ << std::endl;
	//numEvts_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_
	std::cout << " * Number of events with normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso:" << std::endl;
	std::cout << " *             " << numEvts_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_ << std::endl;
	std::cout << " ***** ------------------------------------------ *****" << std::endl;

	//numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_
}


//Methods for application of HEEP cuts ...
std::vector<bool> BstdZeeFirstAnalyser::HEEPCutsWithoutIso(int eleType){
	//Currently a placeholder ...
	std::vector<bool> tmpAnswer;
	tmpAnswer.clear();
	return tmpAnswer;
}


// Method for obtaining MC-matched reconstructed di-electron
tsw::HEEPDiEle* BstdZeeFirstAnalyser::getMcMatchedDiEle(const TLorentzVector& p4_mcEle1, const TLorentzVector& p4_mcEle2, const std::vector<tsw::HEEPEle>& vecOfRecoEles)
{
	if(vecOfRecoEles.size()>1){
		// 1. Calculate dR of each reco'd electron with respect to each of MC eles
		std::vector<double> dR_vsMcEle1, dR_vsMcEle2;
		for(std::vector<tsw::HEEPEle>::const_iterator recoEleIt = vecOfRecoEles.begin();
				recoEleIt != vecOfRecoEles.end(); recoEleIt++){
			dR_vsMcEle1.push_back( p4_mcEle1.DeltaR(recoEleIt->p4()) );
			dR_vsMcEle2.push_back( p4_mcEle2.DeltaR(recoEleIt->p4()) );
		}

		// 2. Determine which reco'd electron is closest to each of MC eles
		unsigned int mcEle1_closestRecoEle_idx = 9999;
		unsigned int mcEle2_closestRecoEle_idx = 9999;
		double mcEle1_closestRecoEle_dR = 9999.9;
		double mcEle2_closestRecoEle_dR = 9999.9;

		for(unsigned int recoEleIdx = 0; recoEleIdx<vecOfRecoEles.size(); recoEleIdx++){
			if( !(vecOfRecoEles.at(recoEleIdx).isEcalDriven()) )
				continue;
			// For MC ele 1 ...
			double dRvsMcEle1 = p4_mcEle1.DeltaR( vecOfRecoEles.at(recoEleIdx).p4() );
			if(dRvsMcEle1<mcEle1_closestRecoEle_dR){
				mcEle1_closestRecoEle_idx = recoEleIdx;
				mcEle1_closestRecoEle_dR = dRvsMcEle1;
			}
			// For MC ele 2 ...
			double dRvsMcEle2 = p4_mcEle2.DeltaR( vecOfRecoEles.at(recoEleIdx).p4() );
			if(dRvsMcEle2<mcEle2_closestRecoEle_dR){
				mcEle2_closestRecoEle_idx = recoEleIdx;
				mcEle2_closestRecoEle_dR = dRvsMcEle2;
			}
		}

		// 3. If either of MC electrons has dR>0.05 wrt closest reco ele, or if closest reco ele is same for two MC eles, then return null pointer
		if( (mcEle1_closestRecoEle_dR>0.05 || mcEle2_closestRecoEle_dR>0.05) || mcEle1_closestRecoEle_idx==mcEle2_closestRecoEle_idx )
			return 0;
		else
			return ( new tsw::HEEPDiEle(vecOfRecoEles.at(mcEle1_closestRecoEle_idx), vecOfRecoEles.at(mcEle2_closestRecoEle_idx) ) );
	}
	else
		return 0;
}

//--------------------------------------------------------//
//---- Methods for filling sets of histograms ...
//--------------------------------------------------------//


void BstdZeeFirstAnalyser::FillHistograms()
{
	const Double_t evtWeight = eventHelper_.totWeight();

	normDiEle_HEEPNoIso_exists_ = false;
	normEles_EB_HEEPNoIso_1stpT_exists_ = false;

	//Currently a placeholder ...
	unsigned int idx_highestEtEle = 0;
	float highestEleEt = -999.9;
	std::vector<bool> normEles_sCutsFlags;           normEles_sCutsFlags.clear();
	std::vector<bool> normEles_HEEPCutsFlags;        normEles_HEEPCutsFlags.clear();
	std::vector<bool> normEles_HEEPCutsNoIsoFlags;   normEles_HEEPCutsNoIsoFlags.clear();
	std::vector<bool> normEles_EB_HEEPCutsNoIsoFlags;normEles_EB_HEEPCutsNoIsoFlags.clear();
	std::vector<bool> normEles_EB_heepIdModIsoFlags;
//	std::vector<bool> normEles_EB_isFidAndEcalDrivenFlags; normEles_EB_isFidAndEcalDrivenFlags.clear();
	std::vector<bool> normEles_EB_isFidEcalDrAndFRPreFlags; normEles_EB_isFidEcalDrAndFRPreFlags.clear();

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
		normEles_reconValidationHistos_.FillHistos( normEles_.at(iEle), evtWeight);
		if ( normEles_.at(iEle).ApplySimpleCuts() ){
			normEles_simpleCuts_reconValHistos_.FillHistos( normEles_.at(iEle), evtWeight);}
		//if ( normEles_.at(iEle).ApplyAllHEEPCuts() ){
		//	normEles_HEEPCuts_reconValHistos_.FillHistos( normEles_.at(iEle) );}

		normEles_sCutsFlags.push_back(   normEles_.at(iEle).ApplySimpleCuts());
		normEles_HEEPCutsFlags.push_back(normEles_.at(iEle).ApplyAllHEEPCuts()  );
		normEles_HEEPCutsNoIsoFlags.push_back(normEles_.at(iEle).ApplyHEEPCutsNoIso()  );
		normEles_EB_HEEPCutsNoIsoFlags.push_back( normEles_.at(iEle).ApplyHEEPCutsNoIso() && normEles_.at(iEle).isEB() );

		normEles_EB_heepIdModIsoFlags.push_back( normEles_.at(iEle).heepIdModIsoCut(eventHelper_) && normEles_.at(iEle).isEB() );

//		normEles_EB_isFidAndEcalDrivenFlags.push_back( normEles_.at(iEle).isFiducialAndEcalDriven() && normEles_.at(iEle).isHEEPEB() );
		normEles_EB_isFidEcalDrAndFRPreFlags.push_back( normEles_.at(iEle).isFidEcalDrAndPassesFRPre() && normEles_.at(iEle).isHEEPEB() );

		//if(normEles_.at(iEle).et()>highestEleEt){
		//	highestEleEt = normEles_.at(iEle).et();
		//	idx_highestEtEle = iEle;
		//}
	}
	normEles_HEEPCuts_reconValHistos_.FillHistos(normEles_, normEles_HEEPCutsFlags, evtWeight);

	idx_highestEtEle = 0;
	highestEleEt = -999.9;

	//	if(readInBstdEles_){
	//		bstdEles_reconValidationHistos_.FillNumberElesHist(bstdGsfEles_number_ , evtWeight);
	//		for(unsigned int iEle=0; iEle<bstdGsfEles_number_; iEle++){
	//			bstdEles_reconValidationHistos_.FillHistos( bstdEles_.at(iEle), evtWeight);
	//			if ( bstdEles_.at(iEle).ApplySimpleCuts() ){
	//				bstdEles_simpleCuts_reconValHistos_.FillHistos( bstdEles_.at(iEle), evtWeight);}
	//			//if ( bstdEles_.at(iEle).ApplyAllHEEPCuts() ){
	//			//	bstdEles_HEEPCuts_reconValHistos_.FillHistos( bstdEles_.at(iEle) );}
	//
	//			bstdEles_sCutsFlags.push_back(         bstdEles_.at(iEle).ApplySimpleCuts()    );
	//			bstdEles_HEEPCutsFlags.push_back(      bstdEles_.at(iEle).ApplyAllHEEPCuts()      );
	//			bstdEles_HEEPCutsNoIsoFlags.push_back( bstdEles_.at(iEle).ApplyHEEPCutsNoIso() );
	//
	//			//if(bstdEles_.at(iEle).et()>highestEleEt){
	//			//	highestEleEt = bstdEles_.at(iEle).et();
	//			//	idx_highestEtEle = iEle;
	//			//}
	//		}
	//		//if(bstdGsfEles_number_>0)
	//		//	bstdEles_reconValidationHistos_.FillHighestEtEleHistos( bstdEles_.at(idx_highestEtEle) );
	//		bstdEles_HEEPCuts_reconValHistos_.FillHistos(bstdEles_, bstdEles_HEEPCutsFlags, evtWeight);
	//	} // End of if(readInBstdEles_)


	if(tsw::NumPassingCuts(normEles_EB_HEEPCutsNoIsoFlags)>0){
		std::vector<unsigned int> idxs_elesPassedCuts = IndicesOfElesPassingCuts(normEles_, normEles_EB_HEEPCutsNoIsoFlags, 0);
		normEles_EB_HEEPNoIso_1stpT_ = normEles_.at( idxs_elesPassedCuts.at(0) );
		normEles_EB_HEEPNoIso_1stpT_exists_ = true;
	}

	/////////////////////////////////////////////////
	// Form the di-electron pairs ...
	tsw::HEEPDiEle normDiEle_EB_fidECALDrFRPre;
	tsw::HEEPDiEle normDiEle_AllHEEP;
	tsw::HEEPDiEle normDiEle_HEEPNoIso;
	tsw::HEEPDiEle normDiEle_EB_HEEPNoIso;

	tsw::HEEPDiEle bstdDiEle_AllHEEP;
	tsw::HEEPDiEle bstdDiEle_HEEPNoIso;

	// ------------------------------
	// HEEP (ALL cuts) di-electrons ...
	if(tsw::NumPassingCuts(normEles_HEEPCutsFlags)>1){
		normDiEle_AllHEEP = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_HEEPCutsFlags );
		normDiEle_AllHEEP_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);
		if( normDiEle_AllHEEP.eleA().isEB() && normDiEle_AllHEEP.eleB().isEB() )
			normDiEle_AllHEEP_EBEB_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);

		if(trg_PathA_decision_)
			normDiEle_AllHEEP_trgA_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);

		// Check if the di-electron mass is in the Z mass range ...
		if(normDiEle_AllHEEP.isInZMassRange()){
			normDiEle_AllHEEP_MZ_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);
			if(trg_PathA_decision_){
				normDiEle_AllHEEP_MZ_trgA_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);
				normDiEle_AllHEEP_MZ_trgA_Rgnl_Histos_.FillHistos(normDiEle_AllHEEP, evtWeight);}
		}
	}

	// ------------------------
	// HEEPNoIso di-electrons ...
	if(tsw::NumPassingCuts(normEles_HEEPCutsNoIsoFlags)>1){
		// Set some of the emu method diele tree branch variablles here ...
		normDiEle_HEEPNoIso = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_HEEPCutsNoIsoFlags );

		normDiEle_HEEPNoIso_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);
		if(trg_PathA_decision_)
			normDiEle_HEEPNoIso_trgA_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);

		// Check if the di-electron mass is in the Z mass range ...
		if(normDiEle_HEEPNoIso.isInZMassRange()){
			normDiEle_HEEPNoIso_MZ_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);
			if(trg_PathA_decision_){
				normDiEle_HEEPNoIso_MZ_trgA_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);
				normDiEle_HEEPNoIso_MZ_trgA_Rgnl_Histos_.FillHistos(normDiEle_HEEPNoIso, evtWeight);}

			mcZboson_stdReconNoIsoEvts_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);

			if( fabs(mcZ_daughterA_p4_.Eta())<=1.5 && fabs(mcZ_daughterB_p4_.Eta())<=1.5 )
				mcZboson_stdReconNoIsoEvts_EBEB_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
			else if( fabs(mcZ_daughterA_p4_.Eta())>=1.56 && fabs(mcZ_daughterB_p4_.Eta())>=1.56 )
				mcZboson_stdReconNoIsoEvts_EEEE_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
		}
	}

	// ------------------------
	// EB HEEPNoIso di-electrons...
	if(tsw::NumPassingCuts(normEles_EB_HEEPCutsNoIsoFlags)>1){
		normDiEle_EB_HEEPNoIso = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_EB_HEEPCutsNoIsoFlags);
		//		//-----------------------------------------------
		//		// Code for testing ApplyDiEleTrkIsolCut method ...
		//		std::cout << std::endl << "   ****------------------------------------------------------------------****" << std::endl;
		//		std::cout << "   * EB HEEPNoIso di-electron formed in this event!!" << std::endl;
		//		normDiEle_EB_HEEPNoIso.PrintOutInfo();
		//		if(normDiEle_EB_HEEPNoIso.ApplyDiEleTrkIsolCut())
		//			std::cout << " ->* This di-electron pair has passed cut in ApplyDiEleTrkIsolCut() method." << std::endl;
		//		else
		//			std::cout << " ->* This di-electron pair has NOT passed cut in ApplyDiEleTrkIsolCut() method." << std::endl;
		////		}
		//		normDiEle_EB_HEEPNoIso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight);
		//		if(trg_PathA_decision_)
		//			normDiEle_EB_HEEPNoIso_trgA_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight);
		//		//-----------------------------------------------

		// Check if the di-electron mass is in the Z mass range ...
//		zCandEffiTree_.FillTree(&normDiEle_EB_HEEPNoIso, mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_, evtWeight);
		if(normDiEle_EB_HEEPNoIso.isInZMassRange()){
			normDiEle_EB_HEEPNoIso_MZ_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);

			noIsoZCandDiEleTree_.FillTree(normDiEle_EB_HEEPNoIso, eventHelper_, trg_PathA_decision_);
			// Apply modified track isolation cut ...
			if(normDiEle_EB_HEEPNoIso.ApplyDiEleTrkIsolCut())
				normDiEle_EB_HEEPNoIso_MZ_trkIso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);
			// Apply modified emHad1 isolation cut ...
			if(normDiEle_EB_HEEPNoIso.ApplyDiEleEmHad1IsolCut()){
				normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);
				// Apply modified track isolation cut ...
				if(normDiEle_EB_HEEPNoIso.ApplyDiEleTrkIsolCut()){
					normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);


					// Check that signal trigger has fired ...
					if(trg_PathA_decision_){
						normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_trgA_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);
						// Set-up tree variables for e-mu method di-ele tree
						normDiEle_HEEPNoIso_ = normDiEle_EB_HEEPNoIso;
						normDiEle_HEEPNoIso_exists_ =  true;
					}
				}// if-End: DiEleTrkIsolCut
			}// if-End: DiEleEmHad1IsoCut
			else
				normDiEle_EB_HEEPNoIso_MZ_notEmHad1Iso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);

			if(normDiEle_EB_HEEPNoIso.eleA().ApplyIsoVarHEEPCuts() && normDiEle_EB_HEEPNoIso.eleB().ApplyIsoVarHEEPCuts())
				normDiEle_EB_HEEPNoIso_MZ_stdHEEPIso_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight, true);

			if(trg_PathA_decision_){
				numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_++;
				normDiEle_EB_HEEPNoIso_MZ_trgA_Histos_.FillHistos(normDiEle_EB_HEEPNoIso, evtWeight);

				// Check if this di-ele is the same as would have been produced when applying HEEPNoIso cut to di-ele formed from two highest Et fiducial ECAL-driven eles AND then applying HEEPNoIso cuts
				unsigned int idx_eleA_HEEPNoIso = IndicesOfElesPassingCuts(normEles_, normEles_EB_HEEPCutsNoIsoFlags, 0).at(0); // Grab indices of two highest Et EB HEEPNoIso eles
				unsigned int idx_eleB_HEEPNoIso = IndicesOfElesPassingCuts(normEles_, normEles_EB_HEEPCutsNoIsoFlags, 0).at(1);
				if(tsw::NumPassingCuts(normEles_EB_isFidEcalDrAndFRPreFlags)>1){
					unsigned int idx_eleA_fidEcalDrAndFRPre = IndicesOfElesPassingCuts(normEles_, normEles_EB_isFidEcalDrAndFRPreFlags, 0).at(0);// Grab indices of two highest Et EB fidECALDrFRPre eles
					unsigned int idx_eleB_fidEcalDrAndFRPre = IndicesOfElesPassingCuts(normEles_, normEles_EB_isFidEcalDrAndFRPreFlags, 0).at(1);

					// The two di-eles previously discussed will not be the same if the eles forming them are different ...
					if( (idx_eleA_HEEPNoIso!=idx_eleA_fidEcalDrAndFRPre) || (idx_eleB_HEEPNoIso!=idx_eleB_fidEcalDrAndFRPre) )
						numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_++;
				}
				else
					numEvts_normDiEle_EB_HEEPNoIso_MZ_trgA_DiffFrom_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_++;

			}
		}
	}

	if( tsw::HEEPDiEle* mcMatchedDiEle = getMcMatchedDiEle(mcZ_daughterA_p4_, mcZ_daughterB_p4_, normEles_) )
	{
		zCandEffiTree_.FillTree(mcMatchedDiEle, mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_);
		delete mcMatchedDiEle;
	}
	else
		zCandEffiTree_.FillTree_NotReconstructed(mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_);

	// ------------------------
	// EB HEEPNoIso di-electrons...
	if(tsw::NumPassingCuts(normEles_EB_heepIdModIsoFlags)>1){
		tsw::HEEPDiEle heepIdModIsoEbDiEle( normEles_, normEles_EB_heepIdModIsoFlags);

		if( heepIdModIsoEbDiEle.isInZMassRange() )
			modIsoZCandDiEleTree_.FillTree(heepIdModIsoEbDiEle, eventHelper_, trg_PathA_decision_);

	}


	//-----------------------------------
	// ABCD method selection code:
	// EB fid, ECAL-driven, fake rate (FR) pre-selection di-electrons
	if( tsw::NumPassingCuts(normEles_EB_isFidEcalDrAndFRPreFlags)>1 ){
		// Form the di-electron
		normDiEle_EB_fidECALDrFRPre = tsw::HEEPDiEle::HEEPDiEle( normEles_, normEles_EB_isFidEcalDrAndFRPreFlags);
		// Require signal trigger to have fired
		if(trg_PathA_decision_){
			std::string trigPathA_name = *trg_PathA_name_ptr_;
//			frPreDiEleTree_.FillTree(&normDiEle_EB_fidECALDrFRPre, trigPathA_name, trg_PathA_decision_, evtWeight);
			//Apply Z mass window cut
			if(normDiEle_EB_fidECALDrFRPre.isInZMassRange()){
				normDiEle_EB_fidECALDrFRPre_trgA_MZ_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				// Setup bools for whether eleA and eleB pass HEEPNoIso ...
				bool eleA_HEEPNoIso = normDiEle_EB_fidECALDrFRPre.eleA().ApplyHEEPCutsNoIso();
				bool eleB_HEEPNoIso = normDiEle_EB_fidECALDrFRPre.eleB().ApplyHEEPCutsNoIso();
				// Setup bools for whether eleA and eleB pass HEEPNoIso+modIso ...
				bool eleA_modAllHEEP = ( eleA_HEEPNoIso && ( normDiEle_EB_fidECALDrFRPre.eleA_modTrkIsolCut() && normDiEle_EB_fidECALDrFRPre.eleA_modEmHad1IsoCut() ) );
				bool eleB_modAllHEEP = ( eleB_HEEPNoIso && ( normDiEle_EB_fidECALDrFRPre.eleB_modTrkIsolCut() && normDiEle_EB_fidECALDrFRPre.eleB_modEmHad1IsoCut() ) );

				// Region A: Both electrons pass HEEP cuts
				if(eleA_HEEPNoIso && eleB_HEEPNoIso){
					numEvts_normDiEle_EB_fidEcal_trgA_MZ_HEEPNoIso_++;
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				}
				else // (Region B in bothModHEEP vs invMass)
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				if(eleA_modAllHEEP && eleB_modAllHEEP)
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				else // (Region B in bothModHEEP vs invMass)
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				// Region B/C: Exactly one of electrons passes HEEP cuts
				if( (eleA_HEEPNoIso && (!eleB_HEEPNoIso)) || (eleB_HEEPNoIso && (!eleA_HEEPNoIso)) )
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				if( (eleA_modAllHEEP && (!eleB_modAllHEEP)) || (eleB_modAllHEEP && (!eleA_modAllHEEP)) )
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				// Region D: Neither electron passes HEEP cuts
				if( (!eleA_HEEPNoIso) && (!eleB_HEEPNoIso))
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				if( (!eleA_modAllHEEP) && (!eleB_modAllHEEP) )
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				// Other QCD plot: >=1 electron passes HEEP cuts
				if( eleA_HEEPNoIso || eleB_HEEPNoIso )
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				if( eleA_modAllHEEP || eleB_modAllHEEP )
					normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

			} // End-if: isInZMassRange
			else if(normDiEle_EB_fidECALDrFRPre.isAboveZMassRange()){
				normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				// Setup bools for whether eleA and eleB pass HEEPNoIso ...
				bool eleA_HEEPNoIso = normDiEle_EB_fidECALDrFRPre.eleA().ApplyHEEPCutsNoIso();
				bool eleB_HEEPNoIso = normDiEle_EB_fidECALDrFRPre.eleB().ApplyHEEPCutsNoIso();
				// Setup bools for whether eleA and eleB pass HEEPNoIso+modIso ...
				bool eleA_modAllHEEP = ( eleA_HEEPNoIso && ( normDiEle_EB_fidECALDrFRPre.eleA_modTrkIsolCut() && normDiEle_EB_fidECALDrFRPre.eleA_modEmHad1IsoCut() ) );
				bool eleB_modAllHEEP = ( eleB_HEEPNoIso && ( normDiEle_EB_fidECALDrFRPre.eleB_modTrkIsolCut() && normDiEle_EB_fidECALDrFRPre.eleB_modEmHad1IsoCut() ) );

				// Region C: Both electrons pass HEEP cuts
				if(eleA_HEEPNoIso && eleB_HEEPNoIso)
					normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				else // (Region D in bothModHEEP vs invMass)
					normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothHEEPNoIso_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);

				if(eleA_modAllHEEP && eleB_modAllHEEP)
					normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
				else // (Region D in bothModHEEP vs invMass)
					normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothModAllHEEP_Histos_.FillHistos(normDiEle_EB_fidECALDrFRPre, evtWeight);
			}
		}// End-if: trg_PathA_decision_
	}

	// QCD estimation (additional) code:
	// EB HEEP + incl. GSF di-electrons
	if( tsw::NumPassingCuts(normEles_EB_isFidEcalDrAndFRPreFlags)>1 && tsw::NumPassingCuts(normEles_EB_HEEPCutsNoIsoFlags)>0 ){
		// Find the indices of the two electrons ...
		unsigned int highestEtHEEPNoIsoEle_idx = IndicesOfElesPassingCuts(normEles_, normEles_EB_HEEPCutsNoIsoFlags, 0).at(0);
		unsigned int highestEtDifffidECALDrFRPreEle_idx = IndicesOfElesPassingCuts(normEles_, normEles_EB_isFidEcalDrAndFRPreFlags, 0).at(0);
		if(highestEtDifffidECALDrFRPreEle_idx==highestEtHEEPNoIsoEle_idx)
			highestEtDifffidECALDrFRPreEle_idx = IndicesOfElesPassingCuts(normEles_, normEles_EB_isFidEcalDrAndFRPreFlags, 0).at(1);

		bool HEEPNoIsoEle_higherEtEleOfPair = (normEles_.at(highestEtHEEPNoIsoEle_idx).et() > normEles_.at(highestEtDifffidECALDrFRPreEle_idx).et());

		// Form the di-electron object
		tsw::HEEPDiEle normDiEle_EB_HEEPNoIsoPLUSgsf = tsw::HEEPDiEle(normEles_.at(highestEtHEEPNoIsoEle_idx),
																						normEles_.at(highestEtDifffidECALDrFRPreEle_idx) );

		// Require signal trigger to have fired
		if(trg_PathA_decision_){
			// Apply Z mass window cut
			if(normDiEle_EB_HEEPNoIsoPLUSgsf.isInZMassRange()){

				normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_Histos_.FillHistos(normDiEle_EB_HEEPNoIsoPLUSgsf, evtWeight);

				// Now check to see if the 'HEEPNoIso' electron also passes the modified HEEP isol. cuts ...
				bool heepNoIsoEle_passesModHEEPIso;
				if(HEEPNoIsoEle_higherEtEleOfPair)
					heepNoIsoEle_passesModHEEPIso = ( normDiEle_EB_HEEPNoIsoPLUSgsf.eleA_modTrkIsolCut() && normDiEle_EB_HEEPNoIsoPLUSgsf.eleA_modEmHad1IsoCut() );
				else
					heepNoIsoEle_passesModHEEPIso = ( normDiEle_EB_HEEPNoIsoPLUSgsf.eleB_modTrkIsolCut() && normDiEle_EB_HEEPNoIsoPLUSgsf.eleB_modEmHad1IsoCut() );

				if(heepNoIsoEle_passesModHEEPIso)
					normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_modHEEPIsoOnHEEPNoIso_Histos_.FillHistos(normDiEle_EB_HEEPNoIsoPLUSgsf, evtWeight);

			} // End-if: isInZMassRange
		} // End-if: trg_PathA_decision_

	}

	//	//Modfied reconstruction di-electron pairs ...
	//	if(readInBstdEles_){
	//		if(tsw::NumPassingCuts(bstdEles_HEEPCutsFlags)>1){
	//			bstdDiEle_AllHEEP = tsw::HEEPDiEle::HEEPDiEle( bstdEles_, bstdEles_HEEPCutsFlags );
	//			bstdDiEle_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
	//			if( bstdDiEle_AllHEEP.invMass()>0.1 ){
	//				bstdDiEle_M0_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);}
	//			else{
	//				bstdDiEle_Mleq0_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);}
	//			if( fabs(bstdDiEle_AllHEEP.scDeltaEta())>0.001 || fabs(bstdDiEle_AllHEEP.scDeltaPhi()>0.001) )
	//				bstdDiEle_diffSC_AllHEEP_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
	//			if( bstdDiEle_AllHEEP.eleA().isEB() && bstdDiEle_AllHEEP.eleB().isEB() )
	//				bstdDiEle_AllHEEP_EBEB_Histos_.FillHistos(bstdDiEle_AllHEEP, evtWeight);
	//		}
	//
	//		if(tsw::NumPassingCuts(bstdEles_HEEPCutsNoIsoFlags)>1){
	//			bstdDiEle_HEEPNoIso = tsw::HEEPDiEle::HEEPDiEle( bstdEles_, bstdEles_HEEPCutsNoIsoFlags );
	//			bstdDiEle_HEEPNoIso_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//			if( bstdDiEle_HEEPNoIso.eleA().isEB() && bstdDiEle_HEEPNoIso.eleB().isEB() )
	//				bstdDiEle_HEEPNoIso_EBEB_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//			if( ( fabs(bstdDiEle_HEEPNoIso.scDeltaEta())>0.001 || fabs(bstdDiEle_HEEPNoIso.scDeltaPhi()>0.001) ) && bstdDiEle_HEEPNoIso.isInZMassRange() ){
	//				bstdDiEle_diffSC_HEEPNoIso_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//
	//				mcZboson_modReconNoIsoEvts_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
	//				if( fabs(mcZ_daughterA_p4_.Eta())<=1.442 && fabs(mcZ_daughterB_p4_.Eta())<=1.442 )
	//					mcZboson_modReconNoIsoEvts_EBEB_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
	//				else if( fabs(mcZ_daughterA_p4_.Eta())>=1.56 && fabs(mcZ_daughterB_p4_.Eta())>=1.56 )
	//					mcZboson_modReconNoIsoEvts_EEEE_Histos_.FillHistos(mcZboson_, mcZ_daughters_dR_, mcZ_daughters_dEta_, mcZ_daughters_dPhi_, mcZ_daughters_openingAngle_, evtWeight);
	//
	//				if(trg_PathA_decision_){
	//					bstdDiEle_diffSC_HEEPNoIso_5e32_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//					if(trg_PathA_highestTrigObjEt_ > 52.0)
	//						bstdDiEle_diffSC_HEEPNoIso_1e33_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//					if(trg_PathA_highestTrigObjEt_ > 65.0)
	//						bstdDiEle_diffSC_HEEPNoIso_2e33_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//					if(trg_PathA_highestTrigObjEt_ > 80.0)
	//						bstdDiEle_diffSC_HEEPNoIso_TrgEt80_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//					if(trg_PathA_highestTrigObjEt_ > 100.0)
	//						bstdDiEle_diffSC_HEEPNoIso_TrgEt100_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//					if(trg_PathA_highestTrigObjEt_ > 120.0)
	//						bstdDiEle_diffSC_HEEPNoIso_TrgEt120_Histos_.FillHistos(bstdDiEle_HEEPNoIso, evtWeight);
	//
	//				}
	//			}
	//		}
	//	} // End of if(readInBstdEles_)
}

//--------------------------------------------------------//
//---- Method for implementing emu analysis ...           //
void BstdZeeFirstAnalyser::DoEMuAnalysis()
{
	const Double_t evtWeight = eventHelper_.totWeight();

	//  Filling histograms with the properties of the highest pT muon (no cuts), and the highest pT tight Muon ...
	if(normMuons_.NumOfMuons()>0)
		normMuons_1stpT_Histos_.FillHistos( normMuons_.HighestPtMuon(), 1.0);
	if(normMuons_tight_.NumOfMuons()>0)
		normMuons_tight_1stpT_Histos_.FillHistos( normMuons_tight_.HighestPtMuon(), 1.0);
	if(normMuons_barrel_tight_.NumOfMuons()>0)
		normMuons_barrel_tight_1stpT_Histos_.FillHistos( normMuons_barrel_tight_.HighestPtMuon(), 1.0 );

	// Forming the emu object (IFF there is at least 1 HEEP non-isol electron, and at least one tight muon) ...
	if(normMuons_barrel_tight_.NumOfMuons()>0 && normEles_EB_HEEPNoIso_1stpT_exists_){
		eleMu_EB_HEEPNoIso_muB_tight_ = tsw::EleMuObject(normEles_EB_HEEPNoIso_1stpT_, normMuons_barrel_tight_.HighestPtMuon());

		// Check that ele-mu object is in Z mass range
		if(eleMu_EB_HEEPNoIso_muB_tight_.isInZMassRange()){
			h_eleMu_EB_HEEPNoIso_muB_tight_MZ_.FillHistos(eleMu_EB_HEEPNoIso_muB_tight_, evtWeight);

			// Check that e-mu trigger has fired
			if(eventHelper_.GetTrigInfo_eMuPath_decision()){
				h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_.FillHistos(eleMu_EB_HEEPNoIso_muB_tight_, evtWeight);
				//Now, setting the branch variables, and adding a new event to the tree ...
				eleMuTreeVar_mass_   = eleMu_EB_HEEPNoIso_muB_tight_.mass();
				eleMuTreeVar_pT_     = eleMu_EB_HEEPNoIso_muB_tight_.pT();
				eleMuTreeVar_weight_ = evtWeight;
				eleMuTree_->Fill();
				// Apply HEEP isol. cuts to electron
				if(eleMu_EB_HEEPNoIso_muB_tight_.GetElectron()->ApplyIsoVarHEEPCuts())
					h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_HEEPIso_.FillHistos(eleMu_EB_HEEPNoIso_muB_tight_, evtWeight);
			}
		}
	}

	//Settting the branch variables for the diEle tree - IFF a di-ele has been constructed ...
	if(normDiEle_HEEPNoIso_exists_){
		diEleTreeVar_mass_ = normDiEle_HEEPNoIso_.invMass();
		diEleTreeVar_pT_ = normDiEle_HEEPNoIso_.pT();
		diEleTreeVar_weight_ = evtWeight;
		diEleTree_->Fill();
	}
}

//--------------------------------------------------------//
//---- Methods for storing histograms in output file...   //
//--------------------------------------------------------//
void BstdZeeFirstAnalyser::SaveReconValidationHistos(TFile* histosFile)
{
	normEles_reconValidationHistos_.WriteHistos(histosFile);
	normEles_simpleCuts_reconValHistos_.WriteHistos(histosFile);
	normEles_HEEPCuts_reconValHistos_.WriteHistos(histosFile);
	bstdEles_reconValidationHistos_.WriteHistos(histosFile);
	bstdEles_simpleCuts_reconValHistos_.WriteHistos(histosFile);
	bstdEles_HEEPCuts_reconValHistos_.WriteHistos(histosFile);

	normDiEle_AllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_AllHEEP_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_AllHEEP_MZ_Histos_.WriteHistos(histosFile);
	normDiEle_AllHEEP_MZ_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_AllHEEP_MZ_trgA_Rgnl_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_MZ_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_MZ_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_HEEPNoIso_MZ_trgA_Rgnl_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_trkIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_EmHad1Iso_trkIso_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_notEmHad1Iso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_trgA_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIso_MZ_stdHEEPIso_Histos_.WriteHistos(histosFile);

	// Di-ele histograms for QCD estimation ...
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_bothModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_NOTbothModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_oneModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_zeroModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_MZ_geqOneModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothHEEPNoIso_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_bothModAllHEEP_Histos_.WriteHistos(histosFile);
	normDiEle_EB_fidECALDrFRPre_trgA_aboveMZ_NOTbothModAllHEEP_Histos_.WriteHistos(histosFile);

	normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_Histos_.WriteHistos(histosFile);
	normDiEle_EB_HEEPNoIsoPLUSgsf_trgA_MZ_modHEEPIsoOnHEEPNoIso_Histos_.WriteHistos(histosFile);

	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_.WriteHistos(histosFile);
	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_.WriteHistos(histosFile);
	h_eleMu_EB_HEEPNoIso_muB_tight_MZ_eMuTrg_HEEPIso_.WriteHistos(histosFile);

	normMuons_1stpT_Histos_.WriteHistos(histosFile);
	normMuons_tight_1stpT_Histos_.WriteHistos(histosFile);
	normMuons_barrel_tight_1stpT_Histos_.WriteHistos(histosFile);

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
	normHEEPEles_ecalEnergyError_ = 0;
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
	normHEEPEles_closestCtfTrk_pt_ = 0;
	normHEEPEles_closestCtfTrk_eta_ = 0;
	normHEEPEles_closestCtfTrk_phi_ = 0;
	normHEEPEles_closestCtfTrk_innerPt_ = 0;
	normHEEPEles_closestCtfTrk_innerEta_ = 0;
	normHEEPEles_closestCtfTrk_innerPhi_ = 0;
	normHEEPEles_closestCtfTrk_outerPt_ = 0;
	normHEEPEles_closestCtfTrk_outerEta_ = 0;
	normHEEPEles_closestCtfTrk_outerPhi_ = 0;

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

	normHEEPEles_SCposn_eta_ = 0;
	normHEEPEles_SCposn_phi_ = 0;
	normHEEPEles_SC_rawEnergy_ = 0;
	normHEEPEles_SC_recHits_Et_  = 0;
	normHEEPEles_SC_recHits_eta_ = 0;
	normHEEPEles_SC_recHits_phi_ = 0;
	normHEEPEles_SC_recHits_isFromEB_ = 0;
	normHEEPEles_SC_totEnergyRecHits_ = 0;
	normHEEPEles_SC_totNumRecHits_ = 0;

	normHEEPEles_gsfTrk_eta_ = 0;
	normHEEPEles_gsfTrk_phi_ = 0;
	normHEEPEles_gsfTrk_vz_ = 0;
	normHEEPEles_innerIsoConeTrks_pt_ = 0;
	normHEEPEles_innerIsoConeTrks_eta_ = 0;
	normHEEPEles_innerIsoConeTrks_phi_ = 0;
	normHEEPEles_innerIsoConeTrks_vz_ = 0;

	normHEEPEles_numMissInnerHits_ = 0;

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
	inputFile_tree_->SetBranchAddress("event", &event_); // inputFile_tree_->SetBranchStatus("event",0); //
	//event_ = &dummyEvent_;
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
	else{
		inputFile_tree_->SetBranchStatus("mcEles_number",0);
		inputFile_tree_->SetBranchStatus("mcEles_HighestEt_*",0);
		inputFile_tree_->SetBranchStatus("mcEles_2ndHighestEt_*",0);
		inputFile_tree_->SetBranchStatus("mcZcandidate_*",0);
		inputFile_tree_->SetBranchStatus("mcZboson_*",0);
	}


	inputFile_tree_->SetBranchStatus("normGsfEles_*",0);
	inputFile_tree_->SetBranchStatus("normGsfEles_number",1);
	inputFile_tree_->SetBranchStatus("normGsfEles_charge",1);
	//Setting up the pointer links for the 'normal' GSF electron branches ...
	inputFile_tree_->SetBranchAddress("normGsfEles_number", &normGsfEles_number_); //unsigned int
	//inputFile_tree_->SetBranchAddress("normGsfEles_p4ptr_", &normGsfEles_OLDp4_  ); // std::vector<ROOT::Math::XYZTVector>*
	inputFile_tree_->SetBranchAddress("normGsfEles_charge", &normGsfEles_charge_); //std::vector<Int_t>*
	if(false){
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
	}
	//*** Branches containing heep::Ele method values for norm eles ***
	//kinematic and geometric methods
	inputFile_tree_->SetBranchAddress("normHEEPEles_et", &normHEEPEles_et_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfEt", &normHEEPEles_gsfEt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scEt", &normHEEPEles_scEt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_energy", &normHEEPEles_energy_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_gsfEnergy",0);									//	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfEnergy", &normHEEPEles_gsfEnergy_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_caloEnergy", &normHEEPEles_caloEnergy_);
//	inputFile_tree_->SetBranchAddress("normHEEPEles_ecalEnergyError", &normHEEPEles_ecalEnergyError_);			/* TEMP v1f/g FIX */
	inputFile_tree_->SetBranchAddress("normHEEPEles_eta", &normHEEPEles_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scEta", &normHEEPEles_scEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_detEta",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_detEta", &normHEEPEles_detEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_detEtaAbs",0);									//	inputFile_tree_->SetBranchAddress("normHEEPEles_detEtaAbs", &normHEEPEles_detEtaAbs_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_phi", &normHEEPEles_phi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_scPhi", &normHEEPEles_scPhi_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_detPhi",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_detPhi", &normHEEPEles_detPhi_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_zVtx",0);											//	inputFile_tree_->SetBranchAddress("normHEEPEles_zVtx", &normHEEPEles_zVtx_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_p4ptr", &normHEEPEles_p4ptr_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfP4ptr", &normHEEPEles_gsfP4ptr_);

	//classification (couldnt they have just named it 'type')
	inputFile_tree_->SetBranchStatus("normHEEPEles_classification",0);							//	inputFile_tree_->SetBranchAddress("normHEEPEles_classification", &normHEEPEles_classification_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEcalDriven", &normHEEPEles_isEcalDriven_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_isTrackerDriven",0);							//	inputFile_tree_->SetBranchAddress("normHEEPEles_isTrackerDriven", &normHEEPEles_isTrackerDriven_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEB", &normHEEPEles_isEB_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isEE", &normHEEPEles_isEE_);

	//track methods
	inputFile_tree_->SetBranchAddress("normHEEPEles_charge", &normHEEPEles_charge_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_trkCharge",0);									//	inputFile_tree_->SetBranchAddress("normHEEPEles_trkCharge", &normHEEPEles_trkCharge_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_pVtx",0);											//	inputFile_tree_->SetBranchAddress("normHEEPEles_pVtx", &normHEEPEles_pVtx_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_pCalo",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_pCalo", &normHEEPEles_pCalo_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_ptVtx", &normHEEPEles_ptVtx_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_ptCalo", &normHEEPEles_ptCalo_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_pt",  &normHEEPEles_closestCtfTrk_pt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_eta", &normHEEPEles_closestCtfTrk_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_phi", &normHEEPEles_closestCtfTrk_phi_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerPt",0);  // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerPt",  &normHEEPEles_closestCtfTrk_innerPt_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerEta",0); // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerEta", &normHEEPEles_closestCtfTrk_innerEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerPhi",0); // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerPhi", &normHEEPEles_closestCtfTrk_innerPhi_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerPt",0);  // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerPt",  &normHEEPEles_closestCtfTrk_outerPt_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerEta",0); // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerEta", &normHEEPEles_closestCtfTrk_outerEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerPhi",0); // inputFile_tree_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerPhi", &normHEEPEles_closestCtfTrk_outerPhi_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	inputFile_tree_->SetBranchAddress("normHEEPEles_hOverE", &normHEEPEles_hOverE_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_dEtaIn", &normHEEPEles_dEtaIn_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_dPhiIn", &normHEEPEles_dPhiIn_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_dPhiOut",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_dPhiOut", &normHEEPEles_dPhiOut_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_epIn", &normHEEPEles_epIn_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_epOut", &normHEEPEles_epOut_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_fbrem",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_fbrem", &normHEEPEles_fbrem_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_bremFrac",0);									//	inputFile_tree_->SetBranchAddress("normHEEPEles_bremFrac", &normHEEPEles_bremFrac_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_invEOverInvP", &normHEEPEles_invEOverInvP_);

	//shower shape variables
	inputFile_tree_->SetBranchStatus("normHEEPEles_sigmaEtaEta",0);								//	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaEtaEta", &normHEEPEles_sigmaEtaEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_sigmaEtaEtaUnCorr",0);						//	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaEtaEtaUnCorr", &normHEEPEles_sigmaEtaEtaUnCorr_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_sigmaIEtaIEta", &normHEEPEles_sigmaIEtaIEta_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_e1x5",0);											//	inputFile_tree_->SetBranchAddress("normHEEPEles_e1x5", &normHEEPEles_e1x5_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_e2x5Max",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_e2x5Max", &normHEEPEles_e2x5Max_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_e5x5",0);											//	inputFile_tree_->SetBranchAddress("normHEEPEles_e5x5", &normHEEPEles_e5x5_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e1x5Over5x5", &normHEEPEles_e1x5Over5x5_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_e2x5MaxOver5x5", &normHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolEm", &normHEEPEles_isolEm_);
	inputFile_tree_->SetBranchStatus("normHEEPEles_isolHad",0);										//	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHad", &normHEEPEles_isolHad_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHadDepth1", &normHEEPEles_isolHadDepth1_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolHadDepth2", &normHEEPEles_isolHadDepth2_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolPtTrks", &normHEEPEles_isolPtTrks_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_isolEmHadDepth1", &normHEEPEles_isolEmHadDepth1_);

//	inputFile_tree_->SetBranchAddress("normHEEPEles_SCposn_eta", &normHEEPEles_SCposn_eta_);								/* TEMP v1f/g FIX */
//	inputFile_tree_->SetBranchAddress("normHEEPEles_SCposn_phi", &normHEEPEles_SCposn_phi_);								/* TEMP v1f/g FIX */
//	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_rawEnergy", &normHEEPEles_SC_rawEnergy_);							/* TEMP v1f/g FIX */
	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_recHits_Et",  &normHEEPEles_SC_recHits_Et_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_recHits_eta", &normHEEPEles_SC_recHits_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_recHits_phi", &normHEEPEles_SC_recHits_phi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_recHits_isFromEB", &normHEEPEles_SC_recHits_isFromEB_);
//	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_totEnergyRecHits", &normHEEPEles_SC_totEnergyRecHits_);		/* TEMP v1f/g FIX */
//	inputFile_tree_->SetBranchAddress("normHEEPEles_SC_totNumRecHits", &normHEEPEles_SC_totNumRecHits_);				/* TEMP v1f/g FIX */

	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfTrk_eta", &normHEEPEles_gsfTrk_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfTrk_phi", &normHEEPEles_gsfTrk_phi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_gsfTrk_vz", &normHEEPEles_gsfTrk_vz_);

	inputFile_tree_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_pt",  &normHEEPEles_innerIsoConeTrks_pt_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_eta", &normHEEPEles_innerIsoConeTrks_eta_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_phi", &normHEEPEles_innerIsoConeTrks_phi_);
	inputFile_tree_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_vz",  &normHEEPEles_innerIsoConeTrks_vz_);

	inputFile_tree_->SetBranchAddress("normHEEPEles_numMissInnerHits", &normHEEPEles_numMissInnerHits_);				/* TEMP v1f/g FIX */

	/////////////////////////////////////////////////////////////////////////////////////////
	//Setting up the pointer links for the special reconstruction GSF electron branches ...
	if(readInBstdEles_){
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
	}
	else
		inputFile_tree_->SetBranchStatus("bstd*",0);

	std::cout << "   pointer links for branch read-out now set up." << std::endl;
}

//==============================================================
void BstdZeeFirstAnalyser::SetupEMuMethodTrees(){
	eleMuTree_ = new TTree("myTree","eleMu data");
	eleMuTree_->SetDirectory(0); // This line is needed as a 'QUICK FIX' to stop the following error when running over very large nos. of events ...
	/* Error is as follows:
	 * Error in <TTree::Fill>: Failed filling branch:myTree.mass, nbytes=-1, entry=3990
	 *  This error is symptomatic of a Tree created as a memory-resident Tree
	 *  Instead of doing:
	 *     TTree *T = new TTree(...)
	 *     TFile *f = new TFile(...)
	 *  you should do:
	 *     TFile *f = new TFile(...)
	 *     TTree *T = new TTree(...)
	 */
	eleMuTree_->Branch("mass",   &eleMuTreeVar_mass_,   "mass/D");
	eleMuTree_->Branch("pT",     &eleMuTreeVar_pT_,     "pT/D");
	eleMuTree_->Branch("weight", &eleMuTreeVar_weight_, "weight/D");

	diEleTree_ = new TTree("myTree","diEle data");
	diEleTree_->SetDirectory(0);// This line is needed as a 'QUICK FIX' to stop the following error when running over very large nos. of events ...
	/* Error is as follows:
	 * Error in <TTree::Fill>: Failed filling branch:myTree.mass, nbytes=-1, entry=3990
	 *  This error is symptomatic of a Tree created as a memory-resident Tree
	 *  Instead of doing:
	 *     TTree *T = new TTree(...)
	 *     TFile *f = new TFile(...)
	 *  you should do:
	 *     TFile *f = new TFile(...)
	 *     TTree *T = new TTree(...)
	 */
	diEleTree_->Branch("mass",   &diEleTreeVar_mass_,   "mass/D");
	diEleTree_->Branch("pT",   &diEleTreeVar_pT_,   "pT/D");
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


//		TString inFilePrefix = "/home/ppd/nnd85574/Work/BstdZee/CMSSW_4_2_8_patch7/src/NTupler/BstdZeeNTupler/";
		TString inFilePrefix = "/opt/ppd/newscratch/williams/Datafiles/NTuples/42Xv1x/ModIsoStudies_2012-05-07/";
		TString outFilePrefix = "/opt/ppd/newscratch/williams/Datafiles/AnaTuples/ModIsoStudies_2012-05-16tests/";
		std::vector<BstdZeeFirstAnalyser> myAnalysers; myAnalysers.clear();
		std::vector<TString> inFile; inFile.clear();
		std::vector<TString> outFile; outFile.clear();
		std::vector<bool> isMCflag; isMCflag.clear();
		std::vector<int>  nEvents; nEvents.clear();
		std::vector<Double_t> intLumiPerEvent; intLumiPerEvent.clear(); // in inv fb !!
		Double_t desiredIntLumi = 3.917;  // in inv fb !!

		std::vector<std::string> giMasses;
		giMasses.push_back("2000");

		for(std::vector<std::string>::const_iterator giMassIt = giMasses.begin(); giMassIt != giMasses.end(); giMassIt++){
			inFile.push_back("mcNTuple_42Xv1x_QstarGI-M"+(*giMassIt)+"--Fall11-PUS6_2012-05-07.root"); isMCflag.push_back(true);
			outFile.push_back("QstarGI-M"+(*giMassIt)+"--Fall11-PUS6_2012-05-16test1");
			nEvents.push_back( -1 );
			intLumiPerEvent.push_back( -1.0 ); // in inv fb
		}

		std::vector<std::string> ciMasses;

		for(std::vector<std::string>::const_iterator ciMassIt = ciMasses.begin(); ciMassIt != ciMasses.end(); ciMassIt++){
			inFile.push_back("mcNTuple_42Xv1x_QstarCI-M"+(*ciMassIt)+"--Fall11-PUS6_2012-05-16test1.root"); isMCflag.push_back(true);
			outFile.push_back("QstarCI-M"+(*ciMassIt)+"--Fall11-PUS6_2012-05-07");
			nEvents.push_back( -1 );
			intLumiPerEvent.push_back( -1.0 ); // in inv fb
		}

		////////////////////
		// Data - Photon

//		// May10ReReco
//		inFile.push_back("42Xv1f/data/Photon-May10ReReco_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-May10ReReco_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Aug05 ReReco
//		inFile.push_back("42Xv1f/data/Photon-Aug05ReReco_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-Aug05ReReco_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Prompt-v4
//		inFile.push_back("42Xv1f/data/Photon-PromptReco-v4_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-PromptReco-v4_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Prompt-v6
//		inFile.push_back("42Xv1f/data/Photon-PromptReco-v6_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-PromptReco-v6_2011-10-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Run11B-Prompt-v1
//		inFile.push_back("42Xv1f/data/Photon-Run11B-PromptReco-v1_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-Run11B-PromptReco-v1_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb

//		// full 42Xv1g Run11B-Prompt-v1
//		inFile.push_back("42Xv1g/data/Photon-Run11B-PromptReco-v1_NTuple-42Xv1g_2011-11-23.root"); isMCflag.push_back(false);
//		outFile.push_back("Photon-FULL-Run11B-PromptReco-v1_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb

//		////////////////////
//		// Data - MuEG
//
//		// May10ReReco
//		inFile.push_back("42Xv1f/data/MuEG-May10ReReco_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-May10ReReco_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Aug05 ReReco
//		inFile.push_back("42Xv1f/data/MuEG-Aug05ReReco_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-Aug05ReReco_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Prompt-v4
//		inFile.push_back("42Xv1f/data/MuEG-PromptReco-v4_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-PromptReco-v4_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Prompt-v6
//		inFile.push_back("42Xv1f/data/MuEG-PromptReco-v6_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-PromptReco-v6_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// Run11B-Prompt-v1
//		inFile.push_back("42Xv1f/data/MuEG-Run11B-PromptReco-v1_NTuple-42Xv1f_2011-10-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-Run11B-PromptReco-v1_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb
//
//		// full 42Xv1g Run11B-Prompt-v1
//		inFile.push_back("42Xv1g/data/MuEG-Run11B-PromptReco-v1_NTuple-42Xv1g_2011-11-23.root"); isMCflag.push_back(false);
//		outFile.push_back("MuEG-FULL-Run11B-PromptReco-v1_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1 ); // in inv fb


//
//		////////////////////
//		// Background MC
//
//		// DYJetsToLL
//		inFile.push_back("42Xv1e/mc/DYJetsToLL_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/2475000.0 ); // in inv fb
//
//		// DYJetsToLL-ee
//		inFile.push_back("42Xv1e/mc/DYJetsToLL-ee_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL-ee_2011-10-23");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/(2475000.0/3.0) ); // in inv fb
//
//		// DYJetsToLL-ee-lowMCZpT
//		inFile.push_back("42Xv1e/mc/DYJetsToLL-ee_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL-ee-lowMCZpT_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/(2475000.0/3.0) ); // in inv fb
//
//		// DYJetsToLL-ZpT100-ee
//		inFile.push_back("42Xv1e/mc/DYJetsToLL-ZpT100-ee_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL-ZpT100-ee_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/(25100.0/3.0) ); // in inv fb
//
//		// DYJetsToLL-ZpT100-ee-highMCZpT
//		inFile.push_back("42Xv1e/mc/DYJetsToLL-ZpT100-ee_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL-ZpT100-ee-highMCZpT_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/(25100.0/3.0) ); // in inv fb

//		// DYToEE-powheg-M20
//		inFile.push_back("42Xv1e/mc/DYToEE-powheg-M20_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYToEE-powheg-M20_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/1628000.0 ); // in inv fb
//
//		// DYJetsToLL-TauTau
//		inFile.push_back("42Xv1e/mc/DYJetsToLL-TauTau_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(true);
//		outFile.push_back("DYJetsToLL-TauTau_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/(2475000.0/3.0) ); // in inv fb
//
//		// TTbar
//		inFile.push_back("42Xv1e/mc/TTbar-pythia_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("TTbar-pythia_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/157500.0 ); // in inv fb // intLumiPerEvent.push_back( 1.0/94000.0 ); // in inv fb
//
//		// TTbarJets
//		inFile.push_back("42Xv1e/mc/TTbarJets_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("TTbarJets_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/157500.0 ); // in inv fb // intLumiPerEvent.push_back( 1.0/94760.0 ); // in inv fb
//
//		// TW
//		inFile.push_back("42Xv1e/mc/TW_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("TW_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/7870.0 ); // in inv fb // intLumiPerEvent.push_back( 1.0/7466.0 ); // in inv fb
//
//		// TbarW
//		inFile.push_back("42Xv1e/mc/TbarW_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("TbarW_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/7870.0 ); // in inv fb // intLumiPerEvent.push_back( 1.0/7460.0 ); // in inv fb
//
//		// WWTo2L2Nu
//		inFile.push_back("42Xv1e/mc/WWTo2L2Nu_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("WWTo2L2Nu_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/2927.0 ); // in inv fb
//
//		// WZTo3LNu
//		inFile.push_back("42Xv1e/mc/WZTo3LNu_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("WZTo3LNu_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/339.4 ); // in inv fb
//
//		// ZGamma
//		inFile.push_back("42Xv1e/mc/ZGamma_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("ZGamma_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/34160.0 ); // in inv fb
//
//		// ZZ
//		inFile.push_back("42Xv1e/mc/ZZ_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("ZZ_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/4287.0 ); // in inv fb
//
//		// QCD-dblEM
//		inFile.push_back("42Xv1e/mc/QCD-dblEM_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("QCD-dblEM_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/43571000.0 ); // in inv fb

//		// WJetsToLNu
//		inFile.push_back("42Xv1e/mc/WJetsToLNu_NTuple-42Xv1e_2011-10-14.root"); isMCflag.push_back(false);
//		outFile.push_back("WJets_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/27770000.0 ); // in inv fb


//		////////////////////
//		// Signal MC
//
//		// SigMC - 0.75TeVq*
//		inFile.push_back("42Xv1f/local/0-75TeVq_NTuple-42Xv1f_2011-11-03.root"); isMCflag.push_back(true);
//		outFile.push_back("0-75TeVq_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/925.0 ); // in inv fb
//
//		// SigMC - 1.00TeVq*
//		inFile.push_back("42Xv1f/local/1-00TeVq_NTuple-42Xv1f_2011-11-03.root"); isMCflag.push_back(true);
//		outFile.push_back("1-00TeVq_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/110.0 ); // in inv fb
//
//		// SigMC - 2.00TeVq*
//		inFile.push_back("42Xv1f/local/2-00TeVq_NTuple-42Xv1f_2011-11-03.root"); isMCflag.push_back(true);
//		outFile.push_back("2-00TeVq_2011-11-22");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( 1.0/1.2 ); // in inv fb


//		// SigMC - 1.00TeVq*
//		inFile.push_back("../nextVersion/1-00TeVq_NTuple-42Xv1g-Pre_2011-11-17.root"); isMCflag.push_back(true);
////		inFile.push_back("testNTuple.root"); isMCflag.push_back(true);
//		outFile.push_back("histos_1-00TeVq_2011-11-17");
//		nEvents.push_back( -1 );
//		intLumiPerEvent.push_back( -1.0 ); // intLumiPerEvent.push_back( 1.0/110.0 ); // in inv fb


//		// SigMC - 0.75TeVu*
//		inFile.push_back("sigMC_0-75TeVu_2e33NTuple_42Xv1d_2kEvts_2011-08-30.root"); isMCflag.push_back(false);
//		outFile.push_back("histos_0-75TeVu_2011-08-30_2e33");
//		nEvents.push_back( 2000 );
//		intLumiPerEvent.push_back( -1.0 ); //intLumiPerEvent.push_back( 1.0/925.0 ); //intLumiPerEvent.push_back( -1.0 ); // in inv fb		// SigMC - 0.75TeVu*
//		inFile.push_back("sigMC_0-75TeVu_2e33noCaloIdL-NTuple_42Xv1d_2kEvts_2011-08-30.root"); isMCflag.push_back(false);
//		outFile.push_back("histos_0-75TeVu_2011-08-30_2e33noCaloIdL");
//		nEvents.push_back( 2000 );
//		intLumiPerEvent.push_back( -1.0 ); //intLumiPerEvent.push_back( 1.0/925.0 ); //intLumiPerEvent.push_back( -1.0 ); // in inv fb
//		// SigMC - 2.00TeVu*
//		inFile.push_back("sigMC_2-00TeVu_2e33NTuple_42Xv1d_2kEvts_2011-08-30.root"); isMCflag.push_back(false);
//		outFile.push_back("histos_2-00TeVu_2011-08-30_2e33");
//		nEvents.push_back( 2000 );
//		intLumiPerEvent.push_back( -1.0 ); //intLumiPerEvent.push_back( 1.0/1.2 ); //intLumiPerEvent.push_back( -1.0 ); // in inv fb



//		// Test
//		BstdZeeFirstAnalyser myAnalyser = BstdZeeFirstAnalyser(0, 10, true, "../../NTupler/BstdZeeNTupler/testNTuple_2-00TeVu.root", "testOut.root", 6, 80, 55, 1100.0);
//		myAnalyser.DoAnalysis(  1.0  );

		// Doing the analysis ...
		for(unsigned int idx=0; idx<inFile.size(); idx++){
			TStopwatch watch;
			watch.Start();
			std::cout << std::endl << "***--- NEW FILE (no " << idx+1 << "/" << inFile.size() << ") ---***" << std::endl;
			std::cout << " * Input:  " << inFile.at(idx) << std::endl;
			std::cout << " * Output: " << outFile.at(idx) << std::endl;
			BstdZeeFirstAnalyser myAnalysers = BstdZeeFirstAnalyser(0, nEvents.at(idx), isMCflag.at(idx), inFilePrefix+inFile.at(idx), outFilePrefix+outFile.at(idx), -1, 80, 55, 1100.0);
			std::cout << " * (" << myAnalysers.GetNumEvtsRunOver() << " events will be run over; " << nEvents.at(idx)  << " requested.)" << std::endl;

			if(outFile.at(idx).Contains("-highMCZpT_")){
				myAnalysers.SkipEvtIfMCZpTBelow(160.0);
				std::cout << "   (Low pT events will be skipped!)" << std::endl;
			}
			if(outFile.at(idx).Contains("-lowMCZpT_")){
				myAnalysers.SkipEvtIfMCZpTAbove(160.0);
				std::cout << "   (High pT events will be skipped!)" << std::endl;
			}
			std::cout << " Running the DoAnalysis method ..." << std::endl;
//			if(inFile.at(idx).find("mcZpTleq100GeV")!=inFile.at(idx).length()){
//				myAnalysers.RemoveEventsWithMcZpTGreaterThan(100.0);
//			}
			if(intLumiPerEvent.at(idx)>=0)
				myAnalysers.DoAnalysis(  desiredIntLumi/( static_cast<double>(myAnalysers.GetNumEvtsRunOver())*intLumiPerEvent.at(idx))  ); //myAnalyser.DoAnalysis( 1.0*(2321000.0/2530516.0) );
			else
				myAnalysers.DoAnalysis(  1.0  ); //myAnalyser.DoAnalysis( 1.0*(2321000.0/2530516.0) );
			//delete myAnalysers;
			//myAnalysers = NULL;
			watch.Stop();
			std::cout << " * " << myAnalysers.GetNumEvtsRunOver() << " events analysed in " << watch.RealTime() << " seconds (" << watch.CpuTime() << " sec of CPU time)." << std::endl;
			std::cout << " * [i.e. " << 1000000.0*watch.RealTime()/static_cast<double>(myAnalysers.GetNumEvtsRunOver()) << "sec/million events. (" << 1000000.0*watch.CpuTime()/static_cast<double>(nEvents.at(idx)) << " sec/event in CPU time.)]" << std::endl;
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


