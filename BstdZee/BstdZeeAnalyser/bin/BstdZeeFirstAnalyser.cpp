
// Standard C/C++ includes
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

// Boost includes
#include "boost/program_options/options_description.hpp"
#include "boost/program_options/variables_map.hpp"
#include "boost/program_options/parsers.hpp"

// ROOT includes

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
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
#include "TSWilliams/BstdZeeAnalyser/interface/tswEventHelper.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswHEEPEle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswReconValidationHistos.h"
#include "TSWilliams/BstdZeeNTupler/interface/tswUsefulFunctions.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswHEEPDiEle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswDiEleDistns.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswDiEleDistnsByRgn.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswMCParticle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswMCParticleDistns.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswMuonDistns.h"

#include "TSWilliams/BstdZeeAnalyser/interface/tswEleMuObject.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswEleMuDistns.h"


#ifdef __MAKECINT__
#pragma link C++ class vector<Int_t>+;
#pragma link C++ class vector< vector<float> >;
#endif

//Functions...


//=====================================================================================================
//-----------------------------------------------------------------------------------------------------
namespace tsw{

	std::vector<std::string> AddPrefixesToFileNamesVec(const std::vector<std::string>& origFileNamesVec){
		std::vector<std::string> newFileNamesVec;

		for(std::vector<std::string>::const_iterator fileNameIt=origFileNamesVec.begin(); fileNameIt<origFileNamesVec.end(); fileNameIt++){
			if(fileNameIt->find("/pnfs/pp.rl.ac.uk")==0)
				newFileNamesVec.push_back("dcap://heplnx203.pp.rl.ac.uk" + (*fileNameIt) );
			else if(fileNameIt->find("/store/")==0)
				newFileNamesVec.push_back("dcap://heplnx203.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms" + (*fileNameIt) );
			else
				newFileNamesVec.push_back(*fileNameIt);
		}

		return newFileNamesVec;
	}

	std::string TimeInfoString(TStopwatch& theStopwatch){
		std::string theString;

		std::ostringstream stream_cpuTime;
		stream_cpuTime << theStopwatch.CpuTime();

		std::ostringstream stream_realTime;
		stream_realTime << theStopwatch.RealTime();

		theString  = stream_cpuTime.str() +" secs (CPU);  ";
		theString += stream_realTime.str()+" secs (real)";
		return theString;
	}

	/*!
	 * \class TreeHandlerBase
	 * \brief Base class for tree handlers that implements all common functionality: saveToFile method, destruction of all trees
	 */
	class TreeHandlerBase{
	public:
		// CTOR & DTOR
		TreeHandlerBase(const std::string& mainTreeName, const std::string& mainTreeDescription, const std::string& fileName);
		virtual ~TreeHandlerBase();

		// PUBLIC METHODS
		void setEventCounter(unsigned int nEvtsTotal){
			totNumEvtsAnalyzed_ = nEvtsTotal;  numEvtsPass_ = mainAnaTree_->GetEntries(); }
		void saveToFile();

	private:
		// PRIVATE MEMBERS
		TFile* outputFile_;
		TTree* eventCountTree_;

		UInt_t numEvtsPass_;
		UInt_t totNumEvtsAnalyzed_;

	protected:
		// PROTECTED MEMBERS
		TTree* mainAnaTree_;
	};

	TreeHandlerBase::TreeHandlerBase(const std::string& mainTreeName, const std::string& mainTreeDescription, const std::string& fileName) :
		outputFile_( new TFile(fileName.c_str(), "RECREATE") ),
		eventCountTree_( new TTree("eventCountTree", "Tree containing orig event counts per file.") ),
		numEvtsPass_(0), totNumEvtsAnalyzed_(0),
		mainAnaTree_( new TTree(mainTreeName.c_str(), mainTreeDescription.c_str()) )
	{
		mainAnaTree_->SetDirectory( outputFile_ );
		eventCountTree_->SetDirectory( outputFile_ );

		eventCountTree_->Branch("nEvtsPass", &numEvtsPass_, "nEvtsPass/i");
		eventCountTree_->Branch("nEvtsRunOver", &totNumEvtsAnalyzed_, "nEvtsRunOver/i");
	}
	TreeHandlerBase::~TreeHandlerBase()
	{
		delete outputFile_; // TTrees themselves should not be deleted since they reside in the file (??)
	}

	void TreeHandlerBase::saveToFile()
	{
		// Firstly fill eventCountTree_
		eventCountTree_->Fill();

		// Then, write TTrees to file ...
		outputFile_->cd();
		mainAnaTree_->Write();
		eventCountTree_->Write();
		outputFile_->Close();
	}


	class ABCDMethodTree : public TreeHandlerBase {
	private:
		// PRIVATE MEMBERS

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
		ABCDMethodTree(const std::string& fileName) :
			TreeHandlerBase("abcdTree","ABCD method di-ele data", fileName)
		{
			// Setting up the event / di-ele branches ...
			mainAnaTree_->Branch("trgDecision", &treeVar_trgDecision_, "trgDecision/O");
			mainAnaTree_->Branch("trgName",     &treeVar_trgNamePtr_);
			treeVar_trgNamePtr_ = &treeVar_trgName_;
			mainAnaTree_->Branch("weight",      &treeVar_weight_, "weight/D");

			treeVar_p4Ptr_ = &treeVar_p4_;  mainAnaTree_->Branch("diEle_p4",   &treeVar_p4Ptr_);
			mainAnaTree_->Branch("diEle_mass", &treeVar_mass_,   "mass/D");
			mainAnaTree_->Branch("diEle_pT",   &treeVar_pT_,     "pT/D");
			mainAnaTree_->Branch("diEle_passHEEPNoIso",  &treeVar_passHEEPNoIso_,  "diEle_passHEEPNoIso/O");
			mainAnaTree_->Branch("diEle_passModHEEPIso", &treeVar_passModHEEPIso_, "diEle_passModHEEPIso/O");
			mainAnaTree_->Branch("diEle_passModTrkIso", &treeVar_passModTrkIso_, "diEle_passModTrkIso/O");
			mainAnaTree_->Branch("diEle_passModEmHad1Iso", &treeVar_passModEmHad1Iso_, "diEle_passModEmHad1Iso/O");

			// Setting up the branches for each electron ...
			treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;  mainAnaTree_->Branch("eleA_p4",     &treeVar_eleA_p4Ptr_);
			mainAnaTree_->Branch("eleA_charge",   &treeVar_eleA_charge_,   "eleA_charge/I"); //Int_t
			mainAnaTree_->Branch("eleA_EOverP", &treeVar_eleA_EOverP_, "eleA_EOverP/D"); //Double_t
			mainAnaTree_->Branch("eleA_passHEEPNoIso", &treeVar_eleA_passHEEPNoIso_, "eleA_passHEEPNoIso/O");

			treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;   mainAnaTree_->Branch("eleB_p4",     &treeVar_eleB_p4Ptr_);
			mainAnaTree_->Branch("eleB_charge",   &treeVar_eleB_charge_,   "eleB_charge/I"); //Int_t
			mainAnaTree_->Branch("eleB_EOverP", &treeVar_eleB_EOverP_, "eleB_EOverP/D"); //Double_t
			mainAnaTree_->Branch("eleB_passHEEPNoIso", &treeVar_eleB_passHEEPNoIso_, "eleB_passHEEPNoIso/O");
		}
		~ABCDMethodTree(){}

		void fillTree(tsw::HEEPDiEle* diEle, std::string trgName, bool trgDecision, double evtWeight){
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
			mainAnaTree_->Fill();
		}
	};

	///////////////////////////////////////////////////////////////////////////////
	// DiEleTree: Simple class for generating trees containing info about di-ele Z boson candidates
	//
	class DiEleTree : public TreeHandlerBase {
	private:
		// PRIVATE MEMBERS

		// For the event & kinematic branches ...
		Double_t treeVar_weight_;
		Double_t treeVar_genWeight_;
		Double_t treeVar_puWeight_;
		Double_t treeVar_xsecWeight_;

		UInt_t treeVar_runNum_;
		UInt_t treeVar_lumiNum_;
		UInt_t treeVar_evtNum_;

		Bool_t treeVar_trgDecision_;
		UInt_t treeVar_nVtx_;

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
		DiEleTree(const std::string& treeName, const std::string& fileName) :
			TreeHandlerBase(treeName, "Tree of Z candidates ("+treeName+")", fileName)
		{
			// Setting up the event / di-ele branches ...
			mainAnaTree_->Branch("weight",     &treeVar_weight_,     "weight/D");
			mainAnaTree_->Branch("genWeight",  &treeVar_genWeight_,  "genWeight/D");
			mainAnaTree_->Branch("puWeight",   &treeVar_puWeight_,   "puWeight/D");
			mainAnaTree_->Branch("xsecWeight", &treeVar_xsecWeight_, "xsecWeight/D");

			mainAnaTree_->Branch("run",    &treeVar_runNum_,   "run/i");
			mainAnaTree_->Branch("lumi",   &treeVar_lumiNum_,   "lumi/i");
			mainAnaTree_->Branch("evtNum", &treeVar_evtNum_,   "evtNum/i");

			mainAnaTree_->Branch("trgDecision", &treeVar_trgDecision_, "trgDecision/O");
			mainAnaTree_->Branch("nVtx", &treeVar_nVtx_, "nVtx/i");

			treeVar_Zp4Ptr_ = &treeVar_Zp4_; mainAnaTree_->Branch("Zp4",   &treeVar_Zp4Ptr_);
			mainAnaTree_->Branch("ZpT",   &treeVar_ZpT_,     "ZpT/D");
			mainAnaTree_->Branch("Zmass", &treeVar_Zmass_,   "Zmass/D");
			mainAnaTree_->Branch("dR",    &treeVar_dR_,      "dR/D");
			mainAnaTree_->Branch("dEta",  &treeVar_dEta_,    "dEta/D");
			mainAnaTree_->Branch("dPhi",  &treeVar_dPhi_,    "dPhi/D");

			treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;  mainAnaTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);
			treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;  mainAnaTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);

			mainAnaTree_->Branch("eleA_modHeepCutCode", &treeVar_eleA_modHeepCutCode_, "eleA_modHeepCutCode/I");
			mainAnaTree_->Branch("eleB_modHeepCutCode", &treeVar_eleB_modHeepCutCode_, "eleB_modHeepCutCode/I");
		}
		~DiEleTree(){}

		void fillTree(const tsw::HEEPDiEle& diEle, const tsw::EventHelper& evtHelper, bool trigDecision)
		{
			treeVar_weight_     = evtHelper.totWeight();
			treeVar_genWeight_  = evtHelper.genWeight();
			treeVar_puWeight_   = evtHelper.puWeight();
			treeVar_xsecWeight_ = evtHelper.xsecWeight();

			treeVar_runNum_ = evtHelper.runNum();
			treeVar_lumiNum_ = evtHelper.lumiSec();
			treeVar_evtNum_ = evtHelper.eventNum();

			treeVar_trgDecision_ = trigDecision;
			treeVar_nVtx_  = evtHelper.GetRecoVtx_nGoodVtxs();

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
			mainAnaTree_->Fill();
		}
	};

	///////////////////////////////////////////////////////////////////////////////
	// EffiCalcTree: Simple class used to generate trees that are used in calculating
	//               the efficiency of cuts on reconstruction of Z bosons
	//
	class EffiCalcTree : public TreeHandlerBase {
	private:
		// PRIVATE MEMBERS

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
		EffiCalcTree(const std::string& fileName) :
			TreeHandlerBase("zBosonEffiTree", "Tree of Z candidate information for signal MC effi calc'ns", fileName)
		{
			// Setting up the event / di-ele branches ...
			mainAnaTree_->Branch("weight",      &treeVar_weight_, "weight/D");
			mainAnaTree_->Branch("mc_numVtx", &treeVar_mc_numVtx_, "mc_numVtx/f");
			mainAnaTree_->Branch("run", &treeVar_runNum_, "run/i");
			mainAnaTree_->Branch("lumi", &treeVar_lumiSec_, "lumi/i");
			mainAnaTree_->Branch("evtNum", &treeVar_evtNum_, "evtNum/i");

			mainAnaTree_->Branch("mcZ_ele1_p4", &treeVar_mcZ_ele1_p4Ptr_);	treeVar_mcZ_ele1_p4Ptr_ = &treeVar_mcZ_ele1_p4_;
			mainAnaTree_->Branch("mcZ_ele2_p4", &treeVar_mcZ_ele2_p4Ptr_);	treeVar_mcZ_ele2_p4Ptr_ = &treeVar_mcZ_ele2_p4_;
			mainAnaTree_->Branch("mcAccept_pt",   &treeVar_ptAcc_, "mcAccept_pt/O");
			mainAnaTree_->Branch("mcAccept_ebeb", &treeVar_ebebAcceptance_, "mcAccept_ebeb/O");
			mainAnaTree_->Branch("mcAccept_ebee", &treeVar_ebeeAcceptance_, "mcAccept_ebee/O");
			mainAnaTree_->Branch("mcAccept_eeee", &treeVar_eeeeAcceptance_, "mcAccept_eeee/O");

			mainAnaTree_->Branch("bothRecod", &treeVar_bothRecod_, "bothRecod/O");

			mainAnaTree_->Branch("ZpT", &treeVar_ZpT_,     "ZpT/D");
			mainAnaTree_->Branch("ZdEta", &treeVar_ZdEta_,     "ZdEta/D");
			mainAnaTree_->Branch("ZdPhi", &treeVar_ZdPhi_,     "ZdPhi/D");
			mainAnaTree_->Branch("ZdR", &treeVar_ZdR_, "ZdR/D");
			mainAnaTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);	treeVar_eleA_p4Ptr_ = &treeVar_eleA_p4_;
			mainAnaTree_->Branch("eleA_dRmc", &treeVar_eleA_dRmc_, "eleA_dRmc/D");
			mainAnaTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);	treeVar_eleB_p4Ptr_ = &treeVar_eleB_p4_;
			mainAnaTree_->Branch("eleB_dRmc", &treeVar_eleB_dRmc_, "eleB_dRmc/D");


			mainAnaTree_->Branch("cut_both_fiducial", &treeVar_cut_both_fiducial_, "cut_both_fiducial/O");
			mainAnaTree_->Branch("cut_both_ecalDriven", &treeVar_cut_both_ecalDriven_, "cut_both_ecalDriven/O");
			mainAnaTree_->Branch("cut_both_dEta", &treeVar_cut_both_dEta_, "cut_both_dEta/O");
			mainAnaTree_->Branch("cut_both_dPhi", &treeVar_cut_both_dPhi_, "cut_both_dPhi/O");
			mainAnaTree_->Branch("cut_both_hOverE", &treeVar_cut_both_hOverE_, "cut_both_hOverE/O");
			mainAnaTree_->Branch("cut_both_showerShape", &treeVar_cut_both_showerShape_, "cut_both_showerShape/O");
			mainAnaTree_->Branch("cut_both_heepId", &treeVar_cut_both_heepId_, "cut_both_heepId/O");

			mainAnaTree_->Branch("cut_eleA_stdTrkIso", &treeVar_cut_eleA_stdTrkIso_, "cut_eleA_stdTrkIso/O");
			mainAnaTree_->Branch("cut_eleB_stdTrkIso", &treeVar_cut_eleB_stdTrkIso_, "cut_eleB_stdTrkIso/O");
			mainAnaTree_->Branch("eleA_stdTrkIso", &treeVar_eleA_stdTrkIso_, "eleA_stdTrkIso/D");
			mainAnaTree_->Branch("eleB_stdTrkIso", &treeVar_eleB_stdTrkIso_, "eleB_stdTrkIso/D");
			mainAnaTree_->Branch("cut_eleA_stdEmH1Iso", &treeVar_cut_eleA_stdEmH1Iso_, "cut_eleA_stdEmH1Iso/O");
			mainAnaTree_->Branch("cut_eleB_stdEmH1Iso", &treeVar_cut_eleB_stdEmH1Iso_, "cut_eleB_stdEmH1Iso/O");
			mainAnaTree_->Branch("eleA_stdEmH1Iso", &treeVar_eleA_stdEmH1Iso_, "eleA_stdEmH1Iso/D");
			mainAnaTree_->Branch("eleB_stdEmH1Iso", &treeVar_eleB_stdEmH1Iso_, "eleB_stdEmH1Iso/D");

			mainAnaTree_->Branch("eleA_nTrksInnerVeto", &treeVar_eleA_nTrksInnerVeto_, "eleA_nTrksInnerVeto/i");
			mainAnaTree_->Branch("eleB_nTrksInnerVeto", &treeVar_eleB_nTrksInnerVeto_, "eleB_nTrksInnerVeto/i");

			mainAnaTree_->Branch("cut_eleA_modTrkIso", &treeVar_cut_eleA_modTrkIso_, "cut_eleA_modTrkIso/O");
			mainAnaTree_->Branch("cut_eleB_modTrkIso", &treeVar_cut_eleB_modTrkIso_, "cut_eleB_modTrkIso/O");
			mainAnaTree_->Branch("eleA_modTrkIso", &treeVar_eleA_modTrkIso_, "eleA_modTrkIso/D");
			mainAnaTree_->Branch("eleB_modTrkIso", &treeVar_eleB_modTrkIso_, "eleB_modTrkIso/D");
			mainAnaTree_->Branch("cut_eleA_scModEmH1Iso", &treeVar_cut_eleA_scModEmH1Iso_, "cut_eleA_scModEmH1Iso/O");
			mainAnaTree_->Branch("cut_eleB_scModEmH1Iso", &treeVar_cut_eleB_scModEmH1Iso_, "cut_eleB_scModEmH1Iso/O");
			mainAnaTree_->Branch("eleA_scModEmH1Iso", &treeVar_eleA_scModEmH1Iso_, "eleA_scModEmH1Iso/D");
			mainAnaTree_->Branch("eleB_scModEmH1Iso", &treeVar_eleB_scModEmH1Iso_, "eleB_scModEmH1Iso/D");

			mainAnaTree_->Branch("combThr_EmH1", &treeVar_combThr_EmH1_, "combThr_EmH1/D");

			mainAnaTree_->Branch("eleA_isoDep_stdTrk",      &treeVar_eleA_isoDep_stdTrk_,  "eleA_isoDep_stdTrk/D");
			mainAnaTree_->Branch("eleB_isoDep_stdTrk",      &treeVar_eleB_isoDep_stdTrk_, "eleB_isoDep_stdTrk/D");
			mainAnaTree_->Branch("eleA_isoDep_stdEmH1",     &treeVar_eleA_isoDep_stdEmH1_, "eleA_isoDep_stdEmH1/D");
			mainAnaTree_->Branch("eleB_isoDep_stdEmH1",     &treeVar_eleB_isoDep_stdEmH1_, "eleB_isoDep_stdEmH1/D");
			mainAnaTree_->Branch("cut_eleA_isoDep_stdEmH1", &treeVar_cut_eleA_isoDep_stdEmH1_, "cut_eleA_isoDep_stdEmH1/D");
			mainAnaTree_->Branch("cut_eleB_isoDep_stdEmH1", &treeVar_cut_eleB_isoDep_stdEmH1_, "cut_eleB_isoDep_stdEmH1/D");
			mainAnaTree_->Branch("eleA_isoDep_inrVetoModTrk",     &treeVar_eleA_isoDep_inrVetoModTrk_, "eleA_isoDep_inrVetoModTrk/D");
			mainAnaTree_->Branch("eleB_isoDep_inrVetoModTrk",     &treeVar_eleB_isoDep_inrVetoModTrk_, "eleB_isoDep_inrVetoModTrk/D");
			mainAnaTree_->Branch("eleA_isoDep_inrVetoModEmH1",    &treeVar_eleA_isoDep_inrVetoModEmH1_, "eleA_isoDep_inrVetoModEmH1/D");
			mainAnaTree_->Branch("eleB_isoDep_inrVetoModEmH1",    &treeVar_eleB_isoDep_inrVetoModEmH1_, "eleB_isoDep_inrVetoModEmH1/D");
			mainAnaTree_->Branch("cut_eleA_isoDep_inrVetoModEmH1", &treeVar_cut_eleA_isoDep_inrVetoModEmH1_, "cut_eleA_isoDep_inrVetoModEmH1/O");
			mainAnaTree_->Branch("cut_eleB_isoDep_inrVetoModEmH1", &treeVar_cut_eleB_isoDep_inrVetoModEmH1_, "cut_eleB_isoDep_inrVetoModEmH1/O");

			mainAnaTree_->Branch("eleA_inrXSVetoModTrk", &treeVar_eleA_inrXSVetoModTrk_, "eleA_inrXSVetoModTrk/D");
			mainAnaTree_->Branch("eleA_inrSVetoModTrk",  &treeVar_eleA_inrSVetoModTrk_,  "eleA_inrSVetoModTrk/D");
			mainAnaTree_->Branch("eleA_inrMVetoModTrk",  &treeVar_eleA_inrMVetoModTrk_,  "eleA_inrMVetoModTrk/D");
			mainAnaTree_->Branch("eleA_inrLVetoModTrk",  &treeVar_eleA_inrLVetoModTrk_,  "eleA_inrLVetoModTrk/D");
			mainAnaTree_->Branch("eleA_inrXLVetoModTrk", &treeVar_eleA_inrXLVetoModTrk_, "eleA_inrXLVetoModTrk/D");

			mainAnaTree_->Branch("eleB_inrXSVetoModTrk", &treeVar_eleB_inrXSVetoModTrk_, "eleB_inrXSVetoModTrk/D");
			mainAnaTree_->Branch("eleB_inrSVetoModTrk",  &treeVar_eleB_inrSVetoModTrk_,  "eleB_inrSVetoModTrk/D");
			mainAnaTree_->Branch("eleB_inrMVetoModTrk",  &treeVar_eleB_inrMVetoModTrk_,  "eleB_inrMVetoModTrk/D");
			mainAnaTree_->Branch("eleB_inrLVetoModTrk",  &treeVar_eleB_inrLVetoModTrk_,  "eleB_inrLVetoModTrk/D");
			mainAnaTree_->Branch("eleB_inrXLVetoModTrk", &treeVar_eleB_inrXLVetoModTrk_, "eleB_inrXLVetoModTrk/D");

			mainAnaTree_->Branch("eleA_inrXSVetoModEmH1", &treeVar_eleA_inrXSVetoModEmH1_, "eleA_inrXSVetoModEmH1/D");
			mainAnaTree_->Branch("eleA_inrSVetoModEmH1",  &treeVar_eleA_inrSVetoModEmH1_,  "eleA_inrSVetoModEmH1/D");
			mainAnaTree_->Branch("eleA_inrMVetoModEmH1",  &treeVar_eleA_inrMVetoModEmH1_,  "eleA_inrMVetoModEmH1/D");
			mainAnaTree_->Branch("eleA_inrLVetoModEmH1",  &treeVar_eleA_inrLVetoModEmH1_,  "eleA_inrLVetoModEmH1/D");
			mainAnaTree_->Branch("eleA_inrXLVetoModEmH1", &treeVar_eleA_inrXLVetoModEmH1_, "eleA_inrXLVetoModEmH1/D");

			mainAnaTree_->Branch("eleB_inrXSVetoModEmH1", &treeVar_eleB_inrXSVetoModEmH1_, "eleB_inrXSVetoModEmH1/D");
			mainAnaTree_->Branch("eleB_inrSVetoModEmH1",  &treeVar_eleB_inrSVetoModEmH1_,  "eleB_inrSVetoModEmH1/D");
			mainAnaTree_->Branch("eleB_inrMVetoModEmH1",  &treeVar_eleB_inrMVetoModEmH1_,  "eleB_inrMVetoModEmH1/D");
			mainAnaTree_->Branch("eleB_inrLVetoModEmH1",  &treeVar_eleB_inrLVetoModEmH1_,  "eleB_inrLVetoModEmH1/D");
			mainAnaTree_->Branch("eleB_inrXLVetoModEmH1", &treeVar_eleB_inrXLVetoModEmH1_, "eleB_inrXLVetoModEmH1/D");

			mainAnaTree_->Branch("cut_eleA_inrMVetoModEmH1", &treeVar_cut_eleA_inrMVetoModEmH1_, "cut_eleA_inrMVetoModEmH1/O");
			mainAnaTree_->Branch("cut_eleB_inrMVetoModEmH1", &treeVar_cut_eleB_inrMVetoModEmH1_, "cut_eleB_inrMVetoModEmH1/O");

			mainAnaTree_->Branch("eleA_nGenHadronsDr04", &treeVar_eleA_nGenHadronsDr04_, "eleA_nGenHadronsDr04/i");
			mainAnaTree_->Branch("eleB_nGenHadronsDr04", &treeVar_eleB_nGenHadronsDr04_, "eleB_nGenHadronsDr04/i");
			mainAnaTree_->Branch("eleA_ptSumGenHadronsDr04", &treeVar_eleA_ptSumGenHadronsDr04_, "eleA_ptSumGenHadronsDr04/D");
			mainAnaTree_->Branch("eleB_ptSumGenHadronsDr04", &treeVar_eleB_ptSumGenHadronsDr04_, "eleB_ptSumGenHadronsDr04/D");

			mainAnaTree_->Branch("eleA_EmH1RhoCorrn", &treeVar_eleA_EmH1RhoCorrn_, "eleA_EmH1RhoCorrn_/D");
			mainAnaTree_->Branch("eleB_EmH1RhoCorrn", &treeVar_eleB_EmH1RhoCorrn_, "eleB_EmH1RhoCorrn_/D");

			mainAnaTree_->Branch("eleA_heepIdModIsoCutCode", &treeVar_eleA_heepIdModIsoCutCode_, "eleA_heepIdModIsoCutCode/I");
			mainAnaTree_->Branch("eleB_heepIdModIsoCutCode", &treeVar_eleB_heepIdModIsoCutCode_, "eleB_heepIdModIsoCutCode/I");
		}
		~EffiCalcTree(){ }

		void fillTree(tsw::HEEPDiEle* diEle, const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2,
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
			mainAnaTree_->Fill();
		}

		void fillTree_NotReconstructed(const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2,
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
			mainAnaTree_->Fill();
		}
	};

}

//-----------------------------------------------------------------------------------------------------
//=========================================== Analyser class ==========================================

class BstdZeeFirstAnalyser{
	public:
		BstdZeeFirstAnalyser(int runMode, int numEvts, bool isMC, const std::vector<std::string>& inFileNamesVec, const std::string& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax);
		~BstdZeeFirstAnalyser();
		void DoAnalysis(const Double_t evtWeight);
		
	private:
		void AnalyseEvent(const Double_t );
		void SetupEleClassVectors();
		void SetupMuonCollection();
		void PrintOutBranchVariables();
		void FinishOffAnalysis(); //Method to be called after all events analysed - i.e. for normalisation of histograms, calculating errors on each bins, etc

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
		void SetupBranchLinks();
		void SetupEMuMethodTrees();

	
	//Member variables...
	private:
		const int runMode_;
		const int vFlg_;
		const int numEvtsRequested_;
		unsigned int numEvtsRunOver_;
		const bool isMC_;

		TStopwatch timer_DoAnalysis_readIn_;
		TStopwatch timer_DoAnalysis_setup_;
		TStopwatch timer_DoAnalysis_FillHistograms_;
		TStopwatch timer_DoAnalysis_DoEMuMethod_;


		//Input and output files...
		const std::vector<std::string> inputFileNamesVec_;
		TChain* inputFilesTChain_;
		const std::string outputFileName_;

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

BstdZeeFirstAnalyser::BstdZeeFirstAnalyser(int runMode, int numEvts, bool isMC, const std::vector<std::string>& inFileNamesVec, const std::string& outFileName, int vFlg, Int_t nbins_mass, Int_t nbins_pt, Double_t ptmax) :
   // Initialise the member variables that are arguments of the CTOR...
   runMode_(runMode),
   vFlg_(vFlg),
	numEvtsRequested_(numEvts),
   isMC_(isMC),
   inputFileNamesVec_( tsw::AddPrefixesToFileNamesVec(inFileNamesVec) ),
   inputFilesTChain_( new TChain("demo/EventDataTree") ),
   outputFileName_( outFileName.rfind(".root")==(outFileName.length()-5) ? outFileName.substr(0, outFileName.length()-5) : outFileName ),
   //
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
	noIsoZCandDiEleTree_("noIsoZBosonTree", outputFileName_ + "_noIsoZCandTree.root"),
	modIsoZCandDiEleTree_("modIsoZBosonTree", outputFileName_ + "_modIsoZCandTree.root"),
	zCandEffiTree_(outputFileName_ + "_zEffiTree.root"),
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

	timer_DoAnalysis_readIn_.Stop(); timer_DoAnalysis_readIn_.Reset();
	timer_DoAnalysis_setup_.Stop(); timer_DoAnalysis_setup_.Reset();
	timer_DoAnalysis_FillHistograms_.Stop(); timer_DoAnalysis_FillHistograms_.Reset();
	timer_DoAnalysis_DoEMuMethod_.Stop(); timer_DoAnalysis_DoEMuMethod_.Reset();
	
	//Set default (typically clearly 'incorrect') values for non-histo variables...
	SetMemberVariablesToDefaultValues();

	//Setting up the emu method trees ...
	SetupEMuMethodTrees();

	// ----- SETTING UP TCHAIN AND ASSOCIATED BRANCH LINKS ----- //
	unsigned int inputTrees_numEvts = 0;
   std::vector<std::string>::const_iterator fileNameIt;

   TStopwatch nEntriesTimer; nEntriesTimer.Start();
   for( fileNameIt = inputFileNamesVec_.begin(); fileNameIt != inputFileNamesVec_.end(); fileNameIt++ ){
      inputFilesTChain_->Add(fileNameIt->c_str());
      TFile* theFile = TFile::Open(fileNameIt->c_str(), "READ");
      if( TTree* treeThisFile = dynamic_cast<TTree*>(theFile->Get("demo/EventDataTree")) ){
      	unsigned int nEvtsThisTree = treeThisFile->GetEntries();
      	std::cout << nEvtsThisTree << ", ";
      	inputTrees_numEvts += nEvtsThisTree;
      }
      else
      	std::cerr << " ERROR : Could not get TTree from file '" << (*fileNameIt) << "'" << std::endl;
      theFile->Close();
   }
   nEntriesTimer.Stop(); nEntriesTimer.Print();

   SetupBranchLinks();
	inputFilesTChain_->SetCacheSize(10000000);
	inputFilesTChain_->AddBranchToCache("*");

	// Determine the number of events being run over ...
	if( numEvtsRequested_<0 || (numEvtsRequested_>static_cast<int>(inputTrees_numEvts)) )
		numEvtsRunOver_ = inputTrees_numEvts;
	else
		numEvtsRunOver_ = numEvtsRequested_;
}

BstdZeeFirstAnalyser::~BstdZeeFirstAnalyser(){
	//Clear any heap variables here...
	delete inputFilesTChain_;
}

void BstdZeeFirstAnalyser::DoAnalysis(const Double_t weightFromXsec)
{
	std::cout << std::endl << "Now, running BstdZeeAnalyser with the following parameters ... " << std::endl
			    << "   - input files: " << std::endl;
	for(unsigned int idx=0; idx<inputFileNamesVec_.size(); idx++)
		std::cout << "        " << inputFileNamesVec_.at(idx) << std::endl;
	std::cout << "   - output prefix: " << outputFileName_ << std::endl
				 << "   - Events are " << (isMC_ ? "Monte-Carlo" : "real DATA") << std::endl
				 << "   - Run mode : " << runMode_ << " ,  Verbosity flag = " << vFlg_ << std::endl
				 << "   - xSectionWeight = " << weightFromXsec << std::endl
				 << "     " << GetNumEvtsRunOver() << " events will be run over ..." << std::endl;


	//Call AnalyseEvent method for each event...
	TStopwatch timer_AnalysingEvents; timer_AnalysingEvents.Start();
	for(unsigned int evtIdx=0; evtIdx<numEvtsRunOver_; evtIdx++){
		if(vFlg_>0){std::cout << std::endl << "        + Loading event no. " << evtIdx << std::endl;}
		//Load in data for the evtIdx'th event...
		timer_DoAnalysis_readIn_.Start(false);
		inputFilesTChain_->GetEntry(evtIdx);
		timer_DoAnalysis_readIn_.Stop();

		if( (vFlg_>-2) && (evtIdx%500000==0 || evtIdx==(numEvtsRunOver_-1)) ){std::cout << "        * Analysing event no. " << evtIdx << std::endl;}

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

	// Now finish off the analysis (incl. opening output files) ...
	TStopwatch timer_FinishOffAnalysis; timer_FinishOffAnalysis.Start();
	FinishOffAnalysis();
	timer_FinishOffAnalysis.Stop();

	// Print out timing info
	std::cout << std::endl
				 << " -=-=- DETAILED TIMING INFO -=-=-" << std::endl
				 << "   1. During event analysis:" << std::endl
				 << "       + Total: " << std::endl
				 << "            " << tsw::TimeInfoString(timer_AnalysingEvents) << std::endl
				 << "       - ReadIn:" << std::endl
				 << "            " << tsw::TimeInfoString(timer_DoAnalysis_readIn_) << std::endl
				 << "       - In Setup:" << std::endl
				 << "            " << tsw::TimeInfoString(timer_DoAnalysis_setup_) << std::endl
				 << "       - In FillHistograms:" << std::endl
				 << "            " << tsw::TimeInfoString(timer_DoAnalysis_FillHistograms_) << std::endl
				 << "       - In DoEMuMethod:" << std::endl
				 << "            " << tsw::TimeInfoString(timer_DoAnalysis_DoEMuMethod_) << std::endl
				 << "   2. During FinishOffAnalysis:" << std::endl
				 << "       + Overall: " << std::endl
				 << "            " << tsw::TimeInfoString(timer_FinishOffAnalysis) << std::endl << std::endl;

//	timer_AnalysingEvents.Print();
//	std::cout << std::endl << "*** For timer_DoAnalysis_readIn_:" << std::endl;
//	timer_DoAnalysis_readIn_.Print();
//	std::cout << "*** For timer_DoAnalysis_setup_:" << std::endl;
//	timer_DoAnalysis_setup_.Print();
//	std::cout << "*** For timer_DoAnalysis_FillHistograms_:" << std::endl;
//	timer_DoAnalysis_FillHistograms_.Print();
//	std::cout << "*** For timer_DoAnalysis_DoEMuMethod_:" << std::endl;
//	timer_DoAnalysis_DoEMuMethod_.Print();
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
	if(vFlg_>0){std::cout << " ... and setting up the new Ele class vectors for this event..." << std::endl;}

	for(unsigned int iEle = 0; iEle < normGsfEles_number_; iEle++){
		if(vFlg_>0){std::cout << "   iEle=" << iEle << std::endl;}
		SetDefaultValuesInEleStruct(ithEleStruct);
		ithHEEPEle = tsw::HEEPEle(ithEleStruct);

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

		ithHEEPEle = tsw::HEEPEle(ithEleStruct);
		normEles_.push_back(ithHEEPEle);
	}

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

void BstdZeeFirstAnalyser::PrintOutBranchVariables()
{
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
		std::cout << "     ->HEEP variables ..." << std::endl;
		normEles_.at(iEle).PrintOutVariables();
	}

}

void BstdZeeFirstAnalyser::FinishOffAnalysis(){

	//Calculate errors on histogram bins ...



	//Open the output file ...
	TFile f_histos((outputFileName_+"_histos.root").c_str(),"RECREATE");
	f_histos.Write();

	//Save all of the histograms from analysis ...
	SaveReconValidationHistos(&f_histos);
	//Close the output file ...
	f_histos.Close();

	// Open up eMu method files, save the appropriate tree, and then
	TFile f_eleMuTree((outputFileName_ + "_eMuTree.root").c_str(),"RECREATE");
	eleMuTree_->Write();
	f_eleMuTree.Close();
	delete eleMuTree_;

	TFile f_diEleTree((outputFileName_ + "_diEleTree.root").c_str(),"RECREATE");
	diEleTree_->Write();
	f_diEleTree.Close();
	delete diEleTree_;

	// Save the ABCD QCD estimation tree ...
//	frPreDiEleTree_.SaveToFile("/opt/ppd/newscratch/williams/Datafiles/abcdDiEleTrees/" + outputFile_name_ + "_abcdTree.root");

	// Set the event counters for the tree handlers, and write trees to file ...
	noIsoZCandDiEleTree_.setEventCounter( this->GetNumEvtsRunOver() );
	noIsoZCandDiEleTree_.saveToFile();
	modIsoZCandDiEleTree_.setEventCounter( this->GetNumEvtsRunOver() );
	modIsoZCandDiEleTree_.saveToFile();
	zCandEffiTree_.setEventCounter( this->GetNumEvtsRunOver() );
	zCandEffiTree_.saveToFile();

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
		normDiEle_AllHEEP = tsw::HEEPDiEle( normEles_, normEles_HEEPCutsFlags );
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
		normDiEle_HEEPNoIso = tsw::HEEPDiEle( normEles_, normEles_HEEPCutsNoIsoFlags );

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
		normDiEle_EB_HEEPNoIso = tsw::HEEPDiEle( normEles_, normEles_EB_HEEPCutsNoIsoFlags);
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

			noIsoZCandDiEleTree_.fillTree(normDiEle_EB_HEEPNoIso, eventHelper_, trg_PathA_decision_);
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
		zCandEffiTree_.fillTree(mcMatchedDiEle, mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_);
		delete mcMatchedDiEle;
	}
	else
		zCandEffiTree_.fillTree_NotReconstructed(mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_);

	// ------------------------
	// EB HEEPNoIso di-electrons...
	if(tsw::NumPassingCuts(normEles_EB_heepIdModIsoFlags)>1){
		tsw::HEEPDiEle heepIdModIsoEbDiEle( normEles_, normEles_EB_heepIdModIsoFlags);

		if( heepIdModIsoEbDiEle.isInZMassRange() )
			modIsoZCandDiEleTree_.fillTree(heepIdModIsoEbDiEle, eventHelper_, trg_PathA_decision_);

	}


	//-----------------------------------
	// ABCD method selection code:
	// EB fid, ECAL-driven, fake rate (FR) pre-selection di-electrons
	if( tsw::NumPassingCuts(normEles_EB_isFidEcalDrAndFRPreFlags)>1 ){
		// Form the di-electron
		normDiEle_EB_fidECALDrFRPre = tsw::HEEPDiEle( normEles_, normEles_EB_isFidEcalDrAndFRPreFlags);
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

}

//================================================================//
void BstdZeeFirstAnalyser::SetupBranchLinks()
{
	inputFilesTChain_->SetBranchAddress("event", &event_); // inputFile_tree_->SetBranchStatus("event",0); //
	//Setting up the pointer links for the event information branches ...
	inputFilesTChain_->SetBranchAddress("evt_runNum",&evt_runNum_);   // unsigned int
	inputFilesTChain_->SetBranchAddress("evt_lumiSec",&evt_lumiSec_); // unsigned int
	inputFilesTChain_->SetBranchAddress("evt_evtNum",&evt_evtNum_);   // unsigned int

	//Setting up the pointer links for the trigger branches ...
	inputFilesTChain_->SetBranchAddress("trg_PathA_decision", &trg_PathA_decision_);
	inputFilesTChain_->SetBranchAddress("trg_PathA_name",     &trg_PathA_name_ptr_);
	inputFilesTChain_->SetBranchAddress("trg_PathA_nameOfLastFilter", &trg_PathA_nameOfLastFilter_ptr_);
	inputFilesTChain_->SetBranchAddress("trg_PathA_highestTrigObjEt", &trg_PathA_highestTrigObjEt_);


	if(isMC_){
		//Setting up the pointer links for the MC electron branches ...
		inputFilesTChain_->SetBranchAddress("mcEles_number",          &mc_numFinalStateEles_      );
		inputFilesTChain_->SetBranchAddress("mcEles_HighestEt_charge",&mcEles_HighestEt_charge_   );
		inputFilesTChain_->SetBranchAddress("mcEles_HighestEt_PDGid", &mcEles_HighestEt_PDGid_    );
		inputFilesTChain_->SetBranchAddress("mcEles_HighestEt_status",&mcEles_HighestEt_status_   );
		inputFilesTChain_->SetBranchAddress("mcEles_HighestEt_pt",    &mcEles_HighestEt_pt_       );
		inputFilesTChain_->SetBranchAddress("mcEles_HighestEt_p4",    &mcEles_HighestEt_OLDp4_ptr_);

		inputFilesTChain_->SetBranchAddress("mcEles_2ndHighestEt_charge", &mcEles_2ndHighestEt_charge_   );
		inputFilesTChain_->SetBranchAddress("mcEles_2ndHighestEt_PDGid",  &mcEles_2ndHighestEt_PDGid_    );
		inputFilesTChain_->SetBranchAddress("mcEles_2ndHighestEt_status", &mcEles_2ndHighestEt_status_   );
		inputFilesTChain_->SetBranchAddress("mcEles_2ndHighestEt_pt",     &mcEles_2ndHighestEt_pt_       );
		inputFilesTChain_->SetBranchAddress("mcEles_2ndHighestEt_p4",     &mcEles_2ndHighestEt_OLDp4_ptr_);

		inputFilesTChain_->SetBranchAddress("mcZcandidate_pt",      &mcZcandidate_pt_  ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_eta",     &mcZcandidate_eta_ ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_phi",     &mcZcandidate_phi_ ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_mass",    &mcZcandidate_mass_); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_p4",           &mcZcandidate_OLDp4_ptr_   ); //ROOT::Math::XYZTLorentzVector
		inputFilesTChain_->SetBranchAddress("mcZcandidate_dEtaEles",     &mcZcandidate_dEtaEles_    ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_dPhiEles",     &mcZcandidate_dPhiEles_    ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_dREles",       &mcZcandidate_dREles_      ); //Double_t
		inputFilesTChain_->SetBranchAddress("mcZcandidate_openingAngle", &mcZcandidate_openingAngle_); //Double_t

		inputFilesTChain_->SetBranchAddress("mcZboson_pdgId",  &mcZ_pdgId_ );
		inputFilesTChain_->SetBranchAddress("mcZboson_status", &mcZ_status_);
		inputFilesTChain_->SetBranchAddress("mcZboson_p4",     &mcZ_p4_ptr_); //ROOT::Math::XYZTVector //NB: These pointers have to be initialised!! (Initialisation to 0 is fine...)
		inputFilesTChain_->SetBranchAddress("mcZboson_numDaughters",   &mcZ_numDaughters_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughters_dR",   &mcZ_daughters_dR_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughters_dEta", &mcZ_daughters_dEta_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughters_dPhi", &mcZ_daughters_dPhi_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughters_openingAngle", &mcZ_daughters_openingAngle_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughterA_p4", &mcZ_daughterA_OLDp4_ptr_);
		inputFilesTChain_->SetBranchAddress("mcZboson_daughterB_p4", &mcZ_daughterB_OLDp4_ptr_);
	}

	inputFilesTChain_->SetBranchStatus("normGsfEles_*",0);
	inputFilesTChain_->SetBranchStatus("normGsfEles_number",1);
	inputFilesTChain_->SetBranchStatus("normGsfEles_charge",1);
	//Setting up the pointer links for the 'normal' GSF electron branches ...
	inputFilesTChain_->SetBranchAddress("normGsfEles_number", &normGsfEles_number_); //unsigned int
	//inputFilesTChain_->SetBranchAddress("normGsfEles_p4ptr_", &normGsfEles_OLDp4_  ); // std::vector<ROOT::Math::XYZTVector>*
	inputFilesTChain_->SetBranchAddress("normGsfEles_charge", &normGsfEles_charge_); //std::vector<Int_t>*
	if(false){
		inputFilesTChain_->SetBranchAddress("normGsfEles_Et",      &normGsfEles_Et_     );           // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_HEEP_Et", &normGsfEles_HEEP_Et_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_Eta",     &normGsfEles_Eta_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_scEta",   &normGsfEles_scEta_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_ecalDriven",     &normGsfEles_ecalDriven_);         // std::vector<bool>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_ecalDrivenSeed", &normGsfEles_ecalDrivenSeed_); // std::vector<bool>*

		inputFilesTChain_->SetBranchAddress("normGsfEles_HEEP_dEtaIn", & normGsfEles_HEEP_dEtaIn_);    // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_HEEP_dPhiIn", &normGsfEles_HEEP_dPhiIn_);     // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_HoverE", &normGsfEles_HoverE_);               // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_sigmaIetaIeta", &normGsfEles_sigmaIetaIeta_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_scSigmaIetaIeta", &normGsfEles_scSigmaIetaIeta_); // std::vector<Double_t>*

		inputFilesTChain_->SetBranchAddress("normGsfEles_dr03EmIsoEt", &normGsfEles_dr03EmIsoEt_);             // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_dr03HadDepth1IsoEt", &normGsfEles_dr03HadDepth1IsoEt_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_dr03HadDepth2IsoEt", &normGsfEles_dr03HadDepth2IsoEt_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_dr03TkIsoPt", &normGsfEles_dr03TkIsoPt_);  // std::vector<Double_t>*

		inputFilesTChain_->SetBranchAddress("normGsfEles_e2x5Max", &normGsfEles_e2x5Max_); // std::vector<Double_t>*
		inputFilesTChain_->SetBranchAddress("normGsfEles_e5x5", &normGsfEles_e5x5_); // std::vector<Double_t>*
	}
	//*** Branches containing heep::Ele method values for norm eles ***
	//kinematic and geometric methods
	inputFilesTChain_->SetBranchAddress("normHEEPEles_et", &normHEEPEles_et_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfEt", &normHEEPEles_gsfEt_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_scEt", &normHEEPEles_scEt_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_energy", &normHEEPEles_energy_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_gsfEnergy",0);									//	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfEnergy", &normHEEPEles_gsfEnergy_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_caloEnergy", &normHEEPEles_caloEnergy_);
//	inputFilesTChain_->SetBranchAddress("normHEEPEles_ecalEnergyError", &normHEEPEles_ecalEnergyError_);			/* TEMP v1f/g FIX */
	inputFilesTChain_->SetBranchAddress("normHEEPEles_eta", &normHEEPEles_eta_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_scEta", &normHEEPEles_scEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_detEta",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_detEta", &normHEEPEles_detEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_detEtaAbs",0);									//	inputFilesTChain_->SetBranchAddress("normHEEPEles_detEtaAbs", &normHEEPEles_detEtaAbs_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_phi", &normHEEPEles_phi_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_scPhi", &normHEEPEles_scPhi_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_detPhi",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_detPhi", &normHEEPEles_detPhi_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_zVtx",0);											//	inputFilesTChain_->SetBranchAddress("normHEEPEles_zVtx", &normHEEPEles_zVtx_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_p4ptr", &normHEEPEles_p4ptr_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfP4ptr", &normHEEPEles_gsfP4ptr_);

	//classification (couldnt they have just named it 'type')
	inputFilesTChain_->SetBranchStatus("normHEEPEles_classification",0);							//	inputFilesTChain_->SetBranchAddress("normHEEPEles_classification", &normHEEPEles_classification_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isEcalDriven", &normHEEPEles_isEcalDriven_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_isTrackerDriven",0);							//	inputFilesTChain_->SetBranchAddress("normHEEPEles_isTrackerDriven", &normHEEPEles_isTrackerDriven_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isEB", &normHEEPEles_isEB_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isEE", &normHEEPEles_isEE_);

	//track methods
	inputFilesTChain_->SetBranchAddress("normHEEPEles_charge", &normHEEPEles_charge_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_trkCharge",0);									//	inputFilesTChain_->SetBranchAddress("normHEEPEles_trkCharge", &normHEEPEles_trkCharge_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_pVtx",0);											//	inputFilesTChain_->SetBranchAddress("normHEEPEles_pVtx", &normHEEPEles_pVtx_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_pCalo",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_pCalo", &normHEEPEles_pCalo_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_ptVtx", &normHEEPEles_ptVtx_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_ptCalo", &normHEEPEles_ptCalo_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_pt",  &normHEEPEles_closestCtfTrk_pt_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_eta", &normHEEPEles_closestCtfTrk_eta_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_phi", &normHEEPEles_closestCtfTrk_phi_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerPt",0);  // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerPt",  &normHEEPEles_closestCtfTrk_innerPt_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerEta",0); // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerEta", &normHEEPEles_closestCtfTrk_innerEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_innerPhi",0); // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_innerPhi", &normHEEPEles_closestCtfTrk_innerPhi_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerPt",0);  // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerPt",  &normHEEPEles_closestCtfTrk_outerPt_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerEta",0); // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerEta", &normHEEPEles_closestCtfTrk_outerEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_closestCtfTrk_outerPhi",0); // inputFilesTChain_->SetBranchAddress("normHEEPEles_closestCtfTrk_outerPhi", &normHEEPEles_closestCtfTrk_outerPhi_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	inputFilesTChain_->SetBranchAddress("normHEEPEles_hOverE", &normHEEPEles_hOverE_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_dEtaIn", &normHEEPEles_dEtaIn_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_dPhiIn", &normHEEPEles_dPhiIn_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_dPhiOut",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_dPhiOut", &normHEEPEles_dPhiOut_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_epIn", &normHEEPEles_epIn_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_epOut", &normHEEPEles_epOut_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_fbrem",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_fbrem", &normHEEPEles_fbrem_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_bremFrac",0);									//	inputFilesTChain_->SetBranchAddress("normHEEPEles_bremFrac", &normHEEPEles_bremFrac_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_invEOverInvP", &normHEEPEles_invEOverInvP_);

	//shower shape variables
	inputFilesTChain_->SetBranchStatus("normHEEPEles_sigmaEtaEta",0);								//	inputFilesTChain_->SetBranchAddress("normHEEPEles_sigmaEtaEta", &normHEEPEles_sigmaEtaEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_sigmaEtaEtaUnCorr",0);						//	inputFilesTChain_->SetBranchAddress("normHEEPEles_sigmaEtaEtaUnCorr", &normHEEPEles_sigmaEtaEtaUnCorr_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_sigmaIEtaIEta", &normHEEPEles_sigmaIEtaIEta_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_e1x5",0);											//	inputFilesTChain_->SetBranchAddress("normHEEPEles_e1x5", &normHEEPEles_e1x5_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_e2x5Max",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_e2x5Max", &normHEEPEles_e2x5Max_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_e5x5",0);											//	inputFilesTChain_->SetBranchAddress("normHEEPEles_e5x5", &normHEEPEles_e5x5_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_e1x5Over5x5", &normHEEPEles_e1x5Over5x5_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_e2x5MaxOver5x5", &normHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolEm", &normHEEPEles_isolEm_);
	inputFilesTChain_->SetBranchStatus("normHEEPEles_isolHad",0);										//	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolHad", &normHEEPEles_isolHad_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolHadDepth1", &normHEEPEles_isolHadDepth1_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolHadDepth2", &normHEEPEles_isolHadDepth2_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolPtTrks", &normHEEPEles_isolPtTrks_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_isolEmHadDepth1", &normHEEPEles_isolEmHadDepth1_);

//	inputFilesTChain_->SetBranchAddress("normHEEPEles_SCposn_eta", &normHEEPEles_SCposn_eta_);								/* TEMP v1f/g FIX */
//	inputFilesTChain_->SetBranchAddress("normHEEPEles_SCposn_phi", &normHEEPEles_SCposn_phi_);								/* TEMP v1f/g FIX */
//	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_rawEnergy", &normHEEPEles_SC_rawEnergy_);							/* TEMP v1f/g FIX */
	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_recHits_Et",  &normHEEPEles_SC_recHits_Et_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_recHits_eta", &normHEEPEles_SC_recHits_eta_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_recHits_phi", &normHEEPEles_SC_recHits_phi_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_recHits_isFromEB", &normHEEPEles_SC_recHits_isFromEB_);
//	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_totEnergyRecHits", &normHEEPEles_SC_totEnergyRecHits_);		/* TEMP v1f/g FIX */
//	inputFilesTChain_->SetBranchAddress("normHEEPEles_SC_totNumRecHits", &normHEEPEles_SC_totNumRecHits_);				/* TEMP v1f/g FIX */

	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfTrk_eta", &normHEEPEles_gsfTrk_eta_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfTrk_phi", &normHEEPEles_gsfTrk_phi_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_gsfTrk_vz", &normHEEPEles_gsfTrk_vz_);

	inputFilesTChain_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_pt",  &normHEEPEles_innerIsoConeTrks_pt_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_eta", &normHEEPEles_innerIsoConeTrks_eta_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_phi", &normHEEPEles_innerIsoConeTrks_phi_);
	inputFilesTChain_->SetBranchAddress("normHEEPEles_innerIsoConeTrks_vz",  &normHEEPEles_innerIsoConeTrks_vz_);

	inputFilesTChain_->SetBranchAddress("normHEEPEles_numMissInnerHits", &normHEEPEles_numMissInnerHits_);				/* TEMP v1f/g FIX */

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


//=====================================================================================================
//-----------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	gROOT->ProcessLine("#include <vector>"); //Without this line, the std::vector<bool> branches are not linked correctly(???), giving 'dictionary for class...' errors
				//upon running of ZEventAnalyser program, and accessing these branches results in nonsensical results (e.g. methods such as size()).

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=- // 
	// PART 1: COMMAND LINE OPTIONS  //

	// Alias the boost::program_options namespace ...
	namespace bProgOpts = boost::program_options;
	
	// Declare variables that store the options
   int maxEvents = 0;
   bool isMC = false;
   double xSection = -1.0;
   double desiredIntLumi = -1.0;
   std::vector<std::string> inputFilesVecFromCmdLine;
   unsigned int numThisJob;
   unsigned int totNumJobs;
	std::string outputFile;
	
	// Declare the supported options
	bProgOpts::options_description cmdLineOptions("Allowed options");
	cmdLineOptions.add_options()
		("help",            "print help message")
		("maxEvents",       bProgOpts::value<int>(&maxEvents)->default_value(-1), "number of events to run over" )
      ("mc",              "flag for specifying that input file is MC")
      ("xSection",        bProgOpts::value<double>(&xSection)->default_value(-1.0), "cross section of sample (in fb)")
      ("scaleToLumi",     bProgOpts::value<double>(&desiredIntLumi)->default_value(-1.0), "integrated lumi for scaling events to (in inverse fb)" )
		("input-file,i",    bProgOpts::value<std::vector<std::string> >(&inputFilesVecFromCmdLine), "names of input files (Either specify as normal / positional options.)")
      ("from-list,l",     "flags that specified 'input file' is a list of files.")
      ("num-this-job",    bProgOpts::value<unsigned int>(&numThisJob)->default_value(1), "Index of this job (valid range: 1 <= jobIndex <= totNumJobs)")
      ("tot-num-jobs",    bProgOpts::value<unsigned int>(&totNumJobs)->default_value(1), "The total number of jobs being run over specified input files")
		("output-file,o",   bProgOpts::value<std::string>(&outputFile), "name of output file")
		;
   
	bProgOpts::positional_options_description cmdLinePositionalOpts;
	cmdLinePositionalOpts.add("input-file", -1);

	// Parse options given by the user
	bProgOpts::variables_map optionValueMap;
	try{
		bProgOpts::store(bProgOpts::command_line_parser(argc, argv).options(cmdLineOptions).positional(cmdLinePositionalOpts).run(),
     			optionValueMap);
      bProgOpts::notify(optionValueMap);
	}
   catch( const bProgOpts::unknown_option& e){
      std::cout << std::endl << "   * ERROR during parsing (and storing) command line options  [Exception thrown.]" 
                << std::endl << "   * The option '" << e.get_option_name() << "' has not been recognised ; program will now exit early ..." 
                << std::endl << std::endl << cmdLineOptions << std::endl;
      return 1;         
   }
   catch( const bProgOpts::multiple_occurrences& e){
      std::cout << std::endl << "   * ERROR during parsing (and storing) command line options  [Exception thrown.]" 
                << std::endl << "   * There were several occurrences of the option '" << e.get_option_name() << "', but only one occurrence is expected !"
                << std::endl << "   * The program will now exit early ..." 
                << std::endl << std::endl << cmdLineOptions << std::endl;
      return 1;         
   }
   catch( bProgOpts::invalid_option_value ){
      std::cout << std::endl << "   * ERROR during parsing (and storing) command line options  [Exception thrown.]" 
                << std::endl << "   * An invalid value was given for one of the options ; the program will now exit early ..." 
                << std::endl << std::endl << cmdLineOptions << std::endl;
      return 1;         
   }
   catch( bProgOpts::error ){
      std::cout << std::endl << "   * ERROR during parsing (and storing) command line options  [Exception thrown.]"
                << std::endl << "   * Program will now exit early ..." 
                << std::endl << std::endl << cmdLineOptions << std::endl;
      return 1;         
   }

	// Take action based on values
	if(optionValueMap.count("help")){
		std::cout << cmdLineOptions << std::endl ;//<< cmdLinePositionalOpts << std::endl;
		return 0;
	}
   
   if(optionValueMap.count("mc")!=0)
      isMC = true;
   
   if(optionValueMap.count("input-file")==0){
      std::cout << "  No input files have been entered; program is exiting early ... " << std::endl;
      return 1;
   }
   
   if(optionValueMap.count("output-file")==0){
      std::cout << "  No name has been given for the output file; program is exiting early ... " << std::endl;
      return 1;
   }
   
      
   if(optionValueMap.count("from-list")){
      std::cout << "  Detected that specified input file is a list of input files" << std::endl;
      if(inputFilesVecFromCmdLine.size()!=1){
         std::cout << "     * ERROR : There is (zero or) > one input file, and yet you used the 'from-list' flag" << std::endl
                   << "               Only one input file list can be used when running the analyser exe at the moment ..." << std::endl
                   << "               The program will now exit early!" << std::endl;
         return 1;
      }
      else{
         std::cout << "  File being read for input ROOT files is '" << inputFilesVecFromCmdLine.at(0) << "'" << std::endl; 
         ifstream theFile(inputFilesVecFromCmdLine.at(0).c_str());
         inputFilesVecFromCmdLine.clear();
         
         std::string lineContent;
         if(theFile.is_open()){
            while( ! theFile.eof() ){
               if(theFile.fail()){
                  std::cout << "     * ERROR : Failure when reading the text file listing the names of the input files. " << std::endl
                            << "               The program will now exit early!" << std::endl;
                  return 1;
               }            
               getline(theFile, lineContent);
               if(lineContent.find("#")!=0 && lineContent.length()>0)
               	inputFilesVecFromCmdLine.push_back(lineContent);
            }
         }
         else{
            std::cout << "     * ERROR : Could not open the specified file" << std::endl
                      << "               The program will now exit early!" << std::endl;
            return 1;
         }
      }// End:if-else: 1 or multiple files specified on cmd line
   }// End: if from-list option used
   
   if(numThisJob==0 || numThisJob>totNumJobs){
   	std::cerr << "  ERROR : The job no (" << numThisJob << ") specified by the user is not valid!" << std::endl
   				 << "          It should be in the range 1 to tot-num-jobs (incl.)" << std::endl
   				 << "          Program will now exit early" << std::endl;
   	return 1;
   }
   
   if(totNumJobs==0){
   	std::cerr << "  ERROR : The total number of jobs (0) specified by the user is not valid!" << std::endl
   				 << "          It must be non-zero of course !! " << std::endl
   				 << "          Program will now exit early" << std::endl;
   	return 1;
   }

   // Now work which of the input files should be run over - given the job number ...
   unsigned int numFilesFromCmdLine = inputFilesVecFromCmdLine.size();
   unsigned int minNumFilesPerJob = numFilesFromCmdLine/totNumJobs;
   unsigned int xtraNumFiles = (numFilesFromCmdLine % totNumJobs);

   unsigned int startFileIdx = 0;
   unsigned int numFilesThisJob = 0;

   for(unsigned int iJob=1; iJob<=totNumJobs; iJob++){
   	unsigned int numFilesInJobI = minNumFilesPerJob;
   	if (xtraNumFiles!=0 && iJob<=xtraNumFiles)
   		numFilesInJobI += 1;
   	if(iJob==numThisJob){
   		numFilesThisJob = numFilesInJobI; break; }
  		startFileIdx += numFilesInJobI;
   }

   std::vector<std::string> inputFilesVecToRunOver;
   for(unsigned int idx=startFileIdx; idx<(startFileIdx+numFilesThisJob) && idx<numFilesFromCmdLine; idx++)
      inputFilesVecToRunOver.push_back( inputFilesVecFromCmdLine.at(idx) );
             
   // Finally repeat the parsed options back to the user ...
   std::cout << " The specified options are ..." << std::endl 
             << "    + maxEvents     = " << maxEvents << std::endl
             << "    + mc            = " << isMC << std::endl
             << "    + xSection      = " << xSection << std::endl
             << "    + scaleToLumi   = " << desiredIntLumi << std::endl
             << "    + from-list     = " << (optionValueMap.count("from-list")!=0 ? "true" : "false") << std::endl
             << "    + num-this-job  = " << numThisJob << std::endl
             << "    + tot-num-jobs  = " << totNumJobs << std::endl
             << "    + output-file   = " << outputFile << std::endl
             << "    + input files specified (incl. those for other jobs):" << std::endl;
   for(unsigned int idx=0; idx<inputFilesVecFromCmdLine.size(); idx++)
      std::cout << "         " << (idx+1) << ")  " << inputFilesVecFromCmdLine.at(idx) << std::endl;
   std::cout << std::endl;

   // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- //
	// PART XX: ACTUALLY RUN THE ANALYSER  //

	BstdZeeFirstAnalyser theAnalyser(0, maxEvents, isMC, inputFilesVecToRunOver, outputFile, -1, 80, 55, 1100.0);

	TStopwatch watch;
	watch.Start();

	if( !isMC || xSection<=0.0 || desiredIntLumi<=0.0 )
		theAnalyser.DoAnalysis(1.0);
	else{
		double scaleFactor = (desiredIntLumi * xSection)/static_cast<double>(theAnalyser.GetNumEvtsRunOver());
		theAnalyser.DoAnalysis(scaleFactor);
	}

	watch.Stop();
	double cpuTimeInMins  = watch.CpuTime()/60.0;
	double realTimeInMins = watch.RealTime()/60.0;
	double millionEvtsRunOver = static_cast<double>( theAnalyser.GetNumEvtsRunOver() )/1000000.0;
	std::cout << " -=-=- TIMING INFO -=-=-" << std::endl
				 << "   Total: " << std::endl
				 << "      " << watch.CpuTime()/60.0 << " mins (CPU);  " << watch.RealTime()/60.0 << " mins (real)" << std::endl
				 << "   Per million events: " << std::endl
				 << "      " << cpuTimeInMins/millionEvtsRunOver << " mins (CPU);  " << realTimeInMins/millionEvtsRunOver << " mins (real)" << std::endl
				 << std::endl;

	// -=-=-=-=-=-=-=-=-=-=-=-=-=-=- //
	// *** ----- OLD CODE  ----- *** //
	
	/*//		TString inFilePrefix = "/home/ppd/nnd85574/Work/BstdZee/CMSSW_4_2_8_patch7/src/NTupler/BstdZeeNTupler/";
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
	}*/

	/*for(unsigned int i=0; i<sigInputFilenames.size()-1; i++){
		std::cout << std::endl << std::endl;
		std::cout << "  ***-------------------------------***" << std::endl;
		std::cout << "  *** Signal analysis #" << i << " ..." << std::endl;
		sigAnalysers.push_back( BstdZeeFirstAnalyser(0, 20000, false, sigInputFilenames.at(i), sigOutputFilenames.at(i), -1) );
		std::cout << " Running the DoAnalysis method ..." << std::endl;
		sigAnalysers.at(i).DoAnalysis(1.0);
	}*/

	
	return 0;
}


