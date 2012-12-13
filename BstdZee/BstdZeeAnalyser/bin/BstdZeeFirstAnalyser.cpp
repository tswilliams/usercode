
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


// BstdZeeNTupler includes
#include "TSWilliams/BstdZeeNTupler/interface/tswUsefulFunctions.h"

// BstdZeeAnalyser includes
#include "TSWilliams/BstdZeeAnalyser/interface/tswAnaFunctions.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswEventHelper.h"

#include "TSWilliams/BstdZeeAnalyser/interface/tswHEEPEle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswHEEPDiEle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswMCParticle.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswEleMuObject.h"


#ifdef __MAKECINT__
#pragma link C++ class vector<Int_t>+;
#pragma link C++ class vector< vector<float> >;
#endif


//=====================================================================================================
//-----------------------------------------------------------------------------------------------------
namespace tsw{

	std::vector<std::string> AddPrefixesToFileNamesVec(const std::vector<std::string>& origFileNamesVec){
		std::vector<std::string> newFileNamesVec;

		for(std::vector<std::string>::const_iterator fileNameIt=origFileNamesVec.begin(); fileNameIt<origFileNamesVec.end(); fileNameIt++){
			if(fileNameIt->find("/pnfs/pp.rl.ac.uk")==0)
				newFileNamesVec.push_back("dcap://dcap.pp.rl.ac.uk" + (*fileNameIt) );
			else if(fileNameIt->find("/store/")==0)
				newFileNamesVec.push_back("dcap://dcap.pp.rl.ac.uk/pnfs/pp.rl.ac.uk/data/cms" + (*fileNameIt) );
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

		std::string fileName() const { return fileName_; }

	protected:
		TTree* addTree(const std::string& treeName, const std::string& treeDescription ); ///< Returns TTree with specified name/description on heap and registers it for saving to file

		static bool isMadgraphDrellYan(const std::string& fileName) { return (fileName.find("DY")!=std::string::npos && (fileName.find("MG")!=std::string::npos || fileName.find("MADGRAPH")!=std::string::npos) ); }

	private:
		// PRIVATE MEMBERS
		std::string fileName_;
		TFile* outputFile_;
		TTree* eventCountTree_;

		UInt_t numEvtsPass_;
		UInt_t totNumEvtsAnalyzed_;

	protected:
		// PROTECTED MEMBERS
		TTree* mainAnaTree_;
		std::vector<TTree*> otherTrees_;
	};

	TreeHandlerBase::TreeHandlerBase(const std::string& mainTreeName, const std::string& mainTreeDescription, const std::string& fileName) :
		fileName_( fileName ),
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
		for(size_t idx=0; idx<otherTrees_.size(); idx++)
			otherTrees_.at(idx)->Write();
		outputFile_->Close();
	}

	TTree* TreeHandlerBase::addTree(const std::string& treeName, const std::string& treeDescription )
	{
		TTree* treePtr = new TTree( treeName.c_str(), treeDescription.c_str() );
		treePtr->SetDirectory( outputFile_ );
		otherTrees_.push_back(treePtr);
		return treePtr;
	}

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

		TLorentzVector* treeVar_lheZp4Ptr_; TLorentzVector treeVar_lheZp4_;

		TLorentzVector* treeVar_Zp4Ptr_; TLorentzVector treeVar_Zp4_;
		Double_t treeVar_ZpT_;
		Double_t treeVar_Zmass_;
		Double_t treeVar_dR_;
		Double_t treeVar_dEta_;
		Double_t treeVar_dPhi_;

		TLorentzVector* treeVar_eleA_p4Ptr_; TLorentzVector treeVar_eleA_p4_;
		TLorentzVector* treeVar_eleB_p4Ptr_; TLorentzVector treeVar_eleB_p4_;

		Int_t treeVar_eleA_charge_;
		Int_t treeVar_eleB_charge_;
		Int_t treeVar_eleA_frPreCutCode_;
		Int_t treeVar_eleB_frPreCutCode_;
		Int_t treeVar_eleA_stdHeepCutCode_;
		Int_t treeVar_eleB_stdHeepCutCode_;
		Int_t treeVar_eleA_modHeepStdThrCutCode_;
		Int_t treeVar_eleB_modHeepStdThrCutCode_;
		Int_t treeVar_eleA_modHeepColThrCutCode_;
		Int_t treeVar_eleB_modHeepColThrCutCode_;

		Double_t treeVar_eleA_modIsoTk_otherEleVetoForSelf_;
		Double_t treeVar_eleB_modIsoTk_otherEleVetoForSelf_;
		Double_t treeVar_eleA_modIsoEcal_otherEleVetoForSelf_;
		Double_t treeVar_eleB_modIsoEcal_otherEleVetoForSelf_;
		Double_t treeVar_eleA_modIsoHcalD1_otherEleVetoForSelf_;
		Double_t treeVar_eleB_modIsoHcalD1_otherEleVetoForSelf_;


	public:
		DiEleTree(const std::string& treeName, const std::string& fileName, bool fullInfoInTree=false) :
			TreeHandlerBase(treeName, "Tree of Z candidates ("+treeName+")", fileName),
			treeVar_lheZp4Ptr_( &treeVar_lheZp4_ ),
			treeVar_Zp4Ptr_( &treeVar_Zp4_ ),
			treeVar_eleA_p4Ptr_( &treeVar_eleA_p4_ ),
			treeVar_eleB_p4Ptr_( &treeVar_eleB_p4_ )
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

			if(isMadgraphDrellYan(fileName))
				mainAnaTree_->Branch("lheZp4", &treeVar_lheZp4Ptr_);

			mainAnaTree_->Branch("Zp4",   &treeVar_Zp4Ptr_);
			mainAnaTree_->Branch("ZpT",   &treeVar_ZpT_,     "ZpT/D");
			mainAnaTree_->Branch("Zmass", &treeVar_Zmass_,   "Zmass/D");
			mainAnaTree_->Branch("dR",    &treeVar_dR_,      "dR/D");
			mainAnaTree_->Branch("dEta",  &treeVar_dEta_,    "dEta/D");
			mainAnaTree_->Branch("dPhi",  &treeVar_dPhi_,    "dPhi/D");

			mainAnaTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);
			mainAnaTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);

			if(fullInfoInTree){
				mainAnaTree_->Branch("eleA_charge", &treeVar_eleA_charge_, "eleA_charge/I");
				mainAnaTree_->Branch("eleB_charge", &treeVar_eleB_charge_, "eleB_charge/I");
				mainAnaTree_->Branch("eleA_frPre", &treeVar_eleA_frPreCutCode_, "eleA_frPre/I");
				mainAnaTree_->Branch("eleB_frPre", &treeVar_eleB_frPreCutCode_, "eleB_frPre/I");
			}
			mainAnaTree_->Branch("eleA_stdHeep", &treeVar_eleA_stdHeepCutCode_, "eleA_stdHeep/I");
			mainAnaTree_->Branch("eleB_stdHeep", &treeVar_eleB_stdHeepCutCode_, "eleB_stdHeep/I");
			mainAnaTree_->Branch("eleA_modHeepStdThr", &treeVar_eleA_modHeepStdThrCutCode_, "eleA_modHeepStdThr/I");
			mainAnaTree_->Branch("eleB_modHeepStdThr", &treeVar_eleB_modHeepStdThrCutCode_, "eleB_modHeepStdThr/I");
			mainAnaTree_->Branch("eleA_modHeepColThr", &treeVar_eleA_modHeepColThrCutCode_, "eleA_modHeepColThr/I");
			mainAnaTree_->Branch("eleB_modHeepColThr", &treeVar_eleB_modHeepColThrCutCode_, "eleB_modHeepColThr/I");
			if(fullInfoInTree){
				mainAnaTree_->Branch("eleA_modIsoTk_otherEleVetoForSelf",     &treeVar_eleA_modIsoTk_otherEleVetoForSelf_,     "eleA_modIsoTk_otherEleVetoForSelf/D");
				mainAnaTree_->Branch("eleB_modIsoTk_otherEleVetoForSelf",     &treeVar_eleB_modIsoTk_otherEleVetoForSelf_,     "eleB_modIsoTk_otherEleVetoForSelf/D");
				mainAnaTree_->Branch("eleA_modIsoEcal_otherEleVetoForSelf",   &treeVar_eleA_modIsoEcal_otherEleVetoForSelf_,   "eleA_modIsoEcal_otherEleVetoForSelf/D");
				mainAnaTree_->Branch("eleB_modIsoEcal_otherEleVetoForSelf",   &treeVar_eleB_modIsoEcal_otherEleVetoForSelf_,   "eleB_modIsoEcal_otherEleVetoForSelf/D");
				mainAnaTree_->Branch("eleA_modIsoHcalD1_otherEleVetoForSelf", &treeVar_eleA_modIsoHcalD1_otherEleVetoForSelf_, "eleA_modIsoHcalD1_otherEleVetoForSelf/D");
				mainAnaTree_->Branch("eleB_modIsoHcalD1_otherEleVetoForSelf", &treeVar_eleB_modIsoHcalD1_otherEleVetoForSelf_, "eleB_modIsoHcalD1_otherEleVetoForSelf/D");
			}
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

			treeVar_lheZp4_ = evtHelper.GetLHEZbosonP4();

			treeVar_Zp4_   = diEle.p4();
			treeVar_ZpT_   = diEle.pT();
			treeVar_Zmass_ = diEle.invMass();
			treeVar_dR_    = diEle.deltaR();
			treeVar_dEta_  = diEle.deltaEta();
			treeVar_dPhi_  = diEle.deltaPhi();

			treeVar_eleA_p4_ = diEle.eleA().p4();
			treeVar_eleB_p4_ = diEle.eleB().p4();

			treeVar_eleA_charge_ = diEle.eleA().charge();
			treeVar_eleB_charge_ = diEle.eleB().charge();

			treeVar_eleA_frPreCutCode_ = diEle.eleA().fakeRatePreSelnCutCode();
			treeVar_eleB_frPreCutCode_ = diEle.eleB().fakeRatePreSelnCutCode();

			treeVar_eleA_stdHeepCutCode_ = diEle.eleA().heepIdStdIsoCutCode(evtHelper);
			treeVar_eleB_stdHeepCutCode_ = diEle.eleB().heepIdStdIsoCutCode(evtHelper);

			treeVar_eleA_modHeepStdThrCutCode_ = diEle.eleA().heepIdModIsoStdThrCutCode(evtHelper);
			treeVar_eleB_modHeepStdThrCutCode_ = diEle.eleB().heepIdModIsoStdThrCutCode(evtHelper);

			treeVar_eleA_modHeepColThrCutCode_ = diEle.eleA().heepIdModIsoColThrCutCode(evtHelper);
			treeVar_eleB_modHeepColThrCutCode_ = diEle.eleB().heepIdModIsoColThrCutCode(evtHelper);

			treeVar_eleA_modIsoTk_otherEleVetoForSelf_     = diEle.eleA().modIsoTk_otherEleVetoForSelf();
			treeVar_eleB_modIsoTk_otherEleVetoForSelf_     = diEle.eleB().modIsoTk_otherEleVetoForSelf();
			treeVar_eleA_modIsoEcal_otherEleVetoForSelf_   = diEle.eleA().modIsoEcal_otherEleVetoForSelf();
			treeVar_eleB_modIsoEcal_otherEleVetoForSelf_   = diEle.eleB().modIsoEcal_otherEleVetoForSelf();
			treeVar_eleA_modIsoHcalD1_otherEleVetoForSelf_ = diEle.eleA().modIsoHcalD1_otherEleVetoForSelf();
			treeVar_eleB_modIsoHcalD1_otherEleVetoForSelf_ = diEle.eleB().modIsoHcalD1_otherEleVetoForSelf();

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

		Double_t treeVar_genWeight_;
		Double_t treeVar_puWeight_;
		float treeVar_mc_numVtx_;

		TLorentzVector* treeVar_lheZp4Ptr_; TLorentzVector treeVar_lheZp4_;

		TLorentzVector* treeVar_mcZ_ele1_p4Ptr_; TLorentzVector treeVar_mcZ_ele1_p4_;
		TLorentzVector* treeVar_mcZ_ele2_p4Ptr_; TLorentzVector treeVar_mcZ_ele2_p4_;
		UInt_t treeVar_detRegion_;
		Bool_t treeVar_bothRecod_;

		TLorentzVector* treeVar_Zp4Ptr_; TLorentzVector treeVar_Zp4_;
		Double_t treeVar_ZpT_;
		Double_t treeVar_ZdEta_;
		Double_t treeVar_ZdPhi_;
		Double_t treeVar_ZdR_;
		TLorentzVector* treeVar_eleA_p4Ptr_; TLorentzVector treeVar_eleA_p4_;
		Double_t treeVar_eleA_scEta_;
		Double_t treeVar_eleA_dRmc_;
		TLorentzVector* treeVar_eleB_p4Ptr_;	TLorentzVector treeVar_eleB_p4_;
		Double_t treeVar_eleB_scEta_;
		Double_t treeVar_eleB_dRmc_;


	public:
		EffiCalcTree(const std::string& fileName) :
			TreeHandlerBase("zBosonEffiTree", "Tree of Z candidate information for signal MC effi calc'ns", fileName),
			treeVar_lheZp4Ptr_( &treeVar_lheZp4_ ),
			treeVar_mcZ_ele1_p4Ptr_( &treeVar_mcZ_ele1_p4_ ),
			treeVar_mcZ_ele2_p4Ptr_( &treeVar_mcZ_ele2_p4_ ),
			treeVar_Zp4Ptr_( &treeVar_Zp4_ ),
			treeVar_eleA_p4Ptr_( &treeVar_eleA_p4_ ),
			treeVar_eleB_p4Ptr_( &treeVar_eleB_p4_ )
		{
			// Setting up the event / di-ele branches ...
			mainAnaTree_->Branch("genWeight",      &treeVar_genWeight_, "genWeight/D");
			mainAnaTree_->Branch("puWeight",      &treeVar_puWeight_, "puWeight/D");
			mainAnaTree_->Branch("mc_numVtx", &treeVar_mc_numVtx_, "mc_numVtx/f");
			mainAnaTree_->Branch("run", &treeVar_runNum_, "run/i");
			mainAnaTree_->Branch("lumi", &treeVar_lumiSec_, "lumi/i");
			mainAnaTree_->Branch("evtNum", &treeVar_evtNum_, "evtNum/i");

			if(isMadgraphDrellYan(fileName))
				mainAnaTree_->Branch("lheZp4", &treeVar_lheZp4_);

			mainAnaTree_->Branch("mcZ_ele1_p4", &treeVar_mcZ_ele1_p4Ptr_);
			mainAnaTree_->Branch("mcZ_ele2_p4", &treeVar_mcZ_ele2_p4Ptr_);
			mainAnaTree_->Branch("mc_detRegion",   &treeVar_detRegion_, "mcAccept_pt/i");
			mainAnaTree_->Branch("bothRecod", &treeVar_bothRecod_, "bothRecod/O");

			mainAnaTree_->Branch("Zp4", &treeVar_Zp4Ptr_);
			mainAnaTree_->Branch("ZpT", &treeVar_ZpT_,     "ZpT/D");
			mainAnaTree_->Branch("ZdEta", &treeVar_ZdEta_,     "ZdEta/D");
			mainAnaTree_->Branch("ZdPhi", &treeVar_ZdPhi_,     "ZdPhi/D");
			mainAnaTree_->Branch("ZdR", &treeVar_ZdR_, "ZdR/D");
			mainAnaTree_->Branch("eleA_p4", &treeVar_eleA_p4Ptr_);
			mainAnaTree_->Branch("eleA_scEta", &treeVar_eleA_scEta_, "eleA_scEta/D");
			mainAnaTree_->Branch("eleA_dRmc", &treeVar_eleA_dRmc_, "eleA_dRmc/D");
			mainAnaTree_->Branch("eleB_p4", &treeVar_eleB_p4Ptr_);
			mainAnaTree_->Branch("eleB_scEta", &treeVar_eleB_scEta_, "eleB_scEta/D");
			mainAnaTree_->Branch("eleB_dRmc", &treeVar_eleB_dRmc_, "eleB_dRmc/D");


			setupBranches("eleA", eleA_vars_);
			setupBranches("eleB", eleB_vars_);

		}
		~EffiCalcTree(){ }

		void fillTree(const tsw::HEEPDiEle& diEle, const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2, const tsw::EventHelper& eventHelper)
		{
//			std::cout << std::endl << " *** Entering tsw::EffiCalcTree::fillTree ***" << std::endl;
//			std::cout << "Event (run " << eventHelper.runNum() << ", lumiSec " << eventHelper.lumiSec() << ", event number " << eventHelper.eventNum() << ")" << std::endl
//						 << " * Di-ele:" << std::endl
//						 << "      pT(ee) = " << diEle.pT() << " GeV  , dR_ee = " << diEle.deltaR() << std::endl;
//			std::cout << " * Leading ele:" << std::endl
//						 << "      Et = " << diEle.eleA().et() << "GeV,  eta = " << diEle.eleA().eta() << ",  phi = " << diEle.eleA().phi() << std::endl
//						 << "      Mod ECAL iso = " << diEle.eleA().isol_inrVetoModEm(tsw::Event::mediumVeto) << "GeV ,  mod HCAL iso = " << diEle.eleA().isol_inrVetoModHadD1(tsw::Event::mediumVeto) << "GeV" << std::endl
//						 << "      Mod track iso = " << diEle.eleA().isol_inrVetoModTrk(tsw::Event::xSmallVeto) << "GeV" << std::endl;
//			std::cout << " * Sub-leading ele:" << std::endl
//						 << "      Et = " << diEle.eleB().et() << "GeV,  eta = " << diEle.eleB().eta() << ",  phi = " << diEle.eleB().phi() << std::endl
//						 << "      Mod ECAL iso = " << diEle.eleB().isol_inrVetoModEm(tsw::Event::mediumVeto) << "GeV ,  mod HCAL iso = " << diEle.eleB().isol_inrVetoModHadD1(tsw::Event::mediumVeto) << "GeV" << std::endl
//						 << "      Mod track iso = " << diEle.eleB().isol_inrVetoModTrk(tsw::Event::xSmallVeto) << "GeV" << std::endl;
//			if( diEle.eleA().heepIdModIsoStdThrCut(eventHelper) && diEle.eleB().heepIdModIsoStdThrCut(eventHelper))
//				std::cout << " Both electrons pass HEEP v4.1 ID & mod iso cuts" << std::endl;
//			else{
//				if( !diEle.eleA().heepIdModIsoStdThrCut(eventHelper) )
//					std::cout << " Leading ele failed selection" << std::endl;
//				if( !diEle.eleB().heepIdModIsoStdThrCut(eventHelper) )
//					std::cout << " Sub-leading ele failed selection" << std::endl;
//			}
			// Set branches with common values between two cases (Z boson reco'd / not )
			fillTree_common(mcZboson_ele1, mcZboson_ele2, eventHelper);

			// ... then set values of other branches
			treeVar_bothRecod_ = true;

			treeVar_Zp4_ = diEle.p4();
			treeVar_ZpT_ = diEle.pT();
			treeVar_ZdEta_ = diEle.deltaEta();
			treeVar_ZdPhi_ = diEle.deltaPhi();
			treeVar_ZdR_   = diEle.deltaR();
			treeVar_eleA_p4_ = diEle.eleA().p4();
			treeVar_eleB_p4_ = diEle.eleB().p4();

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
			treeVar_eleA_scEta_ = diEle.eleA().scEta();
			treeVar_eleB_scEta_ = diEle.eleB().scEta();

			setBranchValues(eleA_vars_, diEle.eleA(), eventHelper, diEle.eleA_modTrkIso(), diEle.eleA_modEmHad1Iso() );
			setBranchValues(eleB_vars_, diEle.eleB(), eventHelper, diEle.eleB_modTrkIso(), diEle.eleB_modEmHad1Iso() );

			// And finally fill the tree ...
			mainAnaTree_->Fill();
		}

		void fillTree_NotReconstructed(const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2, const tsw::EventHelper& eventHelper)
		{
			fillTree_common(mcZboson_ele1, mcZboson_ele2, eventHelper);

			treeVar_bothRecod_ = false;

			treeVar_Zp4_ = mcZboson_ele1+mcZboson_ele2;
			treeVar_ZpT_ = treeVar_Zp4_.Pt();
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
			treeVar_eleA_scEta_ = treeVar_eleA_p4_.Eta();
			treeVar_eleB_scEta_ = treeVar_eleB_p4_.Eta();

			// Default value setting from now on ...
			eleA_vars_.resetToDefault();
			eleB_vars_.resetToDefault();

			// And finally fill the tree ...
			mainAnaTree_->Fill();
		}

	private:
		// PRIVATE STRUCTS //
		struct RecoEleVars{
			Int_t mStdHeepCutCode;
			Int_t mModHeepStdThrCutCode;
			Int_t mModHeepColThrCutCode;

			Double_t mStdTrkIso;
			Double_t mStdEmH1Iso;

			unsigned int mNTrksInnerVeto;

			Double_t mModTrkIso;
			Double_t mScModEmH1Iso;

			Double_t mIsoDep_stdTrk;
			Double_t mIsoDep_stdEmH1;
			Double_t mIsoDep_inrVetoModTrk;
			Double_t mIsoDep_inrVetoModEmH1;

			Double_t mInrXSVetoModTrk;
			Double_t mInrSVetoModTrk;
			Double_t mInrMVetoModTrk;
			Double_t mInrLVetoModTrk;
			Double_t mInrXLVetoModTrk;

			Double_t mInrXSVetoModEmH1;
			Double_t mInrSVetoModEmH1;
			Double_t mInrMVetoModEmH1;
			Double_t mInrLVetoModEmH1;
			Double_t mInrXLVetoModEmH1;

			UInt_t mNGenHadronsDr04;
			Double_t mPtSumGenHadronsDr04;

			Double_t mEmH1RhoCorrn;

			// PUBLIC METHODS
			RecoEleVars() { resetToDefault(); }
			void resetToDefault();
		};


		// PRIVATE METHODS //
		void setupBranches(const std::string& prefix, RecoEleVars& eleVarsStruct);
		void setBranchValues(RecoEleVars& eleVarsStruct, const tsw::HEEPEle& ele, const tsw::EventHelper& eventHelper, const double diEle_modTrkIso, const double diEle_modEmH1Iso);
		void fillTree_common(const TLorentzVector& mcZboson_ele1, const TLorentzVector& mcZboson_ele2, const tsw::EventHelper& eventHelper)
		{
			treeVar_genWeight_  = eventHelper.genWeight();
			treeVar_puWeight_   = eventHelper.puWeight();
			treeVar_mc_numVtx_  = eventHelper.GetMCPU_nVtx();
			treeVar_runNum_     = eventHelper.runNum();
			treeVar_lumiSec_    = eventHelper.lumiSec();
			treeVar_evtNum_     = eventHelper.eventNum();

			treeVar_lheZp4_ = eventHelper.GetLHEZbosonP4();

			treeVar_mcZ_ele1_p4_ = mcZboson_ele1;
			treeVar_mcZ_ele2_p4_ = mcZboson_ele2;

 			bool mcZele1_isEB = fabs(mcZboson_ele1.Eta())<1.4;
			bool mcZele2_isEB = fabs(mcZboson_ele2.Eta())<1.4;
			bool mcZele1_isEE = (fabs(mcZboson_ele1.Eta())>1.6 && fabs(mcZboson_ele1.Eta())<2.45);
			bool mcZele2_isEE = (fabs(mcZboson_ele2.Eta())>1.6 && fabs(mcZboson_ele2.Eta())<2.45);

			if( mcZele1_isEB && mcZele2_isEB )
				treeVar_detRegion_ = 0;
			else if( (mcZele1_isEB && mcZele2_isEE) || (mcZele1_isEE && mcZele2_isEB) )
				treeVar_detRegion_ = 1;
			else if( mcZele1_isEE && mcZele2_isEE )
				treeVar_detRegion_ = 2;
			else
				treeVar_detRegion_ = 3;

		}

		// PRIVATE MEMBERS //
		RecoEleVars eleA_vars_;
		RecoEleVars eleB_vars_;
	};

	void EffiCalcTree::RecoEleVars::resetToDefault()
	{
		mStdHeepCutCode       = 0x01000-1;
		mModHeepStdThrCutCode = 0x01000-1;
		mModHeepColThrCutCode = 0x01000-1;

		mStdTrkIso  = 9999.9;
		mModTrkIso  = 9999.9;
		mStdEmH1Iso = 9999.9;

		mNTrksInnerVeto = 9999.9;
		mScModEmH1Iso   = 9999.9;

		mIsoDep_stdTrk         = 9999.9;
		mIsoDep_stdEmH1        = 9999.9;
		mIsoDep_inrVetoModTrk  = 9999.9;
		mIsoDep_inrVetoModEmH1 = 9999.9;

		mInrXSVetoModTrk  = 9999.9;
		mInrSVetoModTrk   = 9999.9;
		mInrMVetoModTrk   = 9999.9;
		mInrLVetoModTrk   = 9999.9;
		mInrXLVetoModTrk  = 9999.9;

		mInrXSVetoModEmH1 = 9999.9;
		mInrSVetoModEmH1  = 9999.9;
		mInrMVetoModEmH1  = 9999.9;
		mInrLVetoModEmH1  = 9999.9;
		mInrXLVetoModEmH1 = 9999.9;

		mNGenHadronsDr04     = 9999;
		mPtSumGenHadronsDr04 = 9999.9;

		mEmH1RhoCorrn = 0.0;
	}

	void EffiCalcTree::setupBranches(const std::string& prefix, RecoEleVars& eleVarsStruct)
	{

		mainAnaTree_->Branch((prefix+"_stdHeep").c_str(), &(eleVarsStruct.mStdHeepCutCode), (prefix+"_stdHeep/I").c_str());
		mainAnaTree_->Branch((prefix+"_modHeepStdThr").c_str(), &(eleVarsStruct.mModHeepStdThrCutCode), (prefix+"_modHeepStdThr/I").c_str());
		mainAnaTree_->Branch((prefix+"_modHeepColThr").c_str(), &(eleVarsStruct.mModHeepColThrCutCode), (prefix+"_modHeepColThr/I").c_str());


		mainAnaTree_->Branch((prefix+"_stdTrkIso").c_str(),  &(eleVarsStruct.mStdTrkIso),  (prefix+"_stdTrkIso/D").c_str());
		mainAnaTree_->Branch((prefix+"_stdEmH1Iso").c_str(), &(eleVarsStruct.mStdEmH1Iso), (prefix+"_stdEmH1Iso/D").c_str());

		mainAnaTree_->Branch((prefix+"_nTrksInnerVeto").c_str(), &(eleVarsStruct.mNTrksInnerVeto), (prefix+"_nTrksInnerVeto/i").c_str());

		mainAnaTree_->Branch((prefix+"_modTrkIso").c_str(), &(eleVarsStruct.mModTrkIso), (prefix+"_modTrkIso/D").c_str());
		mainAnaTree_->Branch((prefix+"_scModEmH1Iso").c_str(), &(eleVarsStruct.mScModEmH1Iso), (prefix+"_scModEmH1Iso/D").c_str());


		mainAnaTree_->Branch((prefix+"_isoDep_stdTrk").c_str(),      &(eleVarsStruct.mIsoDep_stdTrk),  (prefix+"_isoDep_stdTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_isoDep_stdEmH1").c_str(),     &(eleVarsStruct.mIsoDep_stdEmH1), (prefix+"_isoDep_stdEmH1/D").c_str());
		mainAnaTree_->Branch((prefix+"_isoDep_inrVetoModTrk").c_str(),     &(eleVarsStruct.mIsoDep_inrVetoModTrk), (prefix+"_isoDep_inrVetoModTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_isoDep_inrVetoModEmH1").c_str(),    &(eleVarsStruct.mIsoDep_inrVetoModEmH1), (prefix+"_isoDep_inrVetoModEmH1/D").c_str());

		mainAnaTree_->Branch((prefix+"_inrXSVetoModTrk").c_str(), &(eleVarsStruct.mInrXSVetoModTrk), (prefix+"_inrXSVetoModTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrSVetoModTrk").c_str(),  &(eleVarsStruct.mInrSVetoModTrk),  (prefix+"_inrSVetoModTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrMVetoModTrk").c_str(),  &(eleVarsStruct.mInrMVetoModTrk),  (prefix+"_inrMVetoModTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrLVetoModTrk").c_str(),  &(eleVarsStruct.mInrLVetoModTrk),  (prefix+"_inrLVetoModTrk/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrXLVetoModTrk").c_str(), &(eleVarsStruct.mInrXLVetoModTrk), (prefix+"_inrXLVetoModTrk/D").c_str());

		mainAnaTree_->Branch((prefix+"_inrXSVetoModEmH1").c_str(), &(eleVarsStruct.mInrXSVetoModEmH1), (prefix+"_inrXSVetoModEmH1/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrSVetoModEmH1").c_str(),  &(eleVarsStruct.mInrSVetoModEmH1),  (prefix+"_inrSVetoModEmH1/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrMVetoModEmH1").c_str(),  &(eleVarsStruct.mInrMVetoModEmH1),  (prefix+"_inrMVetoModEmH1/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrLVetoModEmH1").c_str(),  &(eleVarsStruct.mInrLVetoModEmH1),  (prefix+"_inrLVetoModEmH1/D").c_str());
		mainAnaTree_->Branch((prefix+"_inrXLVetoModEmH1").c_str(), &(eleVarsStruct.mInrXLVetoModEmH1), (prefix+"_inrXLVetoModEmH1/D").c_str());

		mainAnaTree_->Branch((prefix+"_nGenHadronsDr04").c_str(), &(eleVarsStruct.mNGenHadronsDr04), (prefix+"_nGenHadronsDr04/i").c_str());
		mainAnaTree_->Branch((prefix+"_ptSumGenHadronsDr04").c_str(), &(eleVarsStruct.mPtSumGenHadronsDr04), (prefix+"_ptSumGenHadronsDr04/D").c_str());

		mainAnaTree_->Branch((prefix+"_EmH1RhoCorrn").c_str(), &(eleVarsStruct.mEmH1RhoCorrn), (prefix+"_EmH1RhoCorrn_/D").c_str());
	}

	void EffiCalcTree::setBranchValues(RecoEleVars& eleVarsStruct, const tsw::HEEPEle& ele, const tsw::EventHelper& eventHelper, const double diEle_modTrkIso, const double diEle_modEmH1Iso)
	{
		eleVarsStruct.mStdHeepCutCode       = ele.heepIdStdIsoCutCode(eventHelper);
		eleVarsStruct.mModHeepStdThrCutCode = ele.heepIdModIsoStdThrCutCode(eventHelper);
		eleVarsStruct.mModHeepColThrCutCode = ele.heepIdModIsoColThrCutCode(eventHelper);

		eleVarsStruct.mStdTrkIso     = ele.isolPtTrks();
		eleVarsStruct.mModTrkIso     = diEle_modTrkIso;
		eleVarsStruct.mStdEmH1Iso    = ele.isolEmHadDepth1();

		eleVarsStruct.mNTrksInnerVeto = ele.numInnerIsoConeTrks();
		eleVarsStruct.mScModEmH1Iso   = diEle_modEmH1Iso;

		// Alternative isolation values, from the isoDeps modules ...
		eleVarsStruct.mIsoDep_stdTrk  = ele.isol_isoDep_stdTrk();
		eleVarsStruct.mIsoDep_stdEmH1 = ele.isol_isoDep_stdEmHadD1();
		eleVarsStruct.mIsoDep_inrVetoModTrk  = ele.isol_isoDep_inrVetoModTrk();
		eleVarsStruct.mIsoDep_inrVetoModEmH1 = ele.isol_isoDep_inrVetoModEmHadD1();

		// Alternative isolation values, from the BstdZee EDProducer ...
		eleVarsStruct.mInrXSVetoModTrk  = ele.isol_inrVetoModTrk(tsw::Event::xSmallVeto);
		eleVarsStruct.mInrSVetoModTrk   = ele.isol_inrVetoModTrk(tsw::Event::smallVeto);
		eleVarsStruct.mInrMVetoModTrk   = ele.isol_inrVetoModTrk(tsw::Event::mediumVeto);
		eleVarsStruct.mInrLVetoModTrk   = ele.isol_inrVetoModTrk(tsw::Event::largeVeto);
		eleVarsStruct.mInrXLVetoModTrk  = ele.isol_inrVetoModTrk(tsw::Event::xLargeVeto);

		eleVarsStruct.mInrXSVetoModEmH1 = ele.isol_inrVetoModEmHadD1(tsw::Event::xSmallVeto);
		eleVarsStruct.mInrSVetoModEmH1  = ele.isol_inrVetoModEmHadD1(tsw::Event::smallVeto);
		eleVarsStruct.mInrMVetoModEmH1  = ele.isol_inrVetoModEmHadD1(tsw::Event::mediumVeto);
		eleVarsStruct.mInrLVetoModEmH1  = ele.isol_inrVetoModEmHadD1(tsw::Event::largeVeto);
		eleVarsStruct.mInrXLVetoModEmH1 = ele.isol_inrVetoModEmHadD1(tsw::Event::xLargeVeto);

		eleVarsStruct.mNGenHadronsDr04     = ele.isol_nGenHadronsDr04();
		eleVarsStruct.mPtSumGenHadronsDr04 = ele.isol_ptSumGenHadronsDr04();

		// Rho PU correction to EmH1 isol'n values
		eleVarsStruct.mEmH1RhoCorrn = ele.isol_rhoCorrnEmH1(eventHelper);
	}


	///////////////////////////////////////////////////////////////////////////////
	// TagProbeTree: Used to generate the trees required for tag & probe studies
	//               [1 tag-probe pair tree & 1 gsf-gsf pair tree (same seln) for QCD background estimations]

	class TagProbeTree : public TreeHandlerBase {

	public:
		TagProbeTree(const std::string& tpTreeName, const std::string& gsfGsfTreeName, const std::string& fileName) :
			TreeHandlerBase(tpTreeName, "Tree of tag-probe pairs ("+tpTreeName+")", fileName),
			tpTreeBranchVars_(),
			gsfGsfTree_( addTree(gsfGsfTreeName, "Tree of GSF-GSF pairs for QCD estimate ("+gsfGsfTreeName+")") ),
			gsfGsfTreeBranchVars_()
		{
			setupBranchLinks(mainAnaTree_, tpTreeBranchVars_);
			setupBranchLinks(gsfGsfTree_, gsfGsfTreeBranchVars_);
		}
		~TagProbeTree(){}

		void fillTagProbeTree(const tsw::HEEPEle& tagEle, const tsw::HEEPEle& probeEle, const tsw::EventHelper& evtHelper, bool trigDecision)
		{
			fillTree(mainAnaTree_, tpTreeBranchVars_, tagEle, probeEle, evtHelper, trigDecision);
		}

		void fillQcdTree(const tsw::HEEPEle& tagCandidate, const tsw::HEEPEle& probeCandidate, const tsw::EventHelper& evtHelper, bool trigDecision)
		{
			fillTree(gsfGsfTree_, gsfGsfTreeBranchVars_, tagCandidate, probeCandidate, evtHelper, trigDecision);
		}

	private:
		// PRIVATE STRUCTS

		struct ModIsoVarsWithPhantom{
			Double_t trk, emHad1;
			Double_t dEta, dPhi, dR;
		};

		struct TreeBranchVars{
			// CTOR
			TreeBranchVars() :
				lheZp4Ptr_(&lheZp4_),
				pair_p4Ptr_(&pair_p4_),  tag_p4Ptr_(&tag_p4_),  prb_p4Ptr_(&prb_p4_)
			{}
			// PUBLIC MEMBERS
			Double_t genWeight_;  Double_t puWeight_;

			UInt_t runNum_;  UInt_t lumiSec_;  UInt_t evtNum_;

			Bool_t trgDecision_;
			UInt_t nVtx_;

			TLorentzVector* lheZp4Ptr_; TLorentzVector lheZp4_;

			TLorentzVector* pair_p4Ptr_; TLorentzVector pair_p4_;
			Double_t pair_dEta_, pair_dPhi_, pair_dR_;

			TLorentzVector* tag_p4Ptr_; TLorentzVector tag_p4_;
			Double_t tag_scEta_;
			Int_t tag_stdHeepCutCode_;
			Int_t tag_modHeepStdThrCutCode_;
			Int_t tag_modHeepColThrCutCode_;
			Int_t tag_fakePreCutCode_;
			Int_t tag_charge_;

			TLorentzVector* prb_p4Ptr_; TLorentzVector prb_p4_;
			Double_t prb_scEta_;
			Int_t prb_stdHeepCutCode_;
			Int_t prb_modHeepStdThrCutCode_;
			Int_t prb_modHeepColThrCutCode_;
			Int_t prb_fakePreCutCode_;
			Int_t prb_charge_;
			Double_t prb_stdTrkIso_;
			Double_t prb_stdEmHad1Iso_;
			Double_t prb_calEaCorr_;

			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr005To010;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr010To015;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr015To020;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr020To025;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr025To030;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr030To035;
			ModIsoVarsWithPhantom prbInrVetoModIsoPhantomDr035To040;
		};

		// PRIVATE METHODS
		void setupBranchLinks(TTree* treePtr, TreeBranchVars& treeBranchVars);
		static void setupBranchLinks(TTree* treePtr, ModIsoVarsWithPhantom& phantomModIsoStruct, const std::string& suffix);
		static void fillTree(TTree* treePtr, TreeBranchVars& treeVars, const tsw::HEEPEle& tagEle, const tsw::HEEPEle& probeEle, const tsw::EventHelper& evtHelper, bool trigDecision);
		static void setBranchValues(ModIsoVarsWithPhantom& withPhantomStruct, const tsw::HEEPEle& ele, const tsw::HEEPEle::ModIsoPhantomEleDrRange drRange);

		// PRIVATE MEMBERS
		TreeBranchVars tpTreeBranchVars_;

		TTree* gsfGsfTree_;
		TreeBranchVars gsfGsfTreeBranchVars_;

	};

	void TagProbeTree::setupBranchLinks(TTree* treePtr, TreeBranchVars& branchVars)
	{
		treePtr->Branch("genWeight",  &(branchVars.genWeight_),  "genWeight/D");
		treePtr->Branch("puWeight",   &(branchVars.puWeight_),   "puWeight/D");

		treePtr->Branch("runNum",  &branchVars.runNum_,   "runNum/i");
		treePtr->Branch("lumiSec", &branchVars.lumiSec_,  "lumiSec/i");
		treePtr->Branch("evtNum",  &branchVars.evtNum_,   "evtNum/i");

		treePtr->Branch("trgDecision", &(branchVars.trgDecision_), "trgDecision/O");
		treePtr->Branch("nVtx", &(branchVars.nVtx_), "nVtx/i");

		if( isMadgraphDrellYan( fileName() ) )
			treePtr->Branch("lheZp4", &(branchVars.lheZp4Ptr_) );

		treePtr->Branch("pair_p4",   &(branchVars.pair_p4Ptr_) );
		treePtr->Branch("pair_dEta", &(branchVars.pair_dEta_), "pair_dEta/D");
		treePtr->Branch("pair_dPhi", &(branchVars.pair_dPhi_), "pair_dPhi/D");
		treePtr->Branch("pair_dR",   &(branchVars.pair_dR_),   "pair_dR/D");

		treePtr->Branch("tag_p4",   &(branchVars.tag_p4Ptr_) );
		treePtr->Branch("tag_scEta",   &(branchVars.tag_scEta_), "tag_scEta/D" );
		treePtr->Branch("tag_stdHeep", &(branchVars.tag_stdHeepCutCode_), "tag_stdHeep/I");
		treePtr->Branch("tag_modHeepStdThr", &(branchVars.tag_modHeepStdThrCutCode_), "tag_modHeepStdThr/I");
		treePtr->Branch("tag_modHeepColThr", &(branchVars.tag_modHeepColThrCutCode_), "tag_modHeepColThr/I");
		treePtr->Branch("tag_fakePreCutCode", &(branchVars.tag_fakePreCutCode_), "tag_fakePreCutCode/I");
		treePtr->Branch("tag_charge", &(branchVars.tag_charge_), "tag_charge/I");

		treePtr->Branch("probe_p4", &(branchVars.prb_p4Ptr_) );
		treePtr->Branch("probe_scEta",   &(branchVars.prb_scEta_), "probe_scEta/D" );
		treePtr->Branch("probe_stdHeep", &(branchVars.prb_stdHeepCutCode_), "probe_stdHeep/I");
		treePtr->Branch("probe_modHeepStdThr", &(branchVars.prb_modHeepStdThrCutCode_), "probe_modHeepStdThr/I");
		treePtr->Branch("probe_modHeepColThr", &(branchVars.prb_modHeepColThrCutCode_), "probe_modHeepColThr/I");
		treePtr->Branch("probe_fakePreCutCode", &(branchVars.prb_fakePreCutCode_), "probe_fakePreCutCode/I");
		treePtr->Branch("probe_charge", &(branchVars.prb_charge_), "probe_charge/I");
		treePtr->Branch("probe_stdTrkIso", &(branchVars.prb_stdTrkIso_), "probe_stdTrkIso/D");
		treePtr->Branch("probe_stdEmHad1Iso", &(branchVars.prb_stdEmHad1Iso_), "probe_stdEmHad1Iso/D");
		treePtr->Branch("probe_calEaCorr", &(branchVars.prb_calEaCorr_), "probe_calEaCorr/D");

		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr005To010, "005To010");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr010To015, "010To015");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr015To020, "015To020");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr020To025, "020To025");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr025To030, "025To030");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr030To035, "030To035");
		setupBranchLinks(treePtr, branchVars.prbInrVetoModIsoPhantomDr035To040, "035To040");
	}

	void TagProbeTree::setupBranchLinks(TTree* treePtr, ModIsoVarsWithPhantom& phantomModIsoStruct, const std::string& suffix)
	{
		treePtr->Branch(("probe_modTrkIsoPhantom_"+suffix).c_str(),    &(phantomModIsoStruct.trk),    ("probe_modTrkIsoPhantom_"+suffix+"/D").c_str() );
		treePtr->Branch(("probe_modEmHad1IsoPhantom_"+suffix).c_str(), &(phantomModIsoStruct.emHad1), ("probe_modEmHad1IsoPhantom_"+suffix+"/D").c_str() );
		treePtr->Branch(("probe_modIsoPhantomDEta_"+suffix).c_str(), &(phantomModIsoStruct.dEta), ("probe_modIsoPhantomDEta_"+suffix+"/D").c_str() );
		treePtr->Branch(("probe_modIsoPhantomDPhi_"+suffix).c_str(), &(phantomModIsoStruct.dPhi), ("probe_modIsoPhantomDPhi_"+suffix+"/D").c_str() );
		treePtr->Branch(("probe_modIsoPhantomDR_"+suffix).c_str(),   &(phantomModIsoStruct.dR),   ("probe_modIsoPhantomDR_"+suffix+"/D").c_str() );
	}

	void TagProbeTree::fillTree(TTree* treePtr, TreeBranchVars& treeVars, const tsw::HEEPEle& tagEle, const tsw::HEEPEle& probeEle, const tsw::EventHelper& evtHelper, bool trigDecision)
	{
		// Set branch variables ...
		treeVars.genWeight_ = evtHelper.genWeight();
		treeVars.puWeight_  = evtHelper.puWeight();

		treeVars.runNum_ = evtHelper.runNum();
		treeVars.lumiSec_ = evtHelper.lumiSec();
		treeVars.evtNum_ = evtHelper.eventNum();

		treeVars.trgDecision_ = trigDecision;
		treeVars.nVtx_  = evtHelper.GetRecoVtx_nGoodVtxs();

		treeVars.lheZp4_ = evtHelper.GetLHEZbosonP4();

		treeVars.pair_p4_ = ( tagEle.p4() + probeEle.p4() );
		treeVars.pair_dEta_ = probeEle.eta() - tagEle.eta();
		treeVars.pair_dPhi_ = probeEle.p4().DeltaPhi(tagEle.p4());
		treeVars.pair_dR_   = sqrt( pow(treeVars.pair_dEta_,2.0) + pow(treeVars.pair_dPhi_,2.0) );

		treeVars.tag_p4_    = tagEle.p4();
		treeVars.tag_scEta_ = tagEle.scEta();
		treeVars.tag_stdHeepCutCode_ = tagEle.heepIdStdIsoCutCode(evtHelper);
		treeVars.tag_modHeepStdThrCutCode_ = tagEle.heepIdModIsoStdThrCutCode(evtHelper);
		treeVars.tag_modHeepColThrCutCode_ = tagEle.heepIdModIsoColThrCutCode(evtHelper);
		treeVars.tag_fakePreCutCode_ = tagEle.fakeRatePreSelnCutCode();
		treeVars.tag_charge_ = tagEle.charge();

		treeVars.prb_p4_        = probeEle.p4();
		treeVars.prb_scEta_     = probeEle.scEta();
		treeVars.prb_stdHeepCutCode_ = probeEle.heepIdStdIsoCutCode(evtHelper);
		treeVars.prb_modHeepStdThrCutCode_ = probeEle.heepIdModIsoStdThrCutCode(evtHelper);
		treeVars.prb_modHeepColThrCutCode_ = probeEle.heepIdModIsoColThrCutCode(evtHelper);
		treeVars.prb_fakePreCutCode_ = probeEle.fakeRatePreSelnCutCode();
		treeVars.prb_charge_ = probeEle.charge();
		treeVars.prb_stdTrkIso_     = probeEle.isolPtTrks();
		treeVars.prb_stdEmHad1Iso_  = probeEle.isolEmHadDepth1();
		treeVars.prb_calEaCorr_     = probeEle.isol_rhoCorrnEmH1(evtHelper);

		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr005To010 , probeEle, tsw::HEEPEle::PHANTOM_DR_005_010 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr010To015 , probeEle, tsw::HEEPEle::PHANTOM_DR_010_015 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr015To020 , probeEle, tsw::HEEPEle::PHANTOM_DR_015_020 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr020To025 , probeEle, tsw::HEEPEle::PHANTOM_DR_020_025 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr025To030 , probeEle, tsw::HEEPEle::PHANTOM_DR_025_030 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr030To035 , probeEle, tsw::HEEPEle::PHANTOM_DR_030_035 );
		setBranchValues(treeVars.prbInrVetoModIsoPhantomDr035To040 , probeEle, tsw::HEEPEle::PHANTOM_DR_035_040 );

		// And finally fill the tree ...
		treePtr->Fill();
	}

	void TagProbeTree::setBranchValues(ModIsoVarsWithPhantom& withPhantomStruct, const tsw::HEEPEle& ele, const tsw::HEEPEle::ModIsoPhantomEleDrRange drRange)
	{
		withPhantomStruct.trk    = ele.modVetoIsoWithPhantomEle(drRange).trk;
		withPhantomStruct.emHad1 = ele.modVetoIsoWithPhantomEle(drRange).ecal + ele.modVetoIsoWithPhantomEle(drRange).hcalD1;
		withPhantomStruct.dEta   = ele.modVetoIsoWithPhantomEle(drRange).dEta;
		withPhantomStruct.dPhi   = ele.modVetoIsoWithPhantomEle(drRange).dPhi;
		withPhantomStruct.dR     = ele.modVetoIsoWithPhantomEle(drRange).dR();
	}


	///////////////////////////////////////////////////////////////////////////////
	// EleMuTree: Simple class for generating trees containing info about ele-mu objects
	//
	class EleMuTree : public TreeHandlerBase {
	private:
		// PRIVATE MEMBERS

		// For the event-level branches ...
		Double_t evt_weight_;
		Double_t evt_genWeight_;
		Double_t evt_puWeight_;
		Double_t evt_xsecWeight_;

		UInt_t evt_runNum_;
		UInt_t evt_lumiNum_;
		UInt_t evt_evtNum_;

		Bool_t evt_trgDecision_;
		UInt_t evt_nVtx_;

		TLorentzVector* lheZp4ptr_; TLorentzVector lheZp4_;

		// Ele-mu pair info
		TLorentzVector* eleMu_p4ptr_;
		TLorentzVector eleMu_p4_;
		Double_t eleMu_dR_;
		Double_t eleMu_dEta_;
		Double_t eleMu_dPhi_;

		// Electron info
		TLorentzVector* ele_p4ptr_;
		TLorentzVector ele_p4_;
		Int_t ele_codeHeepNoIso_;
		Int_t ele_charge_;

		// Muon info
		TLorentzVector* muon_p4ptr_;
		TLorentzVector muon_p4_;
		Int_t muon_codeHighPtId_;
		Int_t muon_charge_;



	public:
		// CTOR & DTOR //
		EleMuTree(const std::string& treeName, const std::string& fileName);
		~EleMuTree()  {}

		// PUBLIC METHODS //
		void fillTree(const tsw::EleMuObject& eleMu, const tsw::EventHelper& evtHelper, bool trigDecision);

	}; //End: EleMuTree

	EleMuTree::EleMuTree(const std::string& treeName, const std::string& fileName) :
		TreeHandlerBase(treeName, "Tree of ele-mu pairs ("+treeName+")", fileName),
		lheZp4ptr_( &lheZp4_ ),
		eleMu_p4ptr_( &eleMu_p4_ ),
		ele_p4ptr_( &ele_p4_ ),  muon_p4ptr_( &muon_p4_ )
	{
		// Setting up the event / di-ele branches ...
		mainAnaTree_->Branch("weight",     &evt_weight_,     "weight/D");
		mainAnaTree_->Branch("genWeight",  &evt_genWeight_,  "genWeight/D");
		mainAnaTree_->Branch("puWeight",   &evt_puWeight_,   "puWeight/D");
		mainAnaTree_->Branch("xsecWeight", &evt_xsecWeight_, "xsecWeight/D");

		mainAnaTree_->Branch("run",    &evt_runNum_,   "run/i");
		mainAnaTree_->Branch("lumi",   &evt_lumiNum_,  "lumi/i");
		mainAnaTree_->Branch("evtNum", &evt_evtNum_,   "evtNum/i");

		mainAnaTree_->Branch("trgDecision", &evt_trgDecision_, "trgDecision/O");
		mainAnaTree_->Branch("nVtx", &evt_nVtx_, "nVtx/i");

		if(isMadgraphDrellYan(fileName))
			mainAnaTree_->Branch("lheZp4", &lheZp4ptr_);

		mainAnaTree_->Branch("eleMu_p4", &eleMu_p4ptr_);
		mainAnaTree_->Branch("eleMu_dR", &eleMu_dR_, "eleMu_dR/D");
		mainAnaTree_->Branch("eleMu_dEta", &eleMu_dEta_, "eleMu_dEta/D");
		mainAnaTree_->Branch("eleMu_dPhi", &eleMu_dPhi_, "eleMu_dPhi/D");

		mainAnaTree_->Branch("ele_p4", &ele_p4_);
		mainAnaTree_->Branch("ele_heepNoIso", &ele_codeHeepNoIso_, "ele_heepNoIso/I" );
		mainAnaTree_->Branch("ele_charge", &ele_charge_, "ele_charge/I");

		mainAnaTree_->Branch("muon_p4", &muon_p4_);
		mainAnaTree_->Branch("muon_highPtId", &muon_codeHighPtId_, "muon_highPtId/I");
		mainAnaTree_->Branch("muon_charge", &muon_charge_, "muon_charge/I");
	}

	void EleMuTree::fillTree(const tsw::EleMuObject& eleMu, const tsw::EventHelper& evtHelper, bool trigDecision)
	{
		evt_weight_     = evtHelper.totWeight();
		evt_genWeight_  = evtHelper.genWeight();
		evt_puWeight_   = evtHelper.puWeight();
		evt_xsecWeight_ = evtHelper.xsecWeight();

		evt_runNum_ = evtHelper.runNum();
		evt_lumiNum_ = evtHelper.lumiSec();
		evt_evtNum_ = evtHelper.eventNum();

		evt_trgDecision_ = trigDecision;
		evt_nVtx_ = evtHelper.GetRecoVtx_nGoodVtxs();

		lheZp4_ = evtHelper.GetLHEZbosonP4();

		eleMu_p4_   = eleMu.p4();
		eleMu_dR_   = eleMu.deltaR();
		eleMu_dEta_ = eleMu.deltaEta();
		eleMu_dPhi_ = eleMu.deltaPhi();

		// Electron info
		ele_p4_            = eleMu.ele().p4();
		ele_codeHeepNoIso_ = eleMu.ele().heepIdNoIsoCutCode();
		ele_charge_        = eleMu.ele().charge();

		// Muon info
		muon_p4_           = eleMu.muon().p4();
		muon_codeHighPtId_ = eleMu.muon().highPtIdCutCode(35.0, 1.44);
		muon_charge_       = eleMu.muon().charge();

		// And finally fill the tree ...
		mainAnaTree_->Fill();
	}


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
		//Methods called by the constructor (e.g. Initialise methods)...
		void SetMemberVariablesToDefaultValues();
		void SetupBranchLinks();
	
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

//		tsw::ABCDMethodTree frPreDiEleTree_;
		tsw::DiEleTree zPrimeDiEleTree_;
		tsw::DiEleTree noIsoZCandDiEleTree_;
		tsw::DiEleTree modIsoZCandDiEleTree_;

		tsw::DiEleTree abcdDiGsfFrPreTree_;

		tsw::EffiCalcTree zCandEffiTree_;

		tsw::TagProbeTree heepTagProbeTree_;

		tsw::EleMuTree eleMuTree_fullId_;

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
	//
	zPrimeDiEleTree_("zPrimeDiEleTree", outputFileName_+"_zPrimeDiEleTree.root"),
	noIsoZCandDiEleTree_("noIsoZBosonTree", outputFileName_ + "_noIsoZCandTree.root", true),
	modIsoZCandDiEleTree_("modIsoZBosonTree", outputFileName_ + "_modIsoZCandTree.root"),
	abcdDiGsfFrPreTree_("abcdDiGsfFrPreTree", outputFileName_ + "_abcdDiGsfTree.root", true),
	zCandEffiTree_(outputFileName_ + "_zEffiTree.root"),
	heepTagProbeTree_("tagProbeTree", "qcdGsfGsfTree", outputFileName_+"_heepTpTree.root"),
	eleMuTree_fullId_("eleMuTree", outputFileName_ + "_eleMuTree.root")
{

	timer_DoAnalysis_readIn_.Stop(); timer_DoAnalysis_readIn_.Reset();
	timer_DoAnalysis_setup_.Stop(); timer_DoAnalysis_setup_.Reset();
	timer_DoAnalysis_FillHistograms_.Stop(); timer_DoAnalysis_FillHistograms_.Reset();
	timer_DoAnalysis_DoEMuMethod_.Stop(); timer_DoAnalysis_DoEMuMethod_.Reset();
	
	//Set default (typically clearly 'incorrect') values for non-histo variables...
	SetMemberVariablesToDefaultValues();

	// ----- SETTING UP TCHAIN AND ASSOCIATED BRANCH LINKS ----- //
	unsigned int inputTrees_numEvts = 0;
   std::vector<std::string>::const_iterator fileNameIt;

   TStopwatch nEntriesTimer; nEntriesTimer.Start();
   for( fileNameIt = inputFileNamesVec_.begin(); fileNameIt != inputFileNamesVec_.end(); fileNameIt++ ){
      inputFilesTChain_->Add(fileNameIt->c_str());
      TFile* theFile = TFile::Open(fileNameIt->c_str(), "READ");
      if( !theFile ){
      	std::cerr << " ERROR : Could not open file '" << (*fileNameIt) << "'" << std::endl;
      	exit( EXIT_FAILURE );
      }
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

		eventHelper_.SetEleStructValues(iEle, ithEleStruct);

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

void BstdZeeFirstAnalyser::FinishOffAnalysis()
{

	// Set the event counters for the tree handlers, and write trees to file ...
	zPrimeDiEleTree_.setEventCounter( this->GetNumEvtsRunOver() );
	zPrimeDiEleTree_.saveToFile();
	noIsoZCandDiEleTree_.setEventCounter( this->GetNumEvtsRunOver() );
	noIsoZCandDiEleTree_.saveToFile();
	modIsoZCandDiEleTree_.setEventCounter( this->GetNumEvtsRunOver() );
	modIsoZCandDiEleTree_.saveToFile();
	abcdDiGsfFrPreTree_.setEventCounter( this->GetNumEvtsRunOver() );
	abcdDiGsfFrPreTree_.saveToFile();
	if(isMC_){
		zCandEffiTree_.setEventCounter( this->GetNumEvtsRunOver() );
		zCandEffiTree_.saveToFile();
	}
	heepTagProbeTree_.setEventCounter( this->GetNumEvtsRunOver() );
	heepTagProbeTree_.saveToFile();
	eleMuTree_fullId_.setEventCounter( this->GetNumEvtsRunOver() );
	eleMuTree_fullId_.saveToFile();

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
	// ------------------------
	// PRELIMINARIES
	// (i.e. Checking which eles pass various cuts)

	std::vector<bool> normEles_HEEPCutsFlags;
	std::vector<bool> normEles_HEEPCutsNoIsoFlags;
	std::vector<bool> normEles_EB_HEEPCutsNoIsoFlags;
	std::vector<bool> normEles_EB_heepIdModIsoColThrFlags;
	std::vector<bool> normEles_EB_isFidEcalDrAndFRPreFlags;

	std::vector<bool> bstdEles_sCutsFlags;
	std::vector<bool> bstdEles_HEEPCutsFlags;
	std::vector<bool> bstdEles_HEEPCutsNoIsoFlags;


	for(std::vector<tsw::HEEPEle>::const_iterator eleIt = normEles_.begin(); eleIt != normEles_.end(); eleIt++){

		normEles_HEEPCutsFlags.push_back( eleIt->heepIdStdIsoCut(eventHelper_) );
		normEles_HEEPCutsNoIsoFlags.push_back( eleIt->heepIdNoIsoCut()  );
		normEles_EB_HEEPCutsNoIsoFlags.push_back( eleIt->heepIdNoIsoCut() && eleIt->isEB() );

		normEles_EB_heepIdModIsoColThrFlags.push_back( eleIt->heepIdModIsoColThrCut(eventHelper_) && eleIt->isEB() );

		normEles_EB_isFidEcalDrAndFRPreFlags.push_back( eleIt->fakeRatePreSelnCut() && eleIt->isHEEPEB() );
	}


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

	// ------------------------
	// HeepStdIso di-electrons for code sync with Z' ...
	if( tsw::NumPassingCuts(normEles_HEEPCutsFlags)>1 ){
		tsw::HEEPDiEle normDiEle_HeepStdIso(normEles_, normEles_HEEPCutsFlags);
		zPrimeDiEleTree_.fillTree(normDiEle_HeepStdIso, eventHelper_, trg_PathA_decision_);
	}


	// ------------------------
	// EB HEEPNoIso di-electrons...
	if(tsw::NumPassingCuts(normEles_EB_HEEPCutsNoIsoFlags)>1){
		tsw::HEEPDiEle normDiEle_EB_HEEPNoIso( normEles_, normEles_EB_HEEPCutsNoIsoFlags);

		// Check if the di-electron mass is in the Z mass range ...
//		zCandEffiTree_.FillTree(&normDiEle_EB_HEEPNoIso, mcZ_daughterA_p4_, mcZ_daughterB_p4_, evt_runNum_, evt_lumiSec_, evt_evtNum_, eventHelper_, evtWeight);
		if(normDiEle_EB_HEEPNoIso.isInZMassRange()){

			noIsoZCandDiEleTree_.fillTree(normDiEle_EB_HEEPNoIso, eventHelper_, trg_PathA_decision_);

		}
	}

	// ------------------------
	// EB Fr-preselection di-GSFs for ABCD method QCD estimate ...
	if(tsw::NumPassingCuts(normEles_EB_isFidEcalDrAndFRPreFlags)>1){
		tsw::HEEPDiEle frPreEbDiEle( normEles_, normEles_EB_isFidEcalDrAndFRPreFlags);

		if( frPreEbDiEle.invMass()>50.0 && frPreEbDiEle.invMass()<140.0 )
			abcdDiGsfFrPreTree_.fillTree(frPreEbDiEle, eventHelper_, trg_PathA_decision_);
	}

	// ------------------------
	// Fill MC-matched Z boson effi tree (iff running over MC)
	if(isMC_){
		const TLorentzVector& mcZ_eleA_p4 = (mcZ_numDaughters_>1 ? mcZ_daughterA_p4_ : mcEles_HighestEt_p4_);
		const TLorentzVector& mcZ_eleB_p4 = (mcZ_numDaughters_>1 ? mcZ_daughterB_p4_ : mcEles_2ndHighestEt_p4_);

		if(mcZ_eleA_p4.Pt()>0.00001 && mcZ_eleB_p4.Pt()>0.00001){
			if( tsw::HEEPDiEle* mcMatchedDiEle = getMcMatchedDiEle(mcZ_eleA_p4, mcZ_eleB_p4, normEles_) )
			{
				zCandEffiTree_.fillTree(*mcMatchedDiEle, mcZ_eleA_p4, mcZ_eleB_p4, eventHelper_);
				delete mcMatchedDiEle;
			}
			else
				zCandEffiTree_.fillTree_NotReconstructed(mcZ_eleA_p4, mcZ_eleB_p4, eventHelper_);
		}
	}

	// ------------------------
	// EB HEEPModIso di-electrons...
	if(tsw::NumPassingCuts(normEles_EB_heepIdModIsoColThrFlags)>1){
		tsw::HEEPDiEle heepIdModIsoEbDiEle( normEles_, normEles_EB_heepIdModIsoColThrFlags);

		if( heepIdModIsoEbDiEle.isInZMassRange() )
			modIsoZCandDiEleTree_.fillTree(heepIdModIsoEbDiEle, eventHelper_, trg_PathA_decision_);
	}

	//-----------------------------------
	// Selection for HEEP (std/mod iso) tag-probe studies:
	//     * Tag-Probe pairs ... TAG: barrel & pass HEEP ID  ;  PROBE: heepFid
	//     * GSF-GSF pairs for QCD estimation ... TAG CANDIDATE: barrel & heepFid  ;  PROBE: heepFid
	// [Since take all T-P/GSF-GSF pairs in each event satisfying these criteria, stricter requirements can be imposed later in T-P effi value calc'ing code.]
	for(std::vector<tsw::HEEPEle>::const_iterator tagEleIt = normEles_.begin(); tagEleIt != normEles_.end(); tagEleIt++){
		for(std::vector<tsw::HEEPEle>::const_iterator probeEleIt = normEles_.begin(); probeEleIt != normEles_.end(); probeEleIt++){
			if(tagEleIt==probeEleIt)
				continue;

			// Tag-probe pair selection
			if( tagEleIt->isHEEPEB() && tagEleIt->heepIdNoIsoCut() && probeEleIt->heepFidCut() )
				heepTagProbeTree_.fillTagProbeTree(*tagEleIt, *probeEleIt, eventHelper_, trg_PathA_decision_);

			// GSF-GSF pair selection
			if( tagEleIt->isHEEPEB() && tagEleIt->heepFidCut() && probeEleIt->heepFidCut() )
				heepTagProbeTree_.fillQcdTree(*tagEleIt, *probeEleIt, eventHelper_, trg_PathA_decision_);
		}
	}

	//-----------------------------------
	// Selection for ABCD tree


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
//---- Method for implementing emu analysis  ...          //
void BstdZeeFirstAnalyser::DoEMuAnalysis()
{
	// Grab highest Et electron passing HEEP V4.1 (w/o iso)
	const tsw::HEEPEle* ptr_ele_EB_HEEPNoIso_1stpT = 0;
	for(std::vector<tsw::HEEPEle>::const_iterator eleIt = normEles_.begin(); eleIt != normEles_.end(); eleIt++){
		if(!ptr_ele_EB_HEEPNoIso_1stpT){
			if( eleIt->heepIdNoIsoCut() )
				ptr_ele_EB_HEEPNoIso_1stpT = &(*eleIt);
		}
		else{
			if( eleIt->et()>ptr_ele_EB_HEEPNoIso_1stpT->et() && eleIt->heepIdNoIsoCut() )
				ptr_ele_EB_HEEPNoIso_1stpT = &(*eleIt);
		}
	}

	// Grab highest Pt muon passing highPt selection
	const tsw::Muon* ptr_muon_barrelHighPtId_1stpT = 0;
	for(std::vector<tsw::Muon>::const_iterator muonIt = normMuons_.vec().begin(); muonIt != normMuons_.vec().end(); muonIt++){
		if( !ptr_muon_barrelHighPtId_1stpT ){
			if( muonIt->highPtIdCut(35.0, 1.44) )
				ptr_muon_barrelHighPtId_1stpT = &(*muonIt);
		}
		else{
			if( muonIt->pT()>ptr_muon_barrelHighPtId_1stpT->pT() && muonIt->highPtIdCut(35.0, 1.44) )
				ptr_muon_barrelHighPtId_1stpT = &(*muonIt);
		}
	}

	// Fill eleMuTree, if ele-mu pair was found & has invariant mass in Z window
	if( !ptr_ele_EB_HEEPNoIso_1stpT || !ptr_muon_barrelHighPtId_1stpT )
		return;

	const tsw::EleMuObject eleMu_EB_HEEPNoIso_muB_highPtId(*ptr_ele_EB_HEEPNoIso_1stpT, *ptr_muon_barrelHighPtId_1stpT);
	if( eleMu_EB_HEEPNoIso_muB_highPtId.isInZMassRange() )
		eleMuTree_fullId_.fillTree( eleMu_EB_HEEPNoIso_muB_highPtId, eventHelper_, eventHelper_.GetTrigInfo_eMuPath_decision() );
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


