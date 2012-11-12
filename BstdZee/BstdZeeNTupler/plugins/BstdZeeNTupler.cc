// -*- C++ -*-
//
// Package:    BstdZeeNTupler
// Class:      BstdZeeNTupler
// 
/**\class BstdZeeNTupler BstdZeeNTupler.cc TSWilliams/BstdZeeNTupler/src/BstdZeeNTupler.cc

 Description: [one line class summary]
sjlkd

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams
//         Created:  Tue Apr 19 16:40:57 BST 2011
// $Id: BstdZeeNTupler.cc,v 1.22 2012/11/12 11:11:13 tsw Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Utilities/interface/Exception.h"

//To use the GenParticle MC-truth collection...
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//...for the GSF electron collection...
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

//...for the calo geometry classes ...
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

// ... for using the recHits collections and information ...
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"

// ... for the general tracks collection
#include "DataFormats/TrackReco/interface/Track.h"

//...in order to use the heep::Ele class ...
#include "SHarper/HEEPAnalyzer/interface/HEEPEle.h"

// ... for accessing the vertices & beamspot ...
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//...for accessing the muons ...
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

//...for accessing the HLT configuration & trigger decisions
#include <HLTrigger/HLTcore/interface/HLTConfigProvider.h>
#include <DataFormats/Common/interface/TriggerResults.h>

//...for accessing trigger objects
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//...for PU re-weighting
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

//... for accessing generator information (GenEventInfoProduct & LHE)
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//...for histograms creation
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

//...for NTuple generation
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

//Must include TROOT.h in order for the gROOT->ProcessLine(...) lines to work.
#include "TROOT.h"
#include "TSystem.h"
#include "TClassTable.h"

#include "TSWilliams/BstdZeeNTupler/interface/tswEvent.h"
#include "TSWilliams/BstdZeeNTupler/interface/tswUsefulFunctions.h"

//Functions...
Double_t CalcOpeningAngle(const ROOT::Math::XYZTVector& , const ROOT::Math::XYZTVector& );

//
// class declaration
//

class BstdZeeNTupler : public edm::EDAnalyzer {
   public:
      explicit BstdZeeNTupler(const edm::ParameterSet&);
      ~BstdZeeNTupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      std::pair<std::string, unsigned int> beginRun_getTriggerNameAndIdx(const std::vector<std::string>& vecOfPossNames, std::string triggerId_debug);

      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		//For setting up the pointer links between the branches and the variables...
		void SetupMCBranches();
		void SetupTriggerBranches();
		void SetupEvtInfoBranches();
		void SetupStdEleBranches();
		void SetupBstdEleBranches();

		//For clearing contents/setting default values of variables that should get new values in each event...
		void ResetEventByEventVariables();

		//For accessing the trigger info...
		void accessTriggerInfo(const edm::Event&, const edm::EventSetup&, bool beVerbose);
		void printDatasetsAndTriggerNames();

		//For reading the event information into the member variables...
		void ReadInEvtInfo(bool, const edm::Event&);
		//For reading vertex info. and doing pile up re-weighting ...
		void ReadInVtxAndPUInfo(bool, const edm::Event&);

		// For reading generator info into tsw::Event member var ...
		void ReadInGenInfo(bool, const edm::Event&);


		//For reading the electron information into the member variables...
		void ReadInNormGsfEles(bool, const edm::Handle<reco::GsfElectronCollection>&, edm::ESHandle<CaloGeometry>&, EcalRecHitMetaCollection*, EcalRecHitMetaCollection*, edm::Handle<reco::TrackCollection>&, const edm::Event&);
		void ReadInBstdGsfEles(bool, const edm::Handle<reco::GsfElectronCollection>&);

		//For reading in the muon information, and dumping it into tsw::Event class ...
		void ReadInMuons(bool, const edm::Event&);

		//For calculating kinematic parameters...
		void SetMCZcandidateVariables(bool);

      // ----------member data ---------------------------
		//TH1D *numZs_status2Hist;
		//TH1D *numZs_status3Hist;
		//TH1D *numZs_statusOtherHist;
		//TH1D *status2Zs_PtHist;
		//TH1D *status2Zs_LogPtHist;
		TH1D *eleHighestEt_EtHist;
		TH1D *ele2ndHighestEt_EtHist;
		TH1D *fracAboveSingleEleEtThreshold_Hist;
		edm::Service<TFileService> fHistos;
		TTree* EventDataTree;

		int numEvts;
		unsigned int numEvtsStored_;
		// Member vars read in from cfg file
		int dyJetsToLL_EventType_;
		bool mcFlag_;
		bool vBool_;
		const bool readInNormReco_;
		const bool readInBstdReco_;
		const bool readInTrigInfo_;
		const bool is2010SignalDataset_;
		const bool useReducedRecHitsCollns_;
		const edm::InputTag normEleInputTag_;
		const edm::InputTag vertexSrc_;

		double histoEtMin;
		double histoEtMax;
		int    histoEtNBins;

		edm::LumiReWeighting lumiWeights_;
		const bool calcPileUpWeights_;
		reco::Vertex recoVtx_highestSumPtVtx_;

		//Variables whose values will be stored as branches...
		tsw::Event* event_;

		unsigned int numFinalStateEles;

		Int_t mcEles_HighestEt_charge;
		Int_t mcEles_HighestEt_PDGid;
		Int_t mcEles_HighestEt_status;
		Double_t mcEles_HighestEt_pt;
		ROOT::Math::XYZTVector mcEles_HighestEt_p4;
		ROOT::Math::XYZTVector* mcEles_HighestEt_p4_ptr;

		Int_t mcEles_2ndHighestEt_charge;
		Int_t mcEles_2ndHighestEt_PDGid;
		Int_t mcEles_2ndHighestEt_status;
		Double_t mcEles_2ndHighestEt_pt;
		ROOT::Math::XYZTVector mcEles_2ndHighestEt_p4;
		ROOT::Math::XYZTVector* mcEles_2ndHighestEt_p4_ptr;

		Double_t mcZcandidate_pt;
		Double_t mcZcandidate_eta;
		Double_t mcZcandidate_phi;
		Double_t mcZcandidate_mass;
		ROOT::Math::XYZTVector mcZcandidate_p4;
		ROOT::Math::XYZTVector* mcZcandidate_p4_ptr;
		Double_t mcZcandidate_dEtaEles;
		Double_t mcZcandidate_dPhiEles;
		Double_t mcZcandidate_dREles;
		Double_t mcZcandidate_openingAngle;

		Int_t mc_numZs_status2_;
		Int_t mc_numZs_status3_;
		Int_t mc_numZs_statusNot2Or3_;

		Int_t mcZboson_pdgId_; // Int_t
		Int_t mcZboson_status_; // Int_t
		ROOT::Math::XYZTVector* mcZboson_p4_ptr_; // ROOT::Math::XYZTLorentzVector*
		ROOT::Math::XYZTVector mcZboson_p4_;
		ROOT::Math::XYZTVector* mcZboson_status2_p4_ptr_; // ROOT::Math::XYZTLorentzVector*
		ROOT::Math::XYZTVector mcZboson_status2_p4_;
		Int_t mcZboson_numDaughters_;
		Double_t mcZboson_daughters_dR_;
		Double_t mcZboson_daughters_dEta_;
		Double_t mcZboson_daughters_dPhi_;
		Double_t mcZboson_daughters_openingAngle_;
		ROOT::Math::XYZTVector  mcZboson_daughterA_p4_;
		ROOT::Math::XYZTVector* mcZboson_daughterA_p4_ptr_;
		ROOT::Math::XYZTVector  mcZboson_daughterB_p4_;
		ROOT::Math::XYZTVector* mcZboson_daughterB_p4_ptr_;

		// Trigger member variables
		Bool_t      trg_PathA_decision_;
		std::string trg_PathA_name_;
		
		edm::InputTag hltResultsTag_;
		edm::InputTag hltEventTag_;
		std::string hltPathA_;
		const std::vector<std::string> hltPathA_possNames_;
		bool hltPathADecision_;
		std::string hltPathA_nameOfLastFilter_;
		float hltPathA_highestTrigObjEt_;
		TH1D* highestEtTrigObjPathA_EtHist_;
		TH1D* highestEtTrigObjPathA_cumEtHist_;
		TH1D* highestEtTrigObjPathA_EtThrDenomHist_;
		TH1D* highestEtTrigObjPathA_EtThrEffiHist_;
		TGraphAsymmErrors* highestEtTrigObjPathA_EtThrAsymmEffi_;
		//Trigger member variables - HLT config helper
		HLTConfigProvider hltConfig_;
  		unsigned hltPathIndex_;
		 // Trigger member variables - emu cross trigger
		const std::vector<std::string> trg_emuPath_possNames_;
		std::string trg_emuPath_name_;
		unsigned int trg_emuPath_idx_;

		//Event information variables...
		unsigned int evt_runNum_;
		unsigned int evt_lumiSec_;
		unsigned int evt_evtNum_;

		//Standard GSF electron variables...
		unsigned int normGsfEles_number_;
		std::vector<ROOT::Math::XYZTVector> normGsfEles_p4_;
		std::vector<ROOT::Math::XYZTVector>* normGsfEles_p4ptr_;
		std::vector<Int_t> normGsfEles_charge_;

		std::vector<Double_t> normGsfEles_Et_;
		std::vector<Double_t> normGsfEles_HEEP_Et_;
		std::vector<Double_t> normGsfEles_Eta_;
		std::vector<Double_t> normGsfEles_scEta_;
		std::vector<bool> normGsfEles_ecalDriven_;
		std::vector<bool> normGsfEles_ecalDrivenSeed_;

		std::vector<Double_t> normGsfEles_HEEP_dEtaIn_;
		std::vector<Double_t> normGsfEles_HEEP_dPhiIn_;
		std::vector<Double_t> normGsfEles_HoverE_;
		std::vector<Double_t> normGsfEles_sigmaIetaIeta_;
		std::vector<Double_t> normGsfEles_scSigmaIetaIeta_;
		
		std::vector<Double_t> normGsfEles_dr03EmIsoEt_;
		std::vector<Double_t> normGsfEles_dr03HadDepth1IsoEt_;
		std::vector<Double_t> normGsfEles_dr03HadDepth2IsoEt_;
		std::vector<Double_t> normGsfEles_dr03TkIsoPt_;

		std::vector<Double_t> normGsfEles_e2x5Max_;
		std::vector<Double_t> normGsfEles_e5x5_;

		//kinematic and geometric methods
		std::vector<float> normHEEPEles_et_;
		std::vector<float> normHEEPEles_gsfEt_;
		std::vector<float> normHEEPEles_scEt_;
		std::vector<float> normHEEPEles_energy_;
		std::vector<float> normHEEPEles_gsfEnergy_;
		std::vector<float> normHEEPEles_caloEnergy_;
		std::vector<float> normHEEPEles_ecalEnergyError_;
		std::vector<float> normHEEPEles_eta_;
		std::vector<float> normHEEPEles_scEta_;
		std::vector<float> normHEEPEles_detEta_;
		std::vector<float> normHEEPEles_detEtaAbs_;
		std::vector<float> normHEEPEles_phi_;
		std::vector<float> normHEEPEles_scPhi_;
		std::vector<float> normHEEPEles_detPhi_;
		std::vector<float> normHEEPEles_zVtx_;
		std::vector<ROOT::Math::XYZTVector> normHEEPEles_p4_;
		std::vector<ROOT::Math::XYZTVector>* normHEEPEles_p4ptr_;
		std::vector<ROOT::Math::XYZTVector> normHEEPEles_gsfP4_;
		std::vector<ROOT::Math::XYZTVector>* normHEEPEles_gsfP4ptr_;

		//classification (couldnt they have just named it 'type')
		std::vector<Int_t>  normHEEPEles_classification_;
		std::vector<Bool_t> normHEEPEles_isEcalDriven_;
		std::vector<Bool_t> normHEEPEles_isTrackerDriven_;
		std::vector<Bool_t> normHEEPEles_isEB_;
		std::vector<Bool_t> normHEEPEles_isEE_;

		//track methods
		std::vector<Int_t> normHEEPEles_charge_;
		std::vector<Int_t> normHEEPEles_trkCharge_;
		std::vector<float> normHEEPEles_pVtx_;
		std::vector<float> normHEEPEles_pCalo_;
		std::vector<float> normHEEPEles_ptVtx_;
		std::vector<float> normHEEPEles_ptCalo_;
		std::vector<float> normHEEPEles_closestCtfTrk_pt_;
		std::vector<float> normHEEPEles_closestCtfTrk_eta_;
		std::vector<float> normHEEPEles_closestCtfTrk_phi_;
		std::vector<float> normHEEPEles_closestCtfTrk_innerPt_;
		std::vector<float> normHEEPEles_closestCtfTrk_innerEta_;
		std::vector<float> normHEEPEles_closestCtfTrk_innerPhi_;
		std::vector<float> normHEEPEles_closestCtfTrk_outerPt_;
		std::vector<float> normHEEPEles_closestCtfTrk_outerEta_;
		std::vector<float> normHEEPEles_closestCtfTrk_outerPhi_;

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		std::vector<float> normHEEPEles_hOverE_;
		std::vector<float> normHEEPEles_dEtaIn_;
		std::vector<float> normHEEPEles_dPhiIn_;
		std::vector<float> normHEEPEles_dPhiOut_;
		std::vector<float> normHEEPEles_epIn_;
		std::vector<float> normHEEPEles_epOut_;
		std::vector<float> normHEEPEles_fbrem_;
		std::vector<float> normHEEPEles_bremFrac_;
		std::vector<float> normHEEPEles_invEOverInvP_;

		//shower shape variables
		std::vector<float> normHEEPEles_sigmaEtaEta_;
		std::vector<float> normHEEPEles_sigmaEtaEtaUnCorr_;
		std::vector<float> normHEEPEles_sigmaIEtaIEta_;
		std::vector<float> normHEEPEles_e1x5_;
		std::vector<float> normHEEPEles_e2x5Max_;
		std::vector<float> normHEEPEles_e5x5_;
		std::vector<float> normHEEPEles_e1x5Over5x5_;
		std::vector<float> normHEEPEles_e2x5MaxOver5x5_;

		//isolation, we use cone of 0.3
		std::vector<float> normHEEPEles_isolEm_;
		std::vector<float> normHEEPEles_isolHad_;
		std::vector<float> normHEEPEles_isolHadDepth1_;
		std::vector<float> normHEEPEles_isolHadDepth2_;
		std::vector<float> normHEEPEles_isolPtTrks_;
		std::vector<float> normHEEPEles_isolEmHadDepth1_;

		std::vector<UInt_t> normHEEPEles_numMissInnerHits_;

		std::vector<float> normHEEPEles_SCposn_eta_;
		std::vector<float> normHEEPEles_SCposn_phi_;
		std::vector<float> normHEEPEles_SC_rawEnergy_;
		std::vector< std::vector<float> > normHEEPEles_SC_recHits_Et_;
		std::vector< std::vector<float> > normHEEPEles_SC_recHits_eta_;
		std::vector< std::vector<float> > normHEEPEles_SC_recHits_phi_;
		std::vector< std::vector<bool> >  normHEEPEles_SC_recHits_isFromEB_;
		std::vector<float> normHEEPEles_SC_totEnergyRecHits_;
		std::vector<unsigned int> normHEEPEles_SC_totNumRecHits_;

		std::vector<float> normHEEPEles_gsfTrk_eta_;
		std::vector<float> normHEEPEles_gsfTrk_phi_;
		std::vector<float> normHEEPEles_gsfTrk_vz_;
		std::vector< std::vector<float> > normHEEPEles_innerIsoConeTrks_pt_;
		std::vector< std::vector<float> > normHEEPEles_innerIsoConeTrks_eta_;
		std::vector< std::vector<float> > normHEEPEles_innerIsoConeTrks_phi_;
		std::vector< std::vector<float> > normHEEPEles_innerIsoConeTrks_vz_;

		//Boosted reco'n electron variables...
		unsigned int bstdGsfEles_number_;
		std::vector<ROOT::Math::XYZTVector> bstdGsfEles_p4_;
		std::vector<ROOT::Math::XYZTVector>* bstdGsfEles_p4ptr_;
		std::vector<Int_t> bstdGsfEles_charge_;

		std::vector<Double_t> bstdGsfEles_Et_;
		std::vector<Double_t> bstdGsfEles_HEEP_Et_;
		std::vector<Double_t> bstdGsfEles_Eta_;
		std::vector<Double_t> bstdGsfEles_scEta_;
		std::vector<bool> bstdGsfEles_ecalDriven_;
		std::vector<bool> bstdGsfEles_ecalDrivenSeed_;

		std::vector<Double_t> bstdGsfEles_HEEP_dEtaIn_;
		std::vector<Double_t> bstdGsfEles_HEEP_dPhiIn_;
		std::vector<Double_t> bstdGsfEles_HoverE_;
		std::vector<Double_t> bstdGsfEles_sigmaIetaIeta_;
		std::vector<Double_t> bstdGsfEles_scSigmaIetaIeta_;
		
		std::vector<Double_t> bstdGsfEles_dr03EmIsoEt_;
		std::vector<Double_t> bstdGsfEles_dr03HadDepth1IsoEt_;
		std::vector<Double_t> bstdGsfEles_dr03HadDepth2IsoEt_;
		std::vector<Double_t> bstdGsfEles_dr03TkIsoPt_;

		std::vector<Double_t> bstdGsfEles_e2x5Max_;
		std::vector<Double_t> bstdGsfEles_e5x5_;
		
		//kinematic and geometric methods
		std::vector<float> bstdHEEPEles_et_;
		std::vector<float> bstdHEEPEles_gsfEt_;
		std::vector<float> bstdHEEPEles_scEt_;
		std::vector<float> bstdHEEPEles_energy_;
		std::vector<float> bstdHEEPEles_gsfEnergy_;
		std::vector<float> bstdHEEPEles_caloEnergy_;
		std::vector<float> bstdHEEPEles_eta_;
		std::vector<float> bstdHEEPEles_scEta_;
		std::vector<float> bstdHEEPEles_detEta_;
		std::vector<float> bstdHEEPEles_detEtaAbs_;
		std::vector<float> bstdHEEPEles_phi_;
		std::vector<float> bstdHEEPEles_scPhi_;
		std::vector<float> bstdHEEPEles_detPhi_;
		std::vector<float> bstdHEEPEles_zVtx_;
		std::vector<ROOT::Math::XYZTVector> bstdHEEPEles_p4_;
		std::vector<ROOT::Math::XYZTVector>* bstdHEEPEles_p4ptr_;
		std::vector<ROOT::Math::XYZTVector> bstdHEEPEles_gsfP4_;
		std::vector<ROOT::Math::XYZTVector>* bstdHEEPEles_gsfP4ptr_;

		//classification (couldnt they have just named it 'type')
		std::vector<Int_t>  bstdHEEPEles_classification_;
		std::vector<Bool_t> bstdHEEPEles_isEcalDriven_;
		std::vector<Bool_t> bstdHEEPEles_isTrackerDriven_;
		std::vector<Bool_t> bstdHEEPEles_isEB_;
		std::vector<Bool_t> bstdHEEPEles_isEE_;

		//track methods
		std::vector<Int_t> bstdHEEPEles_charge_;
		std::vector<Int_t> bstdHEEPEles_trkCharge_;
		std::vector<float> bstdHEEPEles_pVtx_;
		std::vector<float> bstdHEEPEles_pCalo_;
		std::vector<float> bstdHEEPEles_ptVtx_;
		std::vector<float> bstdHEEPEles_ptCalo_;

		//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
		std::vector<float> bstdHEEPEles_hOverE_;
		std::vector<float> bstdHEEPEles_dEtaIn_;
		std::vector<float> bstdHEEPEles_dPhiIn_;
		std::vector<float> bstdHEEPEles_dPhiOut_;
		std::vector<float> bstdHEEPEles_epIn_;
		std::vector<float> bstdHEEPEles_epOut_;
		std::vector<float> bstdHEEPEles_fbrem_;
		std::vector<float> bstdHEEPEles_bremFrac_;
		std::vector<float> bstdHEEPEles_invEOverInvP_;

		//shower shape variables
		std::vector<float> bstdHEEPEles_sigmaEtaEta_;
		std::vector<float> bstdHEEPEles_sigmaEtaEtaUnCorr_;
		std::vector<float> bstdHEEPEles_sigmaIEtaIEta_;
		std::vector<float> bstdHEEPEles_e1x5_;
		std::vector<float> bstdHEEPEles_e2x5Max_;
		std::vector<float> bstdHEEPEles_e5x5_;
		std::vector<float> bstdHEEPEles_e1x5Over5x5_;
		std::vector<float> bstdHEEPEles_e2x5MaxOver5x5_;

		//isolation, we use cone of 0.3
		std::vector<float> bstdHEEPEles_isolEm_;
		std::vector<float> bstdHEEPEles_isolHad_;
		std::vector<float> bstdHEEPEles_isolHadDepth1_;
		std::vector<float> bstdHEEPEles_isolHadDepth2_;
		std::vector<float> bstdHEEPEles_isolPtTrks_;
		std::vector<float> bstdHEEPEles_isolEmHadDepth1_;
};

//
// constants, enums and typedefs

//
// static data member definitions
//

//
// constructors and destructor
//
BstdZeeNTupler::BstdZeeNTupler(const edm::ParameterSet& iConfig):
	numEvts(0),
	numEvtsStored_(0),
	dyJetsToLL_EventType_(iConfig.getUntrackedParameter<int>("dyJetsToLL_EventType",0)), //==0=>Don't select events, ==11=>ele, ==13=>muon, ==15=>tau
	mcFlag_(iConfig.getUntrackedParameter<bool>("isMC",0)),
	vBool_(iConfig.getUntrackedParameter<bool>("printOutInfo",0)),
	readInNormReco_(iConfig.getUntrackedParameter<bool>("readInNormReco",0)),
	readInBstdReco_(iConfig.getUntrackedParameter<bool>("readInBstdReco",0)),
	readInTrigInfo_(iConfig.getUntrackedParameter<bool>("readInTrigInfo",0)),
	is2010SignalDataset_(iConfig.getUntrackedParameter<bool>("is2010SignalDataset",0)),
	useReducedRecHitsCollns_(iConfig.getUntrackedParameter<bool>("useReducedRecHitsCollns",0)),
	normEleInputTag_(iConfig.getUntrackedParameter<edm::InputTag>("eleCollection") ),
	vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")),
	histoEtMin(0.0),
	histoEtMax(80.0),
	histoEtNBins(40),
	calcPileUpWeights_(  (iConfig.exists("puDists_mcFile") && iConfig.exists("puDists_dataFile")) &&
								(iConfig.exists("puDists_mcHistName") && iConfig.exists("puDists_dataHistName"))  ),
	hltResultsTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltResultsTag",edm::InputTag("TriggerResults","","HLT"))),
	hltEventTag_(iConfig.getUntrackedParameter<edm::InputTag>("hltEventTag",edm::InputTag("hltTriggerSummaryAOD","","HLT"))),
  	hltPathA_("*** DEFAULT TRIGGER NAME ***"),
  	hltPathA_possNames_(iConfig.getUntrackedParameter< std::vector<std::string> >("hltPathA_possNames")),
	hltPathADecision_(false),
  	hltPathA_nameOfLastFilter_(iConfig.getUntrackedParameter<std::string>("hltPathA_nameOfLastFilter",std::string("hltEle45CaloIdVTTrkIdTDphiFilter"))),
  	hltPathA_highestTrigObjEt_(-999.9),
	trg_emuPath_possNames_( iConfig.getUntrackedParameter< std::vector<std::string> >("trg_emuPath_possNames") ),
	trg_emuPath_name_("*** DEFAULT TRIGGER NAME ***"),
	trg_emuPath_idx_(999999)

{
   //now do what ever initialization is needed
	std::cout << " * PU re-weighting ... "  << std::endl;
	if(calcPileUpWeights_){
		std::cout << "         -> Instantiating edm::LumiReweigting with specified file/hist names ..." << std::endl
					 << "            ( Seg fault will occur if file/hist names are incorrect! )" << std::endl;
		lumiWeights_ = edm::LumiReWeighting( iConfig.getUntrackedParameter<std::string>("puDists_mcFile"),
							iConfig.getUntrackedParameter<std::string>("puDists_dataFile"),
							iConfig.getUntrackedParameter<std::string>("puDists_mcHistName"),
							iConfig.getUntrackedParameter<std::string>("puDists_dataHistName") );
		std::cout << "         -> DONE !" << std::endl;
	}
	else
		std::cout << "   *** WARNING : PU weights will not be calculated (relevant PU dist'n file & hist names not specified in config) !" << std::endl;


	double histoEtBinWidth = (histoEtMax-histoEtMin)/static_cast<double>(histoEtNBins);
	
	eleHighestEt_EtHist    = fHistos->make<TH1D>("eleHighestEt_EtHist",   "p_{T} distribution of highest p_{T} electron; Electron p_{T} /GeVc^{-1}; Number per 2GeVc^{-1}",    histoEtNBins, histoEtMin, histoEtMax);
	ele2ndHighestEt_EtHist = fHistos->make<TH1D>("ele2ndHighestEt_EtHist","p_{T} distribution of 2nd highest p_{T} electron; Electron p_{T} /GeVc^{-1}; Number per 2GeVc^{-1}",histoEtNBins, histoEtMin, histoEtMax);
	fracAboveSingleEleEtThreshold_Hist = fHistos->make<TH1D>("fracAboveSingleEleEtThreshold_Hist", "Fraction of events passing a single electron Et threshold; Threshold E_{T} /GeV; Fraction", histoEtNBins, histoEtMin-histoEtBinWidth/2.0, histoEtMax-histoEtBinWidth/2.0);
	Int_t trigEtHist_nBins=80; Double_t trigEtHist_min = 0.0; Double_t trigEtHist_max = 200.0; Double_t trigEtHist_binWidth = (trigEtHist_max-trigEtHist_min)/trigEtHist_nBins;
	highestEtTrigObjPathA_EtHist_ = fHistos->make<TH1D>("highestEtTrigObjPathA_EtHist", "E_{T} distribution for the highest E_{T} HLT electron passing the last filter in hltPathA; HLT electron E_{T} /GeV; Number of electrons per GeV", trigEtHist_nBins, trigEtHist_min, trigEtHist_max);
	highestEtTrigObjPathA_cumEtHist_ = fHistos->make<TH1D>("highestEtTrigObjPathA_cumEtHist", "Cumulative E_{T} distribution for the highest E_{T} HLT electron passing the last filter in hltPathA; HLT electron E_{T} threshold, E_{T}^{thr} /GeV; Number of electrons with E_{T} > E_{T}^{thr}", trigEtHist_nBins, trigEtHist_min - trigEtHist_binWidth/2.0, trigEtHist_max - trigEtHist_binWidth/2.0);
	highestEtTrigObjPathA_EtThrDenomHist_ = fHistos->make<TH1D>("highestEtTrigObjPathA_EtThrDenomHist", "Denominator histogram for E_{T} threshold trigger efficiency curve; HLT electron E_{T} threshold, E_{T}^{thr} /GeV; Total number of events run over",                                      trigEtHist_nBins, trigEtHist_min - trigEtHist_binWidth/2.0, trigEtHist_max - trigEtHist_binWidth/2.0);
	highestEtTrigObjPathA_EtThrEffiHist_  = fHistos->make<TH1D>("highestEtTrigObjPathA_EtThrEffiHist", "Trigger efficiency for signal events vs. signal trigger E_{T} threshold; HLT electron E_{T} threshold, E_{T}^{thr} /GeV; Efficiency for signal events",                                trigEtHist_nBins, trigEtHist_min - trigEtHist_binWidth/2.0, trigEtHist_max - trigEtHist_binWidth/2.0);
	//For TTree and TBranch creation...
	EventDataTree = fHistos->make<TTree>("EventDataTree", "Electron event data");
	//EventDataTree = new TTree("EventDataTree", "GSFelectroneventdata");
	
	//Setting up the links between variables and branches...
	event_ = 0;
	EventDataTree->Branch("event","tsw::Event", &event_, 64000, 1); // This line was taken from Jim's tupiliser
	SetupEvtInfoBranches();
	SetupTriggerBranches();
	if(mcFlag_)
		SetupMCBranches();
	SetupStdEleBranches();
	SetupBstdEleBranches();

}


BstdZeeNTupler::~BstdZeeNTupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BstdZeeNTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//Variable declarations...
	//int numZs_status2 = 0;
	//int numZs_status3 = 0;
	//int numZs_statusOther = 0;
	numFinalStateEles = 0;
	int particlePDGid = -999;
	int particleStatus = -999;
	reco::Candidate::LorentzVector particle4mom(0.0, 0.0, 0.0, 0.0);
	unsigned int numDaughters = -999;
	//const reco::GenParticle mcParticle;
	int idx_eleHighestEt = -999;
	double eleHighestEt = -999.9; //Must be set to less than 0.0
	int idx_ele2ndHighestEt = -999; 
	double ele2ndHighestEt = -9999.9; //Should be set to less than eleHighestEt 
	//const reco::GenParticle& mcEle_HighestEt;
	//const reco::GenParticle& mcEle_2ndHighestEt;

	event_ = new tsw::Event();

   using namespace edm;
	
	//Clearing contents/setting default values of variables that should get new values in each event...
	ResetEventByEventVariables();

	//Handles, and grabbing various CMSSW-specific data types from event...
	Handle<reco::GenParticleCollection> genParticles;
	edm::Handle<reco::GsfElectronCollection> normGsfElesH, bstdGsfElesH;
	edm::ESHandle<CaloGeometry> caloGeomH;
	edm::Handle<EcalRecHitCollection> reducedEBRecHitsH, reducedEERecHitsH;
	edm::Handle< reco::TrackCollection > ctfTrackCollectionH;

	if(mcFlag_)
		iEvent.getByLabel("genParticles", genParticles);
	
	//Getting the handle for the standard and special reco'n GSF electron data...
	if(readInNormReco_)
		iEvent.getByLabel(normEleInputTag_, normGsfElesH);
	//iEvent.getByLabel(edm::InputTag("gsfElectrons::BOOSTEDRECO"), bstdGsfElesH);
	if(readInBstdReco_)
		iEvent.getByLabel(edm::InputTag("ecalDrivenGsfElectrons"), bstdGsfElesH);

	// Setting the handle for the calo geometry information ...
	iSetup.get<CaloGeometryRecord>().get(caloGeomH);

	// Setting the handle for the EB and EE recHits collections ...
	if(useReducedRecHitsCollns_){
		iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEB::RECO"), reducedEBRecHitsH);
		iEvent.getByLabel(edm::InputTag("reducedEcalRecHitsEE::RECO"), reducedEERecHitsH);
	}
	else{
		iEvent.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEB:RECO"), reducedEBRecHitsH);
		iEvent.getByLabel(edm::InputTag("ecalRecHit:EcalRecHitsEE:RECO"), reducedEERecHitsH);
	}
	EcalRecHitMetaCollection* ecalBarrelHitsMeta = new EcalRecHitMetaCollection(*reducedEBRecHitsH);
	EcalRecHitMetaCollection* ecalEndcapHitsMeta = new EcalRecHitMetaCollection(*reducedEERecHitsH);

	// Setting the handle for the general CTF track collection ...
	iEvent.getByLabel(edm::InputTag("generalTracks::RECO"), ctfTrackCollectionH);

	//Reading in the event information...
	ReadInEvtInfo(vBool_, iEvent);
	
	// Reading in the generator (i.e. GenEventInfoProduct & LHE) info, if MC ...
	if(mcFlag_)
		ReadInGenInfo(vBool_, iEvent);

	//Getting the trigger information for the event...	
	if(readInTrigInfo_)
		accessTriggerInfo(iEvent,iSetup, vBool_); //This sets the variable hltPathADecision_ as the event pass/fail decision for the HLT path hltPathA_ in the trigger product specified by HLTResultsTag_
	//Storing branches of trigger information...
	trg_PathA_decision_ = hltPathADecision_;
	trg_PathA_name_ = hltPathA_;
	if(vBool_){
		if(hltPathADecision_)
			std::cout << " ->Event passed trigger path A (i.e. " << trg_PathA_name_ << "=" << hltPathA_ << "?)" << std::endl;
		else
			std::cout << " ->Event did *NOT* pass trigger path A (i.e. " << trg_PathA_name_ << "=" << hltPathA_ << "?)" << std::endl;
	}
	//SetTriggerBranchVariables(true);
		
	//std::cout << "Analysing event no. " << iEvent.id() << " ..." << std::endl;

	numEvts++;

	//Reading in the kin variables for the standard and special reco'n GSF electrons... 
	if(readInNormReco_)
		ReadInNormGsfEles(vBool_, normGsfElesH, caloGeomH, ecalBarrelHitsMeta, ecalEndcapHitsMeta, ctfTrackCollectionH, iEvent);
	if(readInBstdReco_)
		ReadInBstdGsfEles(vBool_, bstdGsfElesH);

	delete ecalBarrelHitsMeta;
	delete ecalEndcapHitsMeta;

	ReadInMuons(vBool_, iEvent);

	Int_t zBosonDaughter_PDGid = 0;
	Int_t zBosonDaughterIfSherpa_PDGid = 0;
	if(mcFlag_){
		if(vBool_){std::cout << " ->There are " << genParticles->size() << " generator-particles in this event." << std::endl;}

		//Running over all of the particles in the MC truth...
		for(unsigned int idx=0; idx<genParticles->size(); idx++){
			//Extracting the properties of each of the MC truth particles....
			const reco::GenParticle & mcParticle = (*genParticles)[idx];
			particlePDGid = mcParticle.pdgId();
			particleStatus = mcParticle.status();
			particle4mom = mcParticle.p4();
			numDaughters = mcParticle.numberOfDaughters();
			if(particlePDGid==23 && particleStatus==3 && numDaughters!=0){ // i.e. if it is a 'decaying' Z boson
				//Get the first Z boson daughter ...
				const reco::Candidate* mcDaughter = mcParticle.daughter(0);
				zBosonDaughter_PDGid = mcDaughter->pdgId();
			}
			if(particleStatus==3 && (abs(particlePDGid)==11 || abs(particlePDGid)==13 || abs(particlePDGid)==15) && mcParticle.numberOfMothers()==2 ){
				zBosonDaughterIfSherpa_PDGid = particlePDGid;
			}
			//Finding the Z bosons, and counting up numbers with various statuses...
			if(abs(particlePDGid)==23){
				if(vBool_){
					std::cout << "   Generator-particle no. " << idx << " is a Z boson. (pdgId= " << particlePDGid << "; status=" << particleStatus << ")" << std::endl;
					std::cout << "      px=" << particle4mom.Px() << "GeV?; py=" << particle4mom.Py() << "GeV?; pz=" << particle4mom.Pz() << "GeV?" << std::endl;
					std::cout << "      pt=" << particle4mom.Pt() << "GeV?; eta=" << particle4mom.Eta() << "; phi=" << particle4mom.Phi() << std::endl;
					std::cout << "      It has " << numDaughters << " daughters." << std::endl;
				}
				for(size_t j=0; j<numDaughters; j++){
					const reco::Candidate * mcDaughter = mcParticle.daughter(j);
					if(vBool_){ std::cout << "         Daughter " << j << ": pdgId=" << mcDaughter->pdgId() << "; (pT, eta, phi) = (" << mcDaughter->p4().Pt() << "," << mcDaughter->p4().Eta() << "," << mcDaughter->p4().Phi() << ")" << std::endl; }
				}

				//For Z bosons with status=2 (i.e. that are decaying)
				if(particleStatus==2){
					mc_numZs_status2_++;
					mcZboson_status2_p4_ = particle4mom;
				}
				//For Z bosons with status=3 (i.e. that are in the 'hard part' of the interaction - i.e. used in the matrix element caLculation)
				else if(particleStatus==3){
					mc_numZs_status3_++;
					mcZboson_pdgId_  = particlePDGid; // Int_t
					mcZboson_status_ = particleStatus; // Int_t
					mcZboson_p4_     = particle4mom;

					if(vBool_){
						std::cout << "      Status=3=> written out to branches: pdgId=" << mcZboson_pdgId_ << "; status=" << mcZboson_status_ << "; ";
						std::cout << "(pT, eta, phi)=(" << mcZboson_p4_.Pt() << ", " << mcZboson_p4_.Eta() << ", " << mcZboson_p4_.Phi() << ")" << std::endl;
					}
					mcZboson_numDaughters_ = numDaughters;
					if(mcZboson_numDaughters_>=2){
						Int_t idx_daughterA=-999; Int_t idx_daughterB=-999;
						// If statements to account for the fact that the 'status=3' Z boson always seems to have 3 daughters - the two particles that it decays into, and a 'status=2' Z boson ...
						if(mcParticle.daughter(0)->pdgId()==23){
							idx_daughterA = 1;
							idx_daughterB = 2;
						}
						else if (mcParticle.daughter(1)->pdgId()==23){
							idx_daughterA = 0;
							idx_daughterB = 2;
						}
						else{
							idx_daughterA = 0;
							idx_daughterB = 1;
						}

						const reco::Candidate* daughterA = mcParticle.daughter(idx_daughterA);
						const reco::Candidate* daughterB = mcParticle.daughter(idx_daughterB);

						TLorentzVector daughterA_p4 = ConvertToTLorentzVector( &(daughterA->p4()) );
						TLorentzVector daughterB_p4 = ConvertToTLorentzVector( &(daughterB->p4()) );

						mcZboson_daughters_dR_   = daughterA_p4.DeltaR(daughterB_p4);
						mcZboson_daughters_dEta_ = daughterB_p4.Eta() - daughterA_p4.Eta();
						mcZboson_daughters_dPhi_ = daughterB_p4.Phi() - daughterA_p4.Phi();
						while( mcZboson_daughters_dPhi_<=-1.0*TMath::Pi() ){ mcZboson_daughters_dPhi_ = mcZboson_daughters_dPhi_+2.0*TMath::Pi(); }
						while( mcZboson_daughters_dPhi_>=TMath::Pi()      ){ mcZboson_daughters_dPhi_ = mcZboson_daughters_dPhi_-2.0*TMath::Pi(); }
						mcZboson_daughters_openingAngle_ = daughterA_p4.Angle(daughterB_p4.Vect());
						mcZboson_daughterA_p4_ = daughterA->p4();
						mcZboson_daughterB_p4_ = daughterB->p4();
					}
					if(vBool_){std::cout << "           " << "Daughters' quantities: dR=" << mcZboson_daughters_dR_ << "; dEta=" << mcZboson_daughters_dEta_ << "; dPhi=" << mcZboson_daughters_dPhi_ << "; ang sep'n=" << mcZboson_daughters_openingAngle_ << std::endl;}
				}
				//For Z bosons with other status numbers
				else
					mc_numZs_statusNot2Or3_++;
			}

			//Finding the electrons and positrons...
			if(abs(particlePDGid)==11){
				if(particlePDGid==11 && vBool_)
					std::cout << "   Generator-particle no. " << idx << " is an electron. (pdgId= " << particlePDGid << "; status=" << particleStatus << "; Pt=" << particle4mom.Pt() << ")" << std::endl;
				else if(vBool_)
					std::cout << "   Generator-particle no. " << idx << " is a positron.  (pdgId= " << particlePDGid << "; status=" << particleStatus << "; Pt=" << particle4mom.Pt() << ")" << std::endl;
				if(vBool_){
					std::cout << "      Has " << mcParticle.numberOfMothers() << " mothers;" << std::endl;
					for(unsigned int idxMother=0; idxMother<mcParticle.numberOfMothers(); idxMother++)
						std::cout << "         Mother(" << idxMother << ") has PDG ID " << mcParticle.mother(idxMother)->pdgId() << " and status=" << mcParticle.mother(idxMother)->status() << std::endl;
				}
				if(particleStatus==1){ //If the electron is a final state ele...
				numFinalStateEles++;
					if(particle4mom.Pt()>eleHighestEt){ //If electron has greater Pt than any electron before it...
						idx_ele2ndHighestEt = idx_eleHighestEt; 	// Switch the old highest Pt electron to be the 2nd highest Pt ele...
						ele2ndHighestEt = eleHighestEt;
						idx_eleHighestEt = idx; 						// .. and set this ele as the highest Pt ele.
						eleHighestEt = particle4mom.Pt();
						if(vBool_){std::cout << "      (Now set as highest Pt electron...)" << std::endl;
							std::cout << "      (... and #" << idx_ele2ndHighestEt << " set as 2nd highest Pt ele...)" << std::endl;}
					}
					else if(particle4mom.Pt()>ele2ndHighestEt){ 	// Otherwise, if ele Pt is still larger than previous 2nd highest ele Pt...
						idx_ele2ndHighestEt = idx;						// ...set this ele as the 2nd highest PT ele.
						ele2ndHighestEt = particle4mom.Pt();
						if(vBool_){std::cout << "      (Now set as 2nd highest Pt electron...)" << std::endl;}
					}
				}
			}

			//Finding the daughters of the Z boson, and then reconstructing the Z boson, for the signal 2010 datasets in which the intermediate Feynaman diagram particles - i.e.
			// the u* and the Z boson - are not stored in the signal datasets that were generated at the start of 2011
			if(is2010SignalDataset_){
				// If statement to select those mc particles with two mothers - one of which is a gluon (pdgId=21), and the other is a up quark (PDG id 2)
				if(mcParticle.numberOfMothers()==2){
					if( (abs(mcParticle.mother(0)->pdgId())==2 && abs(mcParticle.mother(1)->pdgId())==21) || (abs(mcParticle.mother(0)->pdgId())==21 && abs(mcParticle.mother(1)->pdgId())==2)  ){
						// If this particle is an electron (pdgId==11), then assign it to be daughter A ...
						if(mcParticle.pdgId()==11)
							mcZboson_daughterA_p4_ = mcParticle.p4();
						// Otherwise, if this particle is an positron (pdgId==-11), then assign it to be daughter B ...
						else if(mcParticle.pdgId()==-11)
							mcZboson_daughterB_p4_ = mcParticle.p4();
					}
				}

			}

		}// End of for loop over MC particles

		/*std::cout << " ->There are " << numZs_status2 << " Z bosons with status=2 in this event." << std::endl;
		numZs_status2Hist->Fill(numZs_status2);
		std::cout << "   There are " << numZs_status3 << " Z bosons with status=3 in this event." << std::endl;
		numZs_status3Hist->Fill(numZs_status3);
		std::cout << "   There are " << numZs_statusOther << " Z bosons with other status numbers in this event." << std::endl;
		numZs_statusOtherHist->Fill(numZs_statusOther);*/

		//Storing branches of information about highest Et ele...
		if(numFinalStateEles>0){
			if(vBool_){std::cout << " ->The highest Et final-state electron in the event was genParticle#" << idx_eleHighestEt << ", with Et=" << eleHighestEt << std::endl;}
			eleHighestEt_EtHist->Fill(eleHighestEt);
			const reco::GenParticle& mcEle_HighestEt = (*genParticles)[idx_eleHighestEt];
			mcEles_HighestEt_charge = mcEle_HighestEt.charge(); //Int_t
			mcEles_HighestEt_PDGid  = mcEle_HighestEt.pdgId(); //Int_t
			mcEles_HighestEt_status = mcEle_HighestEt.status(); //Int_t
			mcEles_HighestEt_pt     = mcEle_HighestEt.p4().Pt(); //Double_t
			mcEles_HighestEt_p4.SetPxPyPzE(mcEle_HighestEt.p4().Px(),mcEle_HighestEt.p4().Py(),mcEle_HighestEt.p4().Pz(),mcEle_HighestEt.p4().E());

			if(vBool_){std::cout << "       charge=" << mcEles_HighestEt_charge << "; pdgId=" << mcEles_HighestEt_PDGid << "; status=" << mcEles_HighestEt_status << std::endl;
				std::cout << "       P=" << mcEles_HighestEt_p4.P() << "; Pt=" << mcEles_HighestEt_pt << std::endl;
				std::cout << "       eta=" << mcEles_HighestEt_p4.Eta() << "; phi=" << mcEles_HighestEt_p4.Phi() << std::endl;}
		}
		else{
			mcEles_HighestEt_charge = -999;
			mcEles_HighestEt_PDGid  = -999;
			mcEles_HighestEt_status = -999;
			mcEles_HighestEt_pt = -999.9;
			mcEles_HighestEt_p4.SetPxPyPzE(0.0,0.0,0.0,-999.9);
		}

		//Storing branches of information about 2nd highest Et ele...
		if(numFinalStateEles>1){
			if(vBool_){std::cout << "   The 2nd highest Et final-state electron in the event was genParticle#" << idx_ele2ndHighestEt << ", with Et=" << ele2ndHighestEt << std::endl;}
			ele2ndHighestEt_EtHist->Fill(ele2ndHighestEt);
			const reco::GenParticle& mcEle_2ndHighestEt = (*genParticles)[idx_ele2ndHighestEt];
			mcEles_2ndHighestEt_charge = mcEle_2ndHighestEt.charge(); //Int_t
			mcEles_2ndHighestEt_PDGid  = mcEle_2ndHighestEt.pdgId(); //Int_t
			mcEles_2ndHighestEt_status = mcEle_2ndHighestEt.status(); //Int_t
			mcEles_2ndHighestEt_pt     = mcEle_2ndHighestEt.p4().Pt(); //Double_t
			mcEles_2ndHighestEt_p4.SetPxPyPzE(mcEle_2ndHighestEt.p4().Px(),mcEle_2ndHighestEt.p4().Py(),mcEle_2ndHighestEt.p4().Pz(),mcEle_2ndHighestEt.p4().E());

			if(vBool_){std::cout << "       charge=" << mcEles_2ndHighestEt_charge << "; pdgId=" << mcEles_2ndHighestEt_PDGid << "; status=" << mcEles_2ndHighestEt_status << std::endl;
				std::cout << "       P=" << mcEles_2ndHighestEt_p4.P() << "; Pt=" << mcEles_2ndHighestEt_pt << std::endl;
				std::cout << "       eta=" << mcEles_2ndHighestEt_p4.Eta() << "; phi=" << mcEles_2ndHighestEt_p4.Phi() << std::endl;}
		}
		else{
			mcEles_2ndHighestEt_charge = -999;
			mcEles_2ndHighestEt_PDGid  = -999;
			mcEles_2ndHighestEt_status = -999;
			mcEles_2ndHighestEt_pt = -999.9;
			mcEles_2ndHighestEt_p4.SetPxPyPzE(0.0, 0.0, 0.0, -999.9);
		}
		//Storing branches of information about the Z candidate ...
		SetMCZcandidateVariables(vBool_);

		// Sorting out the values of the mcZboson variables in the case that am running over the signal datasets that were generated at the start 2011
		if(is2010SignalDataset_){
			// Two new lines put in in 42X on 2011/08/31, since previous code that was used to find Z boson daughters by looking at mother PDG ID's doesn't seem to work anymore ...
			mcZboson_daughterA_p4_ = mcEles_HighestEt_p4;
			mcZboson_daughterB_p4_ = mcEles_2ndHighestEt_p4;

			mcZboson_pdgId_  = 0; // Int_t
			mcZboson_status_ = -100; // Int_t
			mcZboson_p4_     = mcZboson_daughterA_p4_ + mcZboson_daughterB_p4_;
			mcZboson_numDaughters_ = 2;

			TLorentzVector daughterA_p4 = ConvertToTLorentzVector( &mcZboson_daughterA_p4_ );
			TLorentzVector daughterB_p4 = ConvertToTLorentzVector( &mcZboson_daughterB_p4_ );

			mcZboson_daughters_dR_   = daughterA_p4.DeltaR(daughterB_p4);
			mcZboson_daughters_dEta_ = daughterB_p4.Eta() - daughterA_p4.Eta();
			mcZboson_daughters_dPhi_ = daughterB_p4.Phi() - daughterA_p4.Phi();
			while( mcZboson_daughters_dPhi_<=-1.0*TMath::Pi() ){ mcZboson_daughters_dPhi_ = mcZboson_daughters_dPhi_+2.0*TMath::Pi(); }
			while( mcZboson_daughters_dPhi_>=TMath::Pi()      ){ mcZboson_daughters_dPhi_ = mcZboson_daughters_dPhi_-2.0*TMath::Pi(); }
			mcZboson_daughters_openingAngle_ = daughterA_p4.Angle(daughterB_p4.Vect());

			//Printing to screen for debugging ...
			if(vBool_){
				//std::cout << "   The MC Z boson variables have been determined from the electron and positron which each have two mothers - a u* and a u/ubar ... " << std::endl;
				std::cout << "   The MC Z boson variables have been determined from the two highest Et e+- ... " << std::endl;
				std::cout << "      => pdgId=" << mcZboson_pdgId_ << "; status=" << mcZboson_status_ << "; (pT,eta,phi) = (" << mcZboson_p4_.Pt() << ", " << mcZboson_p4_.Eta() << ", " << mcZboson_p4_.Phi() << ")" << std::endl;
				std::cout << "      It has " << mcZboson_numDaughters_ << " daughters" << std::endl;
				std::cout << "           Daughter A: (pT,eta,phi) = (" << mcZboson_daughterA_p4_.Pt() << "," << mcZboson_daughterA_p4_.Eta() << "," << mcZboson_daughterA_p4_.Phi() << ")" << std::endl;
				std::cout << "           Daughter B: (pT,eta,phi) = (" << mcZboson_daughterB_p4_.Pt() << "," << mcZboson_daughterB_p4_.Eta() << "," << mcZboson_daughterB_p4_.Phi() << ")" << std::endl;
				std::cout << "      Relative to each other:" << std::endl;
				std::cout << "           dR=" << mcZboson_daughters_dR_ << "; dEta=" << mcZboson_daughters_dEta_ << "; dPhi=" << mcZboson_daughters_dPhi_ << "; opening angle=" << mcZboson_daughters_openingAngle_ << std::endl;
			}
		}

		if(vBool_){std::cout << std::endl;}
	}
	
	//Reset branch addresses here each time...
	if(zBosonDaughter_PDGid==0)
		zBosonDaughter_PDGid = zBosonDaughterIfSherpa_PDGid;
	if(dyJetsToLL_EventType_==0){
		numEvtsStored_++;
		EventDataTree->Fill();
	}
	else if( (dyJetsToLL_EventType_!=0) && (abs(zBosonDaughter_PDGid)==dyJetsToLL_EventType_) ){
		numEvtsStored_++;
		EventDataTree->Fill();
	}

	delete event_;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
BstdZeeNTupler::beginJob()
{
}

// ------------ method called when starting to processes a run  ------------
void 
BstdZeeNTupler::beginRun(const edm::Run& run, const edm::EventSetup& iSetup){

	if(readInTrigInfo_){
		// HLT setup
		bool changed; //This variable is not used by the init method of hltConfig_ and so does not need to be initialised with a value here
		hltConfig_.init(run, iSetup, hltResultsTag_.process(), changed); //The value of changed now indicates whether the HLT configuration has changed with respect to the previous run or not.
		//printDatasetsAndTriggerNames();

		// Retrieve the di-ele signal trigger name and index ...
		std::pair<std::string, unsigned int> diEleSigTrig_nameIdxPair = beginRun_getTriggerNameAndIdx(hltPathA_possNames_, "di-ele trigger");
		hltPathA_     = diEleSigTrig_nameIdxPair.first;
		hltPathIndex_ = diEleSigTrig_nameIdxPair.second;

		// Retrieve the e-mu trigger name and index ...
		std::pair<std::string, unsigned int> eMuTrigger_nameIdxPair = beginRun_getTriggerNameAndIdx(trg_emuPath_possNames_, "e-mu trigger");
		trg_emuPath_name_ = eMuTrigger_nameIdxPair.first;
		trg_emuPath_idx_ = eMuTrigger_nameIdxPair.second;
	}

}


std::pair<std::string, unsigned int> BstdZeeNTupler::beginRun_getTriggerNameAndIdx(const std::vector<std::string>& vecOfRegexs, std::string triggerId_debug)
{ // N.B: Must be called AFTER hltConfig_.init() ....
	unsigned int pathIdx = 999999;
	std::string pathName = " *** DEFAULT TRIGGER NAME *** ";

	// Run through all of the possible trigger names for this trigger ...
	std::cout << std::endl << " *****************************************************" << std::endl;
	std::cout << " * Looking through the supplied trigger regexs ...   [Trig table name: '" << hltConfig_.tableName() << "' ]" << std::endl;
	for(std::vector<std::string>::const_iterator regexStringIt = vecOfRegexs.begin(); regexStringIt != vecOfRegexs.end(); regexStringIt++){
		std::cout << "     * " << *regexStringIt << std::endl;;
		std::vector<std::string> matchingTriggerNames = HLTConfigProvider::matched( hltConfig_.triggerNames() , *regexStringIt);

		if(matchingTriggerNames.size()>0){
			pathName = matchingTriggerNames.at(0);
			pathIdx = hltConfig_.triggerIndex( pathName );
			if(pathIdx < hltConfig_.size()){
				std::cout << "       ->Have found matching trigger name (\"" << pathName << "\", idx=" << pathIdx << ") in menu! Exiting for loop now ..." << std::endl;
				std::cout << "         ( There are " << hltConfig_.size() << " triggers in the menu. )" << std::endl << std::endl;
				return std::pair<std::string, unsigned int>(pathName,pathIdx);
			}
			else
				std::cout << "       ->This name isn't in menu! Trying next one now ..." << std::endl;
		} // End-if : match found in menu
	}

	// If get here then none trigger names provided could be found in the menu then throw an exception ...
	cms::Exception ex("readingTrigerInfo");
	ex << " [In BstdZeeNTupler::beginRun_getTriggerNameAndIdx ...]" << std::endl;
	ex << "      None of the supplied possible names for trigger '" << triggerId_debug << "' were found in the trigger menu "<< std::endl;
	ex << "      (See FILE: '" << __FILE__ << "'; LINE NO. " << __LINE__ << ")" << std::endl;
	throw ex;
}

// ------------ method called once each job just after ending the event loop  ------------
void BstdZeeNTupler::endJob()
{
	int totalNumHighestEtEles = numEvts;
	double dbl_totalNumHighestEtEles = static_cast<double>(totalNumHighestEtEles); 
	double dbl_numHighestEtElesAboveThreshold = -999.9;
	double effiOfEtThreshold = -999.9;
	double error_effiOfEtThreshold = -999.9;
	
	//Setting errors on histogram bins to be sum of weights for bin entries (i.e. Poisson errors)
	//status2Zs_PtHist->Sumw2();
	//status2Zs_LogPtHist->Sumw2();
	ele2ndHighestEt_EtHist->Sumw2();
	
	//Calculating the fraction of events passing a single electron Et threshold for the numAboveSingleEleEtThreshold histogram...
	std::cout << std::endl << "Program is now calculating the fraction of events passing a single electron Et threshold..." << std::endl;
	for(int binIdx=0; binIdx<=histoEtNBins+1; binIdx++){
		//std::cout << "  *binIdx=" << binIdx << std::endl;
		dbl_numHighestEtElesAboveThreshold = eleHighestEt_EtHist->Integral(binIdx,histoEtNBins+1);
		effiOfEtThreshold = dbl_numHighestEtElesAboveThreshold/dbl_totalNumHighestEtEles;
		error_effiOfEtThreshold = sqrt(effiOfEtThreshold*(1.0-effiOfEtThreshold)/dbl_totalNumHighestEtEles);
		//std::cout << "     dbl_numHighestEtElesAboveThreshold = " << dbl_numHighestEtElesAboveThreshold << std::endl;
		//std::cout << "     dbl_totalNumHighestEtEles = " << dbl_totalNumHighestEtEles << std::endl;
		//std::cout << "     effiOfEtThreshold = " << effiOfEtThreshold << std::endl;
		
		fracAboveSingleEleEtThreshold_Hist->SetBinContent(binIdx,effiOfEtThreshold);
		fracAboveSingleEleEtThreshold_Hist->SetBinError(binIdx,error_effiOfEtThreshold);
	}
	std::cout << "   done." << std::endl;
	
	//std::cout << std::endl <<  "**histoEtBinWidth=" << (histoEtMax-histoEtMin)/static_cast<double>(histoEtNBins) << std::endl;

	// Filling the cumulative histogram: highestEtTrigObjPathA_EtHist_
	Int_t trgEtHist_nBins = highestEtTrigObjPathA_EtHist_->GetNbinsX();
	Double_t trgEtHist_tmpBinContent = 0.0;
	Double_t trgEtHist_tmpEffi = 0.0;
	for(Int_t binIdx=trgEtHist_nBins+1; binIdx>=0; binIdx--){
		if(binIdx==trgEtHist_nBins+1)
			trgEtHist_tmpBinContent = highestEtTrigObjPathA_EtHist_->GetBinContent(binIdx);
		else
			trgEtHist_tmpBinContent = highestEtTrigObjPathA_cumEtHist_->GetBinContent(binIdx+1) + highestEtTrigObjPathA_EtHist_->GetBinContent(binIdx);

		highestEtTrigObjPathA_cumEtHist_->SetBinContent(binIdx, trgEtHist_tmpBinContent);
		highestEtTrigObjPathA_cumEtHist_->SetBinError(binIdx, sqrt(trgEtHist_tmpBinContent));
		highestEtTrigObjPathA_EtThrDenomHist_->SetBinContent(binIdx, static_cast<double>(numEvts));

		trgEtHist_tmpEffi = highestEtTrigObjPathA_cumEtHist_->GetBinContent(binIdx)/highestEtTrigObjPathA_EtThrDenomHist_->GetBinContent(binIdx);
		highestEtTrigObjPathA_EtThrEffiHist_->SetBinContent(binIdx, trgEtHist_tmpEffi);
		highestEtTrigObjPathA_EtThrEffiHist_->SetBinError(binIdx, sqrt( trgEtHist_tmpEffi*(1.0-trgEtHist_tmpEffi)/highestEtTrigObjPathA_EtThrDenomHist_->GetBinContent(binIdx) ) );
	}
	
	highestEtTrigObjPathA_EtThrAsymmEffi_ = fHistos->make<TGraphAsymmErrors>(highestEtTrigObjPathA_cumEtHist_, highestEtTrigObjPathA_EtThrDenomHist_);
	highestEtTrigObjPathA_EtThrAsymmEffi_->SetName("highestEtTrigObjPathA_EtThrAsymmEffi");
	highestEtTrigObjPathA_EtThrAsymmEffi_->SetTitle("Trigger efficiency for signal events vs. signal trigger E_{T} threshold");
	highestEtTrigObjPathA_EtThrAsymmEffi_->GetXaxis()->SetTitle("HLT electron E_{T} threshold, E_{T}^{thr} /GeV");
	highestEtTrigObjPathA_EtThrAsymmEffi_->GetYaxis()->SetTitle("Efficiency for signal events");

	fHistos->cd();
	EventDataTree->Write();

	if(!calcPileUpWeights_)
		std::cout << std::endl << "   *** WARNING : PU weights have *NOT* been calculated (relevant PU dist'n file & hist names not specified in config) !" << std::endl;

	std::cout << std::endl;
	std::cout << " ***************************************" << std::endl;
	std::cout << " * " << numEvtsStored_ << " events have been stored in the NTuple";
	std::cout << std::endl;
}

/*// ------------ method called when ending the processing of a run  ------------
void 
BstdZeeNTupler::endRun(edm::Run const&, edm::EventSetup const&)
{
}*/

/*// ------------ method called when starting to processes a luminosity block  ------------
void 
BstdZeeNTupler::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
BstdZeeNTupler::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}*/

// ------------ methods for setting up the pointer links between the branches and the variables...
void BstdZeeNTupler::SetupMCBranches(){
	
	EventDataTree->Branch("mcEles_number", &numFinalStateEles, "mcEles_number/i"); //UInt_t

	EventDataTree->Branch("mcEles_HighestEt_charge", &mcEles_HighestEt_charge, "mcEles_HighestEt_charge/I"); //Int_t
	EventDataTree->Branch("mcEles_HighestEt_PDGid",  &mcEles_HighestEt_PDGid, "mcEles_HighestEt_charge/I"); //Int_t
	EventDataTree->Branch("mcEles_HighestEt_status", &mcEles_HighestEt_status, "mcEles_HighestEt_status/I"); //Int_t
	EventDataTree->Branch("mcEles_HighestEt_pt",     &mcEles_HighestEt_pt, "mcEles_HighestEt_pt/D"); //Double_t
	EventDataTree->Branch("mcEles_HighestEt_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcEles_HighestEt_p4_ptr);
	mcEles_HighestEt_p4_ptr = &mcEles_HighestEt_p4;

	EventDataTree->Branch("mcEles_2ndHighestEt_charge", &mcEles_2ndHighestEt_charge, "mcEles_2ndHighestEt_charge/I"); //Int_t
	EventDataTree->Branch("mcEles_2ndHighestEt_PDGid",  &mcEles_2ndHighestEt_PDGid,  "mcEles_2ndHighestEt_charge/I"); //Int_t
	EventDataTree->Branch("mcEles_2ndHighestEt_status", &mcEles_2ndHighestEt_status, "mcEles_2ndHighestEt_status/I"); //Int_t
	EventDataTree->Branch("mcEles_2ndHighestEt_pt",     &mcEles_2ndHighestEt_pt,     "mcEles_2ndHighestEt_pt/D"); //Double_t
	EventDataTree->Branch("mcEles_2ndHighestEt_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcEles_2ndHighestEt_p4_ptr);
	mcEles_2ndHighestEt_p4_ptr = &mcEles_2ndHighestEt_p4;

	/*EventDataTree->Branch("mcEles_charge", &mcEles_charge);
	EventDataTree->Branch("mcEles_PDGid",  &mcEles_PDGid);
	EventDataTree->Branch("mcEles_status", &mcEles_status);
	EventDataTree->Branch("mcEles_pt",     &mcEles_pt);
	EventDataTree->Branch("mcEles_p4", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &mcEles_p4);*/

	EventDataTree->Branch("mcZcandidate_pt",   &mcZcandidate_pt,  "mcZcandidate_pt/D"  ); //Double_t
	EventDataTree->Branch("mcZcandidate_eta",  &mcZcandidate_eta, "mcZcandidate_eta/D" ); //Double_t
	EventDataTree->Branch("mcZcandidate_phi",  &mcZcandidate_phi, "mcZcandidate_phi/D" ); //Double_t
	EventDataTree->Branch("mcZcandidate_mass", &mcZcandidate_mass,"mcZcandidate_mass/D"); //Double_t
	EventDataTree->Branch("mcZcandidate_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcZcandidate_p4_ptr); //ROOT::Math::XYZTLorentzVector
	mcZcandidate_p4_ptr = &mcZcandidate_p4;
	EventDataTree->Branch("mcZcandidate_dEtaEles",     &mcZcandidate_dEtaEles,    "mcZcandidate_dEtaEles/D");     //Double_t
	EventDataTree->Branch("mcZcandidate_dPhiEles",     &mcZcandidate_dPhiEles,    "mcZcandidate_dPhiEles/D");     //Double_t
	EventDataTree->Branch("mcZcandidate_dREles",       &mcZcandidate_dREles,      "mcZcandidate_dREles/D");       //Double_t
	EventDataTree->Branch("mcZcandidate_openingAngle", &mcZcandidate_openingAngle,"mcZcandidate_openingAngle/D"); //Double_t

	EventDataTree->Branch("mc_numZs_status2",       &mc_numZs_status2_,       "mc_numZs_status2/I"); // Int_t
	EventDataTree->Branch("mc_numZs_status3",       &mc_numZs_status3_,       "mc_numZs_status3/I"); // Int_t
	EventDataTree->Branch("mc_numZs_statusNot2Or3", &mc_numZs_statusNot2Or3_, "mc_numZs_statusNot2Or3/I"); // Int_t

	EventDataTree->Branch("mcZboson_pdgId",  &mcZboson_pdgId_,  "mcZboson_pdgId/I"); // Int_t
	EventDataTree->Branch("mcZboson_status", &mcZboson_status_, "mcZboson_status/I"); // Int_t
	EventDataTree->Branch("mcZboson_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcZboson_p4_ptr_); // ROOT::Math::XYZTLorentzVector
	mcZboson_p4_ptr_ = &mcZboson_p4_;
	EventDataTree->Branch("mcZboson_status2_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcZboson_status2_p4_ptr_); // ROOT::Math::XYZTLorentzVector
	mcZboson_status2_p4_ptr_ = &mcZboson_status2_p4_;
	EventDataTree->Branch("mcZboson_numDaughters",           &mcZboson_numDaughters_,           "mcZboson_numDaughters/I");
	EventDataTree->Branch("mcZboson_daughters_dR",           &mcZboson_daughters_dR_,           "mcZboson_daughters_dR/D");
	EventDataTree->Branch("mcZboson_daughters_dEta",         &mcZboson_daughters_dEta_,         "mcZboson_daughters_dEta/D");
	EventDataTree->Branch("mcZboson_daughters_dPhi",         &mcZboson_daughters_dPhi_,         "mcZboson_daughters_dPhi/D");
	EventDataTree->Branch("mcZboson_daughters_openingAngle", &mcZboson_daughters_openingAngle_, "mcZboson_daughters_openingAngle/D");
	EventDataTree->Branch("mcZboson_daughterA_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcZboson_daughterA_p4_ptr_); // ROOT::Math::XYZTLorentzVector
	mcZboson_daughterA_p4_ptr_ = &mcZboson_daughterA_p4_;
	EventDataTree->Branch("mcZboson_daughterB_p4", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", &mcZboson_daughterB_p4_ptr_); // ROOT::Math::XYZTLorentzVector
	mcZboson_daughterB_p4_ptr_ = &mcZboson_daughterB_p4_;
}

void BstdZeeNTupler::SetupTriggerBranches(){
	//EventDataTree->Branch("trg_TrigPathDecisions", &trg_TrigPathDecisions_); // std::vector<Bool_t>
	//EventDataTree->Branch("trg_TrigPathNames", &trg_TrigPathNames);         // std::vector<std::string>
	EventDataTree->Branch("trg_PathA_decision", &trg_PathA_decision_);       // Bool_t
	EventDataTree->Branch("trg_PathA_name", &trg_PathA_name_);               // std::string hltPathA_
	EventDataTree->Branch("trg_PathA_nameOfLastFilter", &hltPathA_nameOfLastFilter_); // std::string
	EventDataTree->Branch("trg_PathA_highestTrigObjEt", &hltPathA_highestTrigObjEt_,"trg_PathA_highestTrigObjEt/F"); // float
}

void BstdZeeNTupler::SetupEvtInfoBranches(){
	//Setting up the pointer links between the event information variables and the corresponding branches...
	EventDataTree->Branch("evt_runNum", &evt_runNum_, "evt_runNum/i"); //unsigned int 
	EventDataTree->Branch("evt_lumiSec", &evt_lumiSec_, "evt_lumiSec/i"); //unsigned int 
	EventDataTree->Branch("evt_evtNum", &evt_evtNum_, "evt_evtNum/i"); //unsigned int 
}

void BstdZeeNTupler::SetupStdEleBranches(){
	//Setting up the pointer links between the standard GSF electron variables and the corresponding branches...
	EventDataTree->Branch("normGsfEles_number", &normGsfEles_number_, "normGsfEles_number/i"); //unsigned int
	EventDataTree->Branch("normGsfEles_p4ptr_", &normGsfEles_p4ptr_); //std::vector<ROOT::Math::XYZTVector>*
	normGsfEles_p4ptr_ = &normGsfEles_p4_;	
	EventDataTree->Branch("normGsfEles_charge", &normGsfEles_charge_); //std::vector<Int_t>
	
	EventDataTree->Branch("normGsfEles_Et", &normGsfEles_Et_);
	EventDataTree->Branch("normGsfEles_HEEP_Et", &normGsfEles_HEEP_Et_);
	EventDataTree->Branch("normGsfEles_Eta", &normGsfEles_Eta_);
	EventDataTree->Branch("normGsfEles_scEta", &normGsfEles_scEta_);
	EventDataTree->Branch("normGsfEles_ecalDriven", &normGsfEles_ecalDriven_);
	EventDataTree->Branch("normGsfEles_ecalDrivenSeed", &normGsfEles_ecalDrivenSeed_);
	
	EventDataTree->Branch("normGsfEles_HEEP_dEtaIn", &normGsfEles_HEEP_dEtaIn_);
	EventDataTree->Branch("normGsfEles_HEEP_dPhiIn", &normGsfEles_HEEP_dPhiIn_);
	EventDataTree->Branch("normGsfEles_HoverE", &normGsfEles_HoverE_);
	EventDataTree->Branch("normGsfEles_sigmaIetaIeta", &normGsfEles_sigmaIetaIeta_);
	EventDataTree->Branch("normGsfEles_scSigmaIetaIeta", &normGsfEles_scSigmaIetaIeta_);
	
	EventDataTree->Branch("normGsfEles_dr03EmIsoEt", &normGsfEles_dr03EmIsoEt_);
	EventDataTree->Branch("normGsfEles_dr03HadDepth1IsoEt", &normGsfEles_dr03HadDepth1IsoEt_);
	EventDataTree->Branch("normGsfEles_dr03HadDepth2IsoEt", &normGsfEles_dr03HadDepth2IsoEt_);
	EventDataTree->Branch("normGsfEles_dr03TkIsoPt", &normGsfEles_dr03TkIsoPt_);
	
	EventDataTree->Branch("normGsfEles_e2x5Max", &normGsfEles_e2x5Max_);
	EventDataTree->Branch("normGsfEles_e5x5", &normGsfEles_e5x5_);

	//normGsfEles_tswHEEPElePtr_ = &normGsfEles_tswHEEPEle_;
	//EventDataTree->Branch("normGsfEles_tswHEEPElePtr_", &normGsfEles_tswHEEPElePtr_);

	//kinematic and geometric methods
	EventDataTree->Branch("normHEEPEles_et",         &normHEEPEles_et_);
	EventDataTree->Branch("normHEEPEles_gsfEt",      &normHEEPEles_gsfEt_);
	EventDataTree->Branch("normHEEPEles_scEt",       &normHEEPEles_scEt_);
	EventDataTree->Branch("normHEEPEles_energy",     &normHEEPEles_energy_);
	EventDataTree->Branch("normHEEPEles_gsfEnergy",  &normHEEPEles_gsfEnergy_);
	EventDataTree->Branch("normHEEPEles_caloEnergy", &normHEEPEles_caloEnergy_);
	EventDataTree->Branch("normHEEPEles_ecalEnergyError", &normHEEPEles_ecalEnergyError_);
	EventDataTree->Branch("normHEEPEles_eta",        &normHEEPEles_eta_);
	EventDataTree->Branch("normHEEPEles_scEta",      &normHEEPEles_scEta_);
	EventDataTree->Branch("normHEEPEles_detEta",     &normHEEPEles_detEta_);
	EventDataTree->Branch("normHEEPEles_detEtaAbs",  &normHEEPEles_detEtaAbs_);
	EventDataTree->Branch("normHEEPEles_phi",        &normHEEPEles_phi_);
	EventDataTree->Branch("normHEEPEles_scPhi",      &normHEEPEles_scPhi_);
	EventDataTree->Branch("normHEEPEles_detPhi",     &normHEEPEles_detPhi_);
	EventDataTree->Branch("normHEEPEles_zVtx",       &normHEEPEles_zVtx_);
	EventDataTree->Branch("normHEEPEles_p4ptr",      &normHEEPEles_p4ptr_);
	normHEEPEles_p4ptr_ = &normHEEPEles_p4_;
	EventDataTree->Branch("normHEEPEles_gsfP4ptr",   &normHEEPEles_gsfP4ptr_);
	normHEEPEles_gsfP4ptr_ = &normHEEPEles_gsfP4_;

	//classification (couldnt they have just named it 'type')
	EventDataTree->Branch("normHEEPEles_classification",  &normHEEPEles_classification_);
	EventDataTree->Branch("normHEEPEles_isEcalDriven",    &normHEEPEles_isEcalDriven_);
	EventDataTree->Branch("normHEEPEles_isTrackerDriven", &normHEEPEles_isTrackerDriven_);
	EventDataTree->Branch("normHEEPEles_isEB",            &normHEEPEles_isEB_);
	EventDataTree->Branch("normHEEPEles_isEE",            &normHEEPEles_isEE_);

	//track methods
	EventDataTree->Branch("normHEEPEles_charge",    &normHEEPEles_charge_);
	EventDataTree->Branch("normHEEPEles_trkCharge", &normHEEPEles_trkCharge_);
	EventDataTree->Branch("normHEEPEles_pVtx",      &normHEEPEles_pVtx_);
	EventDataTree->Branch("normHEEPEles_pCalo",     &normHEEPEles_pCalo_);
	EventDataTree->Branch("normHEEPEles_ptVtx",     &normHEEPEles_ptVtx_);
	EventDataTree->Branch("normHEEPEles_ptCalo",    &normHEEPEles_ptCalo_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_pt", &normHEEPEles_closestCtfTrk_pt_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_eta", &normHEEPEles_closestCtfTrk_eta_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_phi", &normHEEPEles_closestCtfTrk_phi_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_innerPt", &normHEEPEles_closestCtfTrk_innerPt_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_innerEta", &normHEEPEles_closestCtfTrk_innerEta_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_innerPhi", &normHEEPEles_closestCtfTrk_innerPhi_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_outerPt", &normHEEPEles_closestCtfTrk_outerPt_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_outerEta", &normHEEPEles_closestCtfTrk_outerEta_);
	EventDataTree->Branch("normHEEPEles_closestCtfTrk_outerPhi", &normHEEPEles_closestCtfTrk_outerPhi_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	EventDataTree->Branch("normHEEPEles_hOverE",       &normHEEPEles_hOverE_);
	EventDataTree->Branch("normHEEPEles_dEtaIn",       &normHEEPEles_dEtaIn_);
	EventDataTree->Branch("normHEEPEles_dPhiIn",       &normHEEPEles_dPhiIn_);
	EventDataTree->Branch("normHEEPEles_dPhiOut",      &normHEEPEles_dPhiOut_);
	EventDataTree->Branch("normHEEPEles_epIn",         &normHEEPEles_epIn_);
	EventDataTree->Branch("normHEEPEles_epOut",        &normHEEPEles_epOut_);
	EventDataTree->Branch("normHEEPEles_fbrem",        &normHEEPEles_fbrem_);
	EventDataTree->Branch("normHEEPEles_bremFrac",     &normHEEPEles_bremFrac_);
	EventDataTree->Branch("normHEEPEles_invEOverInvP", &normHEEPEles_invEOverInvP_);

	//shower shape variables
	EventDataTree->Branch("normHEEPEles_sigmaEtaEta",       &normHEEPEles_sigmaEtaEta_);
	EventDataTree->Branch("normHEEPEles_sigmaEtaEtaUnCorr", &normHEEPEles_sigmaEtaEtaUnCorr_);
	EventDataTree->Branch("normHEEPEles_sigmaIEtaIEta",     &normHEEPEles_sigmaIEtaIEta_);
	EventDataTree->Branch("normHEEPEles_e1x5",              &normHEEPEles_e1x5_);
	EventDataTree->Branch("normHEEPEles_e2x5Max",           &normHEEPEles_e2x5Max_);
	EventDataTree->Branch("normHEEPEles_e5x5",              &normHEEPEles_e5x5_);
	EventDataTree->Branch("normHEEPEles_e1x5Over5x5",       &normHEEPEles_e1x5Over5x5_);
	EventDataTree->Branch("normHEEPEles_e2x5MaxOver5x5",    &normHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	EventDataTree->Branch("normHEEPEles_isolEm",            &normHEEPEles_isolEm_);
	EventDataTree->Branch("normHEEPEles_isolHad",           &normHEEPEles_isolHad_);
	EventDataTree->Branch("normHEEPEles_isolHadDepth1",     &normHEEPEles_isolHadDepth1_);
	EventDataTree->Branch("normHEEPEles_isolHadDepth2",     &normHEEPEles_isolHadDepth2_);
	EventDataTree->Branch("normHEEPEles_isolPtTrks",        &normHEEPEles_isolPtTrks_);
	EventDataTree->Branch("normHEEPEles_isolEmHadDepth1",   &normHEEPEles_isolEmHadDepth1_);

	// Number of missing inner hits
	EventDataTree->Branch("normHEEPEles_numMissInnerHits", &normHEEPEles_numMissInnerHits_);

	// SC recHits information - for correcting calo-based iso variables
	EventDataTree->Branch("normHEEPEles_SCposn_eta", &normHEEPEles_SCposn_eta_);
	EventDataTree->Branch("normHEEPEles_SCposn_phi", &normHEEPEles_SCposn_phi_);
	EventDataTree->Branch("normHEEPEles_SC_rawEnergy", &normHEEPEles_SC_rawEnergy_);
	EventDataTree->Branch("normHEEPEles_SC_recHits_Et",  &normHEEPEles_SC_recHits_Et_);
	EventDataTree->Branch("normHEEPEles_SC_recHits_eta", &normHEEPEles_SC_recHits_eta_);
	EventDataTree->Branch("normHEEPEles_SC_recHits_phi", &normHEEPEles_SC_recHits_phi_);
	EventDataTree->Branch("normHEEPEles_SC_recHits_isFromEB", &normHEEPEles_SC_recHits_isFromEB_);
	EventDataTree->Branch("normHEEPEles_SC_totEnergyRecHits", &normHEEPEles_SC_totEnergyRecHits_);
	EventDataTree->Branch("normHEEPEles_SC_totNumRecHits", &normHEEPEles_SC_totNumRecHits_);

	// Track information for correcting tracker isolation variable
	EventDataTree->Branch("normHEEPEles_gsfTrk_eta", &normHEEPEles_gsfTrk_eta_);
	EventDataTree->Branch("normHEEPEles_gsfTrk_phi", &normHEEPEles_gsfTrk_phi_);
	EventDataTree->Branch("normHEEPEles_gsfTrk_vz", &normHEEPEles_gsfTrk_vz_);
	EventDataTree->Branch("normHEEPEles_innerIsoConeTrks_pt", &normHEEPEles_innerIsoConeTrks_pt_);
	EventDataTree->Branch("normHEEPEles_innerIsoConeTrks_eta", &normHEEPEles_innerIsoConeTrks_eta_);
	EventDataTree->Branch("normHEEPEles_innerIsoConeTrks_phi", &normHEEPEles_innerIsoConeTrks_phi_);
	EventDataTree->Branch("normHEEPEles_innerIsoConeTrks_vz", &normHEEPEles_innerIsoConeTrks_vz_);
}

void BstdZeeNTupler::SetupBstdEleBranches(){
	//Setting up the pointer links between the special reco'n GSF electron variables and the corresponding branches...
	EventDataTree->Branch("bstdGsfEles_number", &bstdGsfEles_number_, "bstdGsfEles_number/i"); //unsigned int
	EventDataTree->Branch("bstdGsfEles_p4ptr_", &bstdGsfEles_p4ptr_); //std::vector<ROOT::Math::XYZTVector>*
	bstdGsfEles_p4ptr_ = &bstdGsfEles_p4_;	
	EventDataTree->Branch("bstdGsfEles_charge", &bstdGsfEles_charge_); //std::vector<Int_t>
	
	EventDataTree->Branch("bstdGsfEles_Et", &bstdGsfEles_Et_);
	EventDataTree->Branch("bstdGsfEles_HEEP_Et", &bstdGsfEles_HEEP_Et_);
	EventDataTree->Branch("bstdGsfEles_Eta", &bstdGsfEles_Eta_);
	EventDataTree->Branch("bstdGsfEles_scEta", &bstdGsfEles_scEta_);
	EventDataTree->Branch("bstdGsfEles_ecalDriven", &bstdGsfEles_ecalDriven_);
	EventDataTree->Branch("bstdGsfEles_ecalDrivenSeed", &bstdGsfEles_ecalDrivenSeed_);
	
	EventDataTree->Branch("bstdGsfEles_HEEP_dEtaIn", &bstdGsfEles_HEEP_dEtaIn_);
	EventDataTree->Branch("bstdGsfEles_HEEP_dPhiIn", &bstdGsfEles_HEEP_dPhiIn_);
	EventDataTree->Branch("bstdGsfEles_HoverE", &bstdGsfEles_HoverE_);
	EventDataTree->Branch("bstdGsfEles_sigmaIetaIeta", &bstdGsfEles_sigmaIetaIeta_);
	EventDataTree->Branch("bstdGsfEles_scSigmaIetaIeta", &bstdGsfEles_scSigmaIetaIeta_);
	
	EventDataTree->Branch("bstdGsfEles_dr03EmIsoEt", &bstdGsfEles_dr03EmIsoEt_);
	EventDataTree->Branch("bstdGsfEles_dr03HadDepth1IsoEt", &bstdGsfEles_dr03HadDepth1IsoEt_);
	EventDataTree->Branch("bstdGsfEles_dr03HadDepth2IsoEt", &bstdGsfEles_dr03HadDepth2IsoEt_);
	EventDataTree->Branch("bstdGsfEles_dr03TkIsoPt", &bstdGsfEles_dr03TkIsoPt_);
	
	EventDataTree->Branch("bstdGsfEles_e2x5Max", &bstdGsfEles_e2x5Max_);
	EventDataTree->Branch("bstdGsfEles_e5x5", &bstdGsfEles_e5x5_);

	//HEEP analyser variables branches ...
	//kinematic and geometric methods
	EventDataTree->Branch("bstdHEEPEles_et",         &bstdHEEPEles_et_);
	EventDataTree->Branch("bstdHEEPEles_gsfEt",      &bstdHEEPEles_gsfEt_);
	EventDataTree->Branch("bstdHEEPEles_scEt",       &bstdHEEPEles_scEt_);
	EventDataTree->Branch("bstdHEEPEles_energy",     &bstdHEEPEles_energy_);
	EventDataTree->Branch("bstdHEEPEles_gsfEnergy",  &bstdHEEPEles_gsfEnergy_);
	EventDataTree->Branch("bstdHEEPEles_caloEnergy", &bstdHEEPEles_caloEnergy_);
	EventDataTree->Branch("bstdHEEPEles_eta",        &bstdHEEPEles_eta_);
	EventDataTree->Branch("bstdHEEPEles_scEta",      &bstdHEEPEles_scEta_);
	EventDataTree->Branch("bstdHEEPEles_detEta",     &bstdHEEPEles_detEta_);
	EventDataTree->Branch("bstdHEEPEles_detEtaAbs",  &bstdHEEPEles_detEtaAbs_);
	EventDataTree->Branch("bstdHEEPEles_phi",        &bstdHEEPEles_phi_);
	EventDataTree->Branch("bstdHEEPEles_scPhi",      &bstdHEEPEles_scPhi_);
	EventDataTree->Branch("bstdHEEPEles_detPhi",     &bstdHEEPEles_detPhi_);
	EventDataTree->Branch("bstdHEEPEles_zVtx",       &bstdHEEPEles_zVtx_);
	EventDataTree->Branch("bstdHEEPEles_p4ptr",      &bstdHEEPEles_p4ptr_);
	bstdHEEPEles_p4ptr_ = &bstdHEEPEles_p4_;
	EventDataTree->Branch("bstdHEEPEles_gsfP4ptr",   &bstdHEEPEles_gsfP4ptr_);
	bstdHEEPEles_gsfP4ptr_ = &bstdHEEPEles_gsfP4_;

	//classification (couldnt they have just named it 'type')
	EventDataTree->Branch("bstdHEEPEles_classification",  &bstdHEEPEles_classification_);
	EventDataTree->Branch("bstdHEEPEles_isEcalDriven",    &bstdHEEPEles_isEcalDriven_);
	EventDataTree->Branch("bstdHEEPEles_isTrackerDriven", &bstdHEEPEles_isTrackerDriven_);
	EventDataTree->Branch("bstdHEEPEles_isEB",            &bstdHEEPEles_isEB_);
	EventDataTree->Branch("bstdHEEPEles_isEE",            &bstdHEEPEles_isEE_);

	//track methods
	EventDataTree->Branch("bstdHEEPEles_charge",    &bstdHEEPEles_charge_);
	EventDataTree->Branch("bstdHEEPEles_trkCharge", &bstdHEEPEles_trkCharge_);
	EventDataTree->Branch("bstdHEEPEles_pVtx",      &bstdHEEPEles_pVtx_);
	EventDataTree->Branch("bstdHEEPEles_pCalo",     &bstdHEEPEles_pCalo_);
	EventDataTree->Branch("bstdHEEPEles_ptVtx",     &bstdHEEPEles_ptVtx_);
	EventDataTree->Branch("bstdHEEPEles_ptCalo",    &bstdHEEPEles_ptCalo_);

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	EventDataTree->Branch("bstdHEEPEles_hOverE",       &bstdHEEPEles_hOverE_);
	EventDataTree->Branch("bstdHEEPEles_dEtaIn",       &bstdHEEPEles_dEtaIn_);
	EventDataTree->Branch("bstdHEEPEles_dPhiIn",       &bstdHEEPEles_dPhiIn_);
	EventDataTree->Branch("bstdHEEPEles_dPhiOut",      &bstdHEEPEles_dPhiOut_);
	EventDataTree->Branch("bstdHEEPEles_epIn",         &bstdHEEPEles_epIn_);
	EventDataTree->Branch("bstdHEEPEles_epOut",        &bstdHEEPEles_epOut_);
	EventDataTree->Branch("bstdHEEPEles_fbrem",        &bstdHEEPEles_fbrem_);
	EventDataTree->Branch("bstdHEEPEles_bremFrac",     &bstdHEEPEles_bremFrac_);
	EventDataTree->Branch("bstdHEEPEles_invEOverInvP", &bstdHEEPEles_invEOverInvP_);

	//shower shape variables
	EventDataTree->Branch("bstdHEEPEles_sigmaEtaEta",       &bstdHEEPEles_sigmaEtaEta_);
	EventDataTree->Branch("bstdHEEPEles_sigmaEtaEtaUnCorr", &bstdHEEPEles_sigmaEtaEtaUnCorr_);
	EventDataTree->Branch("bstdHEEPEles_sigmaIEtaIEta",     &bstdHEEPEles_sigmaIEtaIEta_);
	EventDataTree->Branch("bstdHEEPEles_e1x5",              &bstdHEEPEles_e1x5_);
	EventDataTree->Branch("bstdHEEPEles_e2x5Max",           &bstdHEEPEles_e2x5Max_);
	EventDataTree->Branch("bstdHEEPEles_e5x5",              &bstdHEEPEles_e5x5_);
	EventDataTree->Branch("bstdHEEPEles_e1x5Over5x5",       &bstdHEEPEles_e1x5Over5x5_);
	EventDataTree->Branch("bstdHEEPEles_e2x5MaxOver5x5",    &bstdHEEPEles_e2x5MaxOver5x5_);

	//isolation, we use cone of 0.3
	EventDataTree->Branch("bstdHEEPEles_isolEm",            &bstdHEEPEles_isolEm_);
	EventDataTree->Branch("bstdHEEPEles_isolHad",           &bstdHEEPEles_isolHad_);
	EventDataTree->Branch("bstdHEEPEles_isolHadDepth1",     &bstdHEEPEles_isolHadDepth1_);
	EventDataTree->Branch("bstdHEEPEles_isolHadDepth2",     &bstdHEEPEles_isolHadDepth2_);
	EventDataTree->Branch("bstdHEEPEles_isolPtTrks",        &bstdHEEPEles_isolPtTrks_);
	EventDataTree->Branch("bstdHEEPEles_isolEmHadDepth1",   &bstdHEEPEles_isolEmHadDepth1_);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BstdZeeNTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BstdZeeNTupler);


//------------ method for clearing contents/setting default values of variables that should get new values in each event -------------
void
BstdZeeNTupler::ResetEventByEventVariables(){
	
	//Resetting the values of the event information variables...
	evt_runNum_  = 0;
	evt_lumiSec_ = 0;
	evt_evtNum_  = 0;
	
	// Resetting the values of the varaibles counting the number of Z bosons per event ...
	mc_numZs_status2_ = 0;
	mc_numZs_status3_ = 0;
	mc_numZs_statusNot2Or3_ = 0;

	mcZboson_pdgId_  = -999; // Int_t
	mcZboson_status_ = -999; // Int_t
	mcZboson_p4_.SetPxPyPzE(99999.9,0.0,0.0,99999.9);
	mcZboson_status2_p4_.SetPxPyPzE(99999.9,0.0,0.0,99999.9);
	mcZboson_numDaughters_ = -999;
	mcZboson_daughters_dR_   = -999.9;
	mcZboson_daughters_dEta_ = -999.9;
	mcZboson_daughters_dPhi_ = -999.9;
	mcZboson_daughters_openingAngle_  = -999.9;
	mcZboson_daughterA_p4_.SetPxPyPzE(99999.9,0.0,0.0,99999.9);
	mcZboson_daughterB_p4_.SetPxPyPzE(99999.9,0.0,0.0,99999.9);

	//Clearing contents of standard GSF electron vectors...
	normGsfEles_number_ = 0;
	normGsfEles_p4_.clear();
	normGsfEles_charge_.clear();

	normGsfEles_Et_.clear();
	normGsfEles_HEEP_Et_.clear();
	normGsfEles_Eta_.clear();
	normGsfEles_scEta_.clear();
	normGsfEles_ecalDriven_.clear();
	normGsfEles_ecalDrivenSeed_.clear();

	normGsfEles_HEEP_dEtaIn_.clear();
	normGsfEles_HEEP_dPhiIn_.clear();
	normGsfEles_HoverE_.clear();
	normGsfEles_sigmaIetaIeta_.clear();
	normGsfEles_scSigmaIetaIeta_.clear();

	normGsfEles_dr03EmIsoEt_.clear();
	normGsfEles_dr03HadDepth1IsoEt_.clear();
	normGsfEles_dr03HadDepth2IsoEt_.clear();
	normGsfEles_dr03TkIsoPt_.clear();

	normGsfEles_e2x5Max_.clear();
	normGsfEles_e5x5_.clear();
	
	//normGsfEles_tswHEEPEle_.clear();
	//Clearing the contents of the heep::Ele variable vectors ...
	//kinematic and geometric methods
	normHEEPEles_et_.clear();
	normHEEPEles_gsfEt_.clear();
	normHEEPEles_scEt_.clear();
	normHEEPEles_energy_.clear();
	normHEEPEles_gsfEnergy_.clear();
	normHEEPEles_caloEnergy_.clear();
	normHEEPEles_ecalEnergyError_.clear();
	normHEEPEles_eta_.clear();
	normHEEPEles_scEta_.clear();
	normHEEPEles_detEta_.clear();
	normHEEPEles_detEtaAbs_.clear();
	normHEEPEles_phi_.clear();
	normHEEPEles_scPhi_.clear();
	normHEEPEles_detPhi_.clear();
	normHEEPEles_zVtx_.clear();
	normHEEPEles_p4_.clear();
	normHEEPEles_gsfP4_.clear();

	//classification (couldnt they have just named it 'type')
	 normHEEPEles_classification_.clear();
	normHEEPEles_isEcalDriven_.clear();
	normHEEPEles_isTrackerDriven_.clear();
	normHEEPEles_isEB_.clear();
	normHEEPEles_isEE_.clear();

	//track methods
	normHEEPEles_charge_.clear();
	normHEEPEles_trkCharge_.clear();
	normHEEPEles_pVtx_.clear();
	normHEEPEles_pCalo_.clear();
	normHEEPEles_ptVtx_.clear();
	normHEEPEles_ptCalo_.clear();
	normHEEPEles_closestCtfTrk_pt_.clear();
	normHEEPEles_closestCtfTrk_eta_.clear();
	normHEEPEles_closestCtfTrk_phi_.clear();
	normHEEPEles_closestCtfTrk_innerPt_.clear();
	normHEEPEles_closestCtfTrk_innerEta_.clear();
	normHEEPEles_closestCtfTrk_innerPhi_.clear();
	normHEEPEles_closestCtfTrk_outerPt_.clear();
	normHEEPEles_closestCtfTrk_outerEta_.clear();
	normHEEPEles_closestCtfTrk_outerPhi_.clear();

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	normHEEPEles_hOverE_.clear();
	normHEEPEles_dEtaIn_.clear();
	normHEEPEles_dPhiIn_.clear();
	normHEEPEles_dPhiOut_.clear();
	normHEEPEles_epIn_.clear();
	normHEEPEles_epOut_.clear();
	normHEEPEles_fbrem_.clear();
	normHEEPEles_bremFrac_.clear();
	normHEEPEles_invEOverInvP_.clear();

	//shower shape variables
	normHEEPEles_sigmaEtaEta_.clear();
	normHEEPEles_sigmaEtaEtaUnCorr_.clear();
	normHEEPEles_sigmaIEtaIEta_.clear();
	normHEEPEles_e1x5_.clear();
	normHEEPEles_e2x5Max_.clear();
	normHEEPEles_e5x5_.clear();
	normHEEPEles_e1x5Over5x5_.clear();
	normHEEPEles_e2x5MaxOver5x5_.clear();

	//isolation, we use cone of 0.3
	normHEEPEles_isolEm_.clear();
	normHEEPEles_isolHad_.clear();
	normHEEPEles_isolHadDepth1_.clear();
	normHEEPEles_isolHadDepth2_.clear();
	normHEEPEles_isolPtTrks_.clear();
	normHEEPEles_isolEmHadDepth1_.clear();

	normHEEPEles_numMissInnerHits_.clear();

	normHEEPEles_SCposn_eta_.clear();
	normHEEPEles_SCposn_phi_.clear();
	normHEEPEles_SC_rawEnergy_.clear();
	normHEEPEles_SC_recHits_Et_.clear();
	normHEEPEles_SC_recHits_eta_.clear();
	normHEEPEles_SC_recHits_phi_.clear();
	normHEEPEles_SC_recHits_isFromEB_.clear();
	normHEEPEles_SC_totEnergyRecHits_.clear();
	normHEEPEles_SC_totNumRecHits_.clear();

	normHEEPEles_gsfTrk_eta_.clear();
	normHEEPEles_gsfTrk_phi_.clear();
	normHEEPEles_gsfTrk_vz_.clear();
	normHEEPEles_innerIsoConeTrks_pt_.clear();
	normHEEPEles_innerIsoConeTrks_eta_.clear();
	normHEEPEles_innerIsoConeTrks_phi_.clear();
	normHEEPEles_innerIsoConeTrks_vz_.clear();

	///////////////////////////////////////////////////////
	//Clearing contents of standard GSF electron vectors...
	bstdGsfEles_number_ = 0;
	bstdGsfEles_p4_.clear();
	bstdGsfEles_charge_.clear();

	bstdGsfEles_Et_.clear();
	bstdGsfEles_HEEP_Et_.clear();
	bstdGsfEles_Eta_.clear();
	bstdGsfEles_scEta_.clear();
	bstdGsfEles_ecalDriven_.clear();
	bstdGsfEles_ecalDrivenSeed_.clear();

	bstdGsfEles_HEEP_dEtaIn_.clear();
	bstdGsfEles_HEEP_dPhiIn_.clear();
	bstdGsfEles_HoverE_.clear();
	bstdGsfEles_sigmaIetaIeta_.clear();
	bstdGsfEles_scSigmaIetaIeta_.clear();

	bstdGsfEles_dr03EmIsoEt_.clear();
	bstdGsfEles_dr03HadDepth1IsoEt_.clear();
	bstdGsfEles_dr03HadDepth2IsoEt_.clear();
	bstdGsfEles_dr03TkIsoPt_.clear();

	bstdGsfEles_e2x5Max_.clear();
	bstdGsfEles_e5x5_.clear();

	//normGsfEles_tswHEEPEle_.clear();
	//Clearing the contents of the heep::Ele variable vectors ...
	//kinematic and geometric methods
	bstdHEEPEles_et_.clear();
	bstdHEEPEles_gsfEt_.clear();
	bstdHEEPEles_scEt_.clear();
	bstdHEEPEles_energy_.clear();
	bstdHEEPEles_gsfEnergy_.clear();
	bstdHEEPEles_caloEnergy_.clear();
	bstdHEEPEles_eta_.clear();
	bstdHEEPEles_scEta_.clear();
	bstdHEEPEles_detEta_.clear();
	bstdHEEPEles_detEtaAbs_.clear();
	bstdHEEPEles_phi_.clear();
	bstdHEEPEles_scPhi_.clear();
	bstdHEEPEles_detPhi_.clear();
	bstdHEEPEles_zVtx_.clear();
	bstdHEEPEles_p4_.clear();
	bstdHEEPEles_gsfP4_.clear();

	//classification (couldnt they have just named it 'type')
	 bstdHEEPEles_classification_.clear();
	bstdHEEPEles_isEcalDriven_.clear();
	bstdHEEPEles_isTrackerDriven_.clear();
	bstdHEEPEles_isEB_.clear();
	bstdHEEPEles_isEE_.clear();

	//track methods
	bstdHEEPEles_charge_.clear();
	bstdHEEPEles_trkCharge_.clear();
	bstdHEEPEles_pVtx_.clear();
	bstdHEEPEles_pCalo_.clear();
	bstdHEEPEles_ptVtx_.clear();
	bstdHEEPEles_ptCalo_.clear();

	//abreviations of overly long GsfElectron methods, I'm sorry but if you cant figure out what hOverE() means, you shouldnt be using this class
	bstdHEEPEles_hOverE_.clear();
	bstdHEEPEles_dEtaIn_.clear();
	bstdHEEPEles_dPhiIn_.clear();
	bstdHEEPEles_dPhiOut_.clear();
	bstdHEEPEles_epIn_.clear();
	bstdHEEPEles_epOut_.clear();
	bstdHEEPEles_fbrem_.clear();
	bstdHEEPEles_bremFrac_.clear();
	bstdHEEPEles_invEOverInvP_.clear();

	//shower shape variables
	bstdHEEPEles_sigmaEtaEta_.clear();
	bstdHEEPEles_sigmaEtaEtaUnCorr_.clear();
	bstdHEEPEles_sigmaIEtaIEta_.clear();
	bstdHEEPEles_e1x5_.clear();
	bstdHEEPEles_e2x5Max_.clear();
	bstdHEEPEles_e5x5_.clear();
	bstdHEEPEles_e1x5Over5x5_.clear();
	bstdHEEPEles_e2x5MaxOver5x5_.clear();

	//isolation, we use cone of 0.3
	bstdHEEPEles_isolEm_.clear();
	bstdHEEPEles_isolHad_.clear();
	bstdHEEPEles_isolHadDepth1_.clear();
	bstdHEEPEles_isolHadDepth2_.clear();
	bstdHEEPEles_isolPtTrks_.clear();
	bstdHEEPEles_isolEmHadDepth1_.clear();
}

//------------ method for accessing the trigger information for each event -------------
void
BstdZeeNTupler::accessTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool beVerbose){
	bool hltBit(false);
	// get HLT results
  	edm::Handle<edm::TriggerResults> HLTR;
  	iEvent.getByLabel(hltResultsTag_, HLTR);

  	if (HLTR.isValid()) {
   	// triggerIndex must be less than the size of HLTR or you get a CMSException: _M_range_check
    	if (hltPathIndex_ < HLTR->size())
			hltBit = HLTR->accept(hltPathIndex_); 
  	}
	hltPathADecision_ = hltBit;

	// Get the e-mu trigger information, and put it in event object
	bool trg_emuPath_decision = false;
  	if (HLTR.isValid()) {
   	// triggerIndex must be less than the size of HLTR or you get a CMSException: _M_range_check
    	if (trg_emuPath_idx_ < HLTR->size())
    		trg_emuPath_decision = HLTR->accept(trg_emuPath_idx_);
  	}
  	event_->SetEMuTriggerInfo(trg_emuPath_name_, trg_emuPath_decision);

  	if(vBool_){
		if(trg_emuPath_decision)
			std::cout << " ->Event passed e-mu trigger (i.e. " << trg_emuPath_name_ << "=" << hltConfig_.triggerName(trg_emuPath_idx_) << "?)" << std::endl;
		else
			std::cout << " ->Event did *NOT* pass e-mu trigger (i.e. " << trg_emuPath_name_ << "=" << hltConfig_.triggerName(trg_emuPath_idx_) << "?)" << std::endl;
	}

	//----------------------------------------------
	//-- Now, looking at trigger object info.
  	// Getting a handle to the TriggerEvent object ...
  	edm::Handle<trigger::TriggerEvent> hltEventH;
  	iEvent.getByLabel(hltEventTag_, hltEventH);

	//Get the Et values for the objects passing the last filter of this path ...
	hltPathA_highestTrigObjEt_ = -999.9;
	//int nrJet300=0;

	//it is important to specify the right HLT process for the filter, not doing this is a common bug
	trigger::size_type filterIndex = hltEventH->filterIndex(edm::InputTag(hltPathA_nameOfLastFilter_,"",hltEventTag_.process()));
	if(filterIndex<hltEventH->sizeFilters()){
		if(beVerbose){std::cout << " ->Accessing the objects that passed the filter " <<  hltPathA_nameOfLastFilter_ << std::endl;}
	   const trigger::Keys& trigKeys = hltEventH->filterKeys(filterIndex);
	   const trigger::TriggerObjectCollection & trigObjColl(hltEventH->getObjects());
	   //now loop of the trigger objects passing filter
	   for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){
	      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	      if(beVerbose){std::cout << "    Object's Et = " << obj.et();}
	      if(obj.et()>hltPathA_highestTrigObjEt_){
	      	hltPathA_highestTrigObjEt_ = obj.et();
	      }
	      if(beVerbose){std::cout << "; highest obj Et = " << hltPathA_highestTrigObjEt_ << " now." << std::endl;}
	   }
	}//end filter size check
	else{
		if(beVerbose){std::cout << " ->FAILLED basic if statement => CANNOT access the objects that passed the filter " <<  hltPathA_nameOfLastFilter_ << std::endl;
			std::cout << "     [This if fine iff the event failled hltPathA.]" << std::endl;}
	}

	//Now fill up the histogram of the highest Et trig obj Et values ...
	if(hltPathA_highestTrigObjEt_>0.0){
		if(beVerbose){std::cout << "       => Et of highest Et obj was: " << hltPathA_highestTrigObjEt_ << std::endl;
			std::cout << "          Filling histogram with this value ..." << std::endl;}
		highestEtTrigObjPathA_EtHist_->Fill(hltPathA_highestTrigObjEt_);
	}
	else{
		if(beVerbose){std::cout << "   There were no trigger objects passing the filter => Et of highest Et object histogram was filled with -999.9." << std::endl;}
		highestEtTrigObjPathA_EtHist_->Fill(-999.9);
	}
}

//------------ method for printing out all of the trigger names... -------------
void
BstdZeeNTupler::printDatasetsAndTriggerNames(){
	const std::vector<std::string> DatasetNames = hltConfig_.datasetNames();
	const std::vector< std::vector<std::string> > TriggerNames = hltConfig_.datasetContents();
	int numTriggerNames = 0;

	std::cout << "****************" << std::endl;
	std::cout << "The names of the " << DatasetNames.size() << " datasets in this run are:" << std::endl;
	for(unsigned int idx_Name = 0; idx_Name < DatasetNames.size(); idx_Name++)
	{
		std::cout << "    " << DatasetNames.at(idx_Name) << std::endl;
	}

	for(unsigned int i = 0; i < TriggerNames.size(); i++)
		numTriggerNames += TriggerNames.at(i).size();
	std::cout << "The names of the " << numTriggerNames << " triggers in this run are:" << std::endl;
	for(unsigned int idxA = 0; idxA < TriggerNames.size(); idxA++)
	{
		for(unsigned int idxB = 0; idxB < TriggerNames.at(idxA).size(); idxB++)
			std::cout << "     " << TriggerNames.at(idxA).at(idxB) << std::endl;
	}
	std::cout << std::endl << "...or from triggerNames() method " << ", the names of the " << hltConfig_.triggerNames().size() << " triggers are:" << std::endl;
	for(unsigned int idx = 0; idx < hltConfig_.triggerNames().size(); idx++)
		std::cout << "     " << hltConfig_.triggerNames().at(idx) << std::endl;
	std::cout << std::endl;

//	// Now, printing out the names of the modules for the trigger ...
//	const std::vector<std::string> hltPathA_moduleNames = hltConfig_.moduleLabels(hltPathA_);
//	std::cout << "The modules of the " << hltPathA_ << " trigger path are:" << std::endl;
//	for(unsigned int i = 0; i<hltPathA_moduleNames.size(); i++)
//		std::cout << "    " << hltPathA_moduleNames.at(i) << std::endl;

	std::cout << "****************" << std::endl;

}

//----------------------------------------------------------------------------------------------------
//------------ method for reading in the event information (run no., lumi sec etc...)  ---------------
void BstdZeeNTupler::ReadInEvtInfo(bool beVerbose, const edm::Event& edmEventObject){

	evt_runNum_  = edmEventObject.id().run();
	evt_lumiSec_ = edmEventObject.id().luminosityBlock();
	evt_evtNum_  = edmEventObject.id().event();
	
	event_->SetBasicEventInformation(evt_runNum_, evt_lumiSec_, evt_evtNum_);

	if(beVerbose){
		std::cout << " ->Run: " << evt_runNum_ << ", Lumi block: " << evt_lumiSec_ << ", Evt: " << evt_evtNum_ << std::endl;
	}

	// Read the vertex and pile up re-weighting info into branches
	ReadInVtxAndPUInfo(beVerbose, edmEventObject);
}

//----------------------------------------------------------------------------------------------------
//------------ method for reading in the vertex information, ...                       ---------------
//------------ and calculating the pile up re-weighting factors                        ---------------
void BstdZeeNTupler::ReadInVtxAndPUInfo(bool beVerbose, const edm::Event& edmEventRef)
{
	// 1) Read in the info about the MC PU vertices, and store PU re-weighting factor [IFF running over MC ...]
	if(mcFlag_){
		edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
		edmEventRef.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			int BX = PVI->getBunchCrossing();
			if(BX == 0) {
				float trueNumPUVertices = PVI->getTrueNumInteractions();
				unsigned int numPUVertices = PVI->getPU_NumInteractions();
				std::vector<float> puVertices_zPosns = PVI->getPU_zpositions();
				double myWeight = 1.0;
				if(calcPileUpWeights_)
					myWeight = lumiWeights_.weight( trueNumPUVertices );

				event_->SetMCPUInfo(trueNumPUVertices, numPUVertices, puVertices_zPosns, myWeight);
			}
		}
	}

	// 2) Read in the info about the reconstructed vertices

	edm::Handle< std::vector<reco::Vertex> > handle_RecoVtxs;
	edmEventRef.getByLabel(vertexSrc_, handle_RecoVtxs);

	unsigned int numVtxs = 0;
	unsigned int numGoodVtxs = 0;
	// Count total num of vertices, and num of good vertices
	for (std::vector<reco::Vertex>::const_iterator itv = handle_RecoVtxs->begin(); itv != handle_RecoVtxs->end(); ++itv) {
		numVtxs += 1;
		// require that the vertex meets certain criteria
		if(itv->ndof()<5) continue;
		if(fabs(itv->z())>50.0) continue;
		if(fabs(itv->position().rho())>2.0) continue;
		numGoodVtxs += 1;
	}
	event_->SetRecoVtxInfo(numVtxs, numGoodVtxs);

	// Print out this information if in 'verbose' mode ...
	if(beVerbose)
		event_->PrintPUVtxInfo();

	// Store the pointer to the first (i.e. highest sum pt) vertex, for future use ...
	recoVtx_highestSumPtVtx_ = handle_RecoVtxs->front();

	// 3) Read in the 'rho' jet recon 'PU activity' parameter - to be used for PU correction of calo isol values
	edm::Handle<double> rhoHandle;
	edmEventRef.getByLabel(edm::InputTag("kt6PFJets:rho"),rhoHandle);
	const double rho = *(rhoHandle.product());
	event_->SetPURho(rho);
	if(beVerbose)
		event_->PrintPURho();
}

void BstdZeeNTupler::ReadInGenInfo(bool beVerbose, const edm::Event& edmEventRef)
{
	// Read in the mc Event weight from the GenEventInfoProduct ...
	edm::Handle<GenEventInfoProduct> handle_GenInfoProduct;
	edmEventRef.getByLabel("generator", handle_GenInfoProduct);

	event_->SetMCGenWeight(handle_GenInfoProduct->weight());

	if(beVerbose)
		event_->PrintMCGenWeight();

	// Read in the LHE Z boson 4-momentum (iff the file contains LHE info) ...
	edm::Handle<LHEEventProduct> h_lheEvent;
	edmEventRef.getByLabel("source",h_lheEvent);
	if( !(h_lheEvent.failedToGet()) ){
		const lhef::HEPEUP eventInfo = h_lheEvent->hepeup();

		// Now, loop over the particles in the LHE event ...
		int numParticles = eventInfo.NUP;
		for(int idx=0; idx<numParticles; idx++){
			// Check if the particle is a Z boson; if so, then extract it's 4 momentum & read into tsw::Event
			if(eventInfo.IDUP.at(idx)==23){
				lhef::HEPEUP::FiveVector p5 =  eventInfo.PUP.at(idx);

				ROOT::Math::XYZTVector Zboson_p4;
				Zboson_p4.SetPxPyPzE(p5[0], p5[1], p5[2], p5[3]);
				event_->SetLHEZbosonInfo(Zboson_p4);
			}//end: if PID==23
		}//END: for idx<numParticles
		if(beVerbose)
			event_->PrintLHEZbosonInfo();
	}
	else if(beVerbose)
		std::cout << " * Couldn't find any LHE information *" << std::endl;


}


//------------ method for reading in the values of the standard GSF electron variables ---------------
void BstdZeeNTupler::ReadInNormGsfEles(bool beVerbose, const edm::Handle<reco::GsfElectronCollection>& handle_normGsfEles, edm::ESHandle<CaloGeometry>& handle_caloGeom, EcalRecHitMetaCollection* recHitsMeta_EB, EcalRecHitMetaCollection* recHitsMeta_EE, edm::Handle<reco::TrackCollection>& handle_ctfTracks, const edm::Event& iEvent){
	
	// Setting up value map handles for electron isol values re-calculated from IsoDeposits
	edm::Handle< edm::ValueMap<double> > h_stdEleIsoFromDeps_Tk;
	edm::Handle< edm::ValueMap<double> > h_modEleIsoFromDeps_Tk;
	iEvent.getByLabel("stdEleIsoFromDepsTk",h_stdEleIsoFromDeps_Tk);
	iEvent.getByLabel("modEleIsoFromDepsTk",h_modEleIsoFromDeps_Tk);

	edm::Handle< edm::ValueMap<double> > h_stdEleIsoFromDeps_Ecal;
	edm::Handle< edm::ValueMap<double> > h_modEleIsoFromDeps_Ecal;
	iEvent.getByLabel("stdEleIsoFromDepsEcal",h_stdEleIsoFromDeps_Ecal);
	iEvent.getByLabel("modEleIsoFromDepsEcal",h_modEleIsoFromDeps_Ecal);

	edm::Handle< edm::ValueMap<double> > h_stdEleIsoFromDeps_HcalD1;
	edm::Handle< edm::ValueMap<double> > h_modEleIsoFromDeps_HcalD1;
	iEvent.getByLabel("stdEleIsoFromDepsHcalD1",h_stdEleIsoFromDeps_HcalD1);
	iEvent.getByLabel("modEleIsoFromDepsHcalD1",h_modEleIsoFromDeps_HcalD1);

	edm::Handle< edm::ValueMap<double> > h_inrXSVetoEleIso_Tk;
	edm::Handle< edm::ValueMap<double> > h_inrSVetoEleIso_Tk;
	edm::Handle< edm::ValueMap<double> > h_inrMVetoEleIso_Tk;
	edm::Handle< edm::ValueMap<double> > h_inrLVetoEleIso_Tk;
	edm::Handle< edm::ValueMap<double> > h_inrXLVetoEleIso_Tk;
	iEvent.getByLabel("innerXSVetoModEleIso","track", h_inrXSVetoEleIso_Tk);
	iEvent.getByLabel("innerSVetoModEleIso","track",  h_inrSVetoEleIso_Tk);
	iEvent.getByLabel("innerMVetoModEleIso","track",  h_inrMVetoEleIso_Tk);
	iEvent.getByLabel("innerLVetoModEleIso","track",  h_inrLVetoEleIso_Tk);
	iEvent.getByLabel("innerXLVetoModEleIso","track", h_inrXLVetoEleIso_Tk);

	edm::Handle< edm::ValueMap<double> > h_inrXSVetoEleIso_Ecal;
	edm::Handle< edm::ValueMap<double> > h_inrSVetoEleIso_Ecal;
	edm::Handle< edm::ValueMap<double> > h_inrMVetoEleIso_Ecal;
	edm::Handle< edm::ValueMap<double> > h_inrLVetoEleIso_Ecal;
	edm::Handle< edm::ValueMap<double> > h_inrXLVetoEleIso_Ecal;
	iEvent.getByLabel("innerXSVetoModEleIso","ecal", h_inrXSVetoEleIso_Ecal);
	iEvent.getByLabel("innerSVetoModEleIso","ecal",  h_inrSVetoEleIso_Ecal);
	iEvent.getByLabel("innerMVetoModEleIso","ecal",  h_inrMVetoEleIso_Ecal);
	iEvent.getByLabel("innerLVetoModEleIso","ecal",  h_inrLVetoEleIso_Ecal);
	iEvent.getByLabel("innerXLVetoModEleIso","ecal", h_inrXLVetoEleIso_Ecal);

	edm::Handle< edm::ValueMap<double> > h_inrXSVetoEleIso_HcalD1;
	edm::Handle< edm::ValueMap<double> > h_inrSVetoEleIso_HcalD1;
	edm::Handle< edm::ValueMap<double> > h_inrMVetoEleIso_HcalD1;
	edm::Handle< edm::ValueMap<double> > h_inrLVetoEleIso_HcalD1;
	edm::Handle< edm::ValueMap<double> > h_inrXLVetoEleIso_HcalD1;
	iEvent.getByLabel("innerXSVetoModEleIso","hcalDepth1", h_inrXSVetoEleIso_HcalD1);
	iEvent.getByLabel("innerSVetoModEleIso","hcalDepth1",  h_inrSVetoEleIso_HcalD1);
	iEvent.getByLabel("innerMVetoModEleIso","hcalDepth1",  h_inrMVetoEleIso_HcalD1);
	iEvent.getByLabel("innerLVetoModEleIso","hcalDepth1",  h_inrLVetoEleIso_HcalD1);
	iEvent.getByLabel("innerXLVetoModEleIso","hcalDepth1", h_inrXLVetoEleIso_HcalD1);


	std::vector<std::string> phantomEleValMapLabels;
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr005To010");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr010To015");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr015To020");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr020To025");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr025To030");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr030To035");
	phantomEleValMapLabels.push_back("vetoAreaModEleIsoPhantomDr035To040");
	size_t nrPhantomEleValMaps = phantomEleValMapLabels.size();
	std::vector< edm::Handle< edm::ValueMap<double> > > hVec_inrVetoModIsosPhantomEle_dEta(nrPhantomEleValMaps),
				hVec_inrVetoModIsosPhantomEle_dPhi(nrPhantomEleValMaps),
				hVec_inrVetoModTkIsosPhantomEle(nrPhantomEleValMaps),
				hVec_inrVetoModEmIsosPhantomEle(nrPhantomEleValMaps),
				hVec_inrVetoModH1IsosPhantomEle(nrPhantomEleValMaps);
	for(size_t i=0; i<nrPhantomEleValMaps; i++){
		const std::string& phantomEleValMapLabel = phantomEleValMapLabels.at(i);
		iEvent.getByLabel(phantomEleValMapLabel, "dEtaPhantomEle", hVec_inrVetoModIsosPhantomEle_dEta.at(i));
		iEvent.getByLabel(phantomEleValMapLabel, "dPhiPhantomEle", hVec_inrVetoModIsosPhantomEle_dPhi.at(i));
		iEvent.getByLabel(phantomEleValMapLabel, "track",          hVec_inrVetoModTkIsosPhantomEle.at(i));
		iEvent.getByLabel(phantomEleValMapLabel, "ecal",           hVec_inrVetoModEmIsosPhantomEle.at(i));
		iEvent.getByLabel(phantomEleValMapLabel, "hcalDepth1",     hVec_inrVetoModH1IsosPhantomEle.at(i));
	}

	edm::Handle<reco::VertexCollection> h_vertices;
	iEvent.getByLabel("offlinePrimaryVertices", h_vertices);
	if( ! h_vertices.isValid() )
		edm::LogWarning("BstdZeeNTupler") << "WARNING!! Error will occur soon since handle to primary vertices is not valid!\n";
	const reco::Vertex& firstPrimaryVertex = h_vertices->front() ;

	//Setting the values of the standard reconstruction GSF electron variables...
	normGsfEles_number_ = handle_normGsfEles.product()->size();
	if(beVerbose){std::cout << " ->There are " << normGsfEles_number_ << " standard GSF electrons in this event."<< std::endl;}
	
	for(unsigned int eleIdx = 0; eleIdx < normGsfEles_number_; eleIdx++){
		
		reco::GsfElectron ithGsfEle = handle_normGsfEles.product()->at(eleIdx);
		heep::Ele ithHEEPEle(ithGsfEle);

		normGsfEles_p4ptr_->push_back(            ithGsfEle.p4()     );
		normGsfEles_charge_.push_back(            ithGsfEle.charge() );
 	
		normGsfEles_Et_.push_back(                ithGsfEle.et() );
		normGsfEles_HEEP_Et_.push_back(           ithGsfEle.caloEnergy()*sin(ithGsfEle.p4().theta()) );
		normGsfEles_Eta_.push_back(               ithGsfEle.eta()              );
		normGsfEles_scEta_.push_back(             ithGsfEle.caloPosition().eta() );

		normGsfEles_ecalDriven_.push_back(        ithGsfEle.ecalDriven()       );
		normGsfEles_ecalDrivenSeed_.push_back(    ithGsfEle.ecalDrivenSeed()   );

		normGsfEles_HEEP_dEtaIn_.push_back(       ithGsfEle.deltaEtaSuperClusterTrackAtVtx() );
		normGsfEles_HEEP_dPhiIn_.push_back(       ithGsfEle.deltaPhiSuperClusterTrackAtVtx() );
		normGsfEles_HoverE_.push_back(            ithGsfEle.hadronicOverEm()  );
		normGsfEles_sigmaIetaIeta_.push_back(     ithGsfEle.sigmaIetaIeta()   );
		normGsfEles_scSigmaIetaIeta_.push_back(   ithGsfEle.scSigmaIEtaIEta() );

		normGsfEles_dr03EmIsoEt_.push_back(        ithGsfEle.dr03EcalRecHitSumEt()      );
		normGsfEles_dr03HadDepth1IsoEt_.push_back( ithGsfEle.dr03HcalDepth1TowerSumEt() );
		normGsfEles_dr03HadDepth2IsoEt_.push_back( ithGsfEle.dr03HcalDepth2TowerSumEt() );
		normGsfEles_dr03TkIsoPt_.push_back(        ithGsfEle.dr03TkSumPt()              );

		normGsfEles_e2x5Max_.push_back(            ithGsfEle.e2x5Max() );
		normGsfEles_e5x5_.push_back(               ithGsfEle.e5x5()    );

		//Variables storing the heep::Ele method values ...
		normHEEPEles_et_.push_back(        ithHEEPEle.et() );
	   normHEEPEles_gsfEt_.push_back(     ithHEEPEle.gsfEt() );
		normHEEPEles_scEt_.push_back(      ithHEEPEle.scEt() );
		normHEEPEles_energy_.push_back(    ithHEEPEle.energy() );
		normHEEPEles_gsfEnergy_.push_back( ithHEEPEle.gsfEnergy() );
		normHEEPEles_caloEnergy_.push_back(ithHEEPEle.caloEnergy() );
		normHEEPEles_ecalEnergyError_.push_back(ithHEEPEle.gsfEle().ecalEnergyError() );
		normHEEPEles_eta_.push_back(       ithHEEPEle.eta() );
		normHEEPEles_scEta_.push_back(     ithHEEPEle.scEta() );
		normHEEPEles_detEta_.push_back(    ithHEEPEle.detEta() );
		normHEEPEles_detEtaAbs_.push_back( ithHEEPEle.detEtaAbs() );
		normHEEPEles_phi_.push_back(       ithHEEPEle.phi() );
		normHEEPEles_scPhi_.push_back(     ithHEEPEle.scPhi() );
		normHEEPEles_detPhi_.push_back(    ithHEEPEle.detPhi() );
		normHEEPEles_zVtx_.push_back(      ithHEEPEle.zVtx() );
		normHEEPEles_p4_.push_back(        ithHEEPEle.p4() );
		normHEEPEles_gsfP4_.push_back(     ithHEEPEle.gsfP4() );

		//Variables storing the heep::Ele method values - Classification...
		normHEEPEles_classification_.push_back( ithHEEPEle.classification() );
		normHEEPEles_isEcalDriven_.push_back(   ithHEEPEle.isEcalDriven() );
		normHEEPEles_isTrackerDriven_.push_back(ithHEEPEle.isTrackerDriven() );
		normHEEPEles_isEB_.push_back(           ithHEEPEle.isEB() );
		normHEEPEles_isEE_.push_back(           ithHEEPEle.isEE() );

		//Variables storing the heep::Ele method values - track methods ...
		normHEEPEles_charge_.push_back(   ithHEEPEle.charge() );
		normHEEPEles_trkCharge_.push_back(ithHEEPEle.trkCharge() );
		normHEEPEles_pVtx_.push_back(     ithHEEPEle.pVtx() );
		normHEEPEles_pCalo_.push_back(    ithHEEPEle.pCalo() );
		normHEEPEles_ptVtx_.push_back(    ithHEEPEle.ptVtx() );
		normHEEPEles_ptCalo_.push_back(   ithHEEPEle.ptCalo() );
		if(ithHEEPEle.gsfEle().closestCtfTrackRef().get()==0){
			if(beVerbose){std::cout << std::endl << "    ***** closestCtfTrackRef is NULL *****" << std::endl << std::endl;}
			normHEEPEles_closestCtfTrk_pt_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_eta_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_phi_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_innerPt_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_innerEta_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_innerPhi_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_outerPt_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_outerEta_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_outerPhi_.push_back( -999.9 );
		}
		else{
			normHEEPEles_closestCtfTrk_pt_.push_back(  ithHEEPEle.gsfEle().closestCtfTrackRef().get()->pt() );
			normHEEPEles_closestCtfTrk_eta_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->eta() );
			normHEEPEles_closestCtfTrk_phi_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->phi() );
//			normHEEPEles_closestCtfTrk_innerPt_.push_back(  ithHEEPEle.gsfEle().closestCtfTrackRef().get()->innerMomentum().Rho() );
//			normHEEPEles_closestCtfTrk_innerEta_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->innerMomentum().Eta() );
//			normHEEPEles_closestCtfTrk_innerPhi_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->innerMomentum().Phi() );
//			normHEEPEles_closestCtfTrk_outerPt_.push_back(  ithHEEPEle.gsfEle().closestCtfTrackRef().get()->outerPt() );
//			normHEEPEles_closestCtfTrk_outerEta_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->outerEta() );
//			normHEEPEles_closestCtfTrk_outerPhi_.push_back( ithHEEPEle.gsfEle().closestCtfTrackRef().get()->outerPhi() );
			// Had to switch to just using pt eta and phi methods of closest CTF track -  and filling other branches with dummy values for the moment - as am currently running over AOD
			normHEEPEles_closestCtfTrk_innerPt_.push_back(  -999.9 );
			normHEEPEles_closestCtfTrk_innerEta_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_innerPhi_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_outerPt_.push_back(  -999.9 );
			normHEEPEles_closestCtfTrk_outerEta_.push_back( -999.9 );
			normHEEPEles_closestCtfTrk_outerPhi_.push_back( -999.9 );
		}

		//Variables storing the heep::Ele method values ...
		normHEEPEles_hOverE_.push_back(      ithHEEPEle.hOverE() );
		normHEEPEles_dEtaIn_.push_back(      ithHEEPEle.dEtaIn() );
		normHEEPEles_dPhiIn_.push_back(      ithHEEPEle.dPhiIn() );
		normHEEPEles_dPhiOut_.push_back(     ithHEEPEle.dPhiOut() );
		normHEEPEles_epIn_.push_back(        ithHEEPEle.epIn() );
		normHEEPEles_epOut_.push_back(       ithHEEPEle.epOut() );
	   normHEEPEles_fbrem_.push_back(       ithHEEPEle.fbrem() );
	   normHEEPEles_bremFrac_.push_back(    ithHEEPEle.bremFrac() );
	   normHEEPEles_invEOverInvP_.push_back(ithHEEPEle.invEOverInvP() );

	   //Variables storing the heep::Ele method values - shower shape variables ...
	   normHEEPEles_sigmaEtaEta_.push_back(      ithHEEPEle.sigmaEtaEta() );
	   normHEEPEles_sigmaEtaEtaUnCorr_.push_back(ithHEEPEle.sigmaEtaEtaUnCorr() );
	   normHEEPEles_sigmaIEtaIEta_.push_back(    ithHEEPEle.sigmaIEtaIEta() );
	   normHEEPEles_e1x5_.push_back(             ithHEEPEle.e1x5() );
	   normHEEPEles_e2x5Max_.push_back(          ithHEEPEle.e2x5Max() );
	   normHEEPEles_e5x5_.push_back(             ithHEEPEle.e5x5() );
	   normHEEPEles_e1x5Over5x5_.push_back(      ithHEEPEle.e1x5Over5x5() );
	   normHEEPEles_e2x5MaxOver5x5_.push_back(   ithHEEPEle.e2x5MaxOver5x5() );

	   //Variables storing the heep::Ele method values - isolation variables ...
	   normHEEPEles_isolEm_.push_back(         ithHEEPEle.isolEm() );
	   normHEEPEles_isolHad_.push_back(        ithHEEPEle.isolHad() );
	   normHEEPEles_isolHadDepth1_.push_back(  ithHEEPEle.isolHadDepth1() );
	   normHEEPEles_isolHadDepth2_.push_back(  ithHEEPEle.isolHadDepth2() );
	   normHEEPEles_isolPtTrks_.push_back(     ithHEEPEle.isolPtTrks() );
	   normHEEPEles_isolEmHadDepth1_.push_back(ithHEEPEle.isolEmHadDepth1() );

	   normHEEPEles_numMissInnerHits_.push_back( ithHEEPEle.nrMissHits() );


	   ithHEEPEle.setEvtPrimVertexPos( firstPrimaryVertex.position() );
	   event_->AddStdEleInfo_dxy( ithHEEPEle.dxy() );

	   // For getting the caloGeometry information ... es.get<CaloGeometryRecord>().get(eventSetupData_->caloGeom);
	   // Variables storing the recHits information ...
	   std::vector<float> recHits_EtValues; recHits_EtValues.clear();
	   std::vector<float> recHits_etaValues; recHits_etaValues.clear();
	   std::vector<float> recHits_phiValues; recHits_phiValues.clear();
	   std::vector<bool>  recHits_isFromEBFlags; recHits_isFromEBFlags.clear();
//	   std::cout << "         -=-=-" << std::endl;
//	   std::cout << "         -=-=-" << std::endl;
//	   std::cout << "       SC INFORMATION ..." << std::endl;
	   reco::SuperClusterRef eleSC = ithGsfEle.superCluster();
	   GlobalPoint scPosition( eleSC->position().x(), eleSC->position().y(), eleSC->position().z());
	   float eleSC_eta = scPosition.eta();
	   float eleSC_phi = scPosition.phi();
	   float eleSC_rawEnergy = eleSC->rawEnergy();
	   float recHits_totEnergy = 0.0;
	   unsigned int recHits_totNum = 0;
	   for(reco::CaloCluster_iterator eleBCIt=eleSC->clustersBegin(); eleBCIt != eleSC->clustersEnd(); ++eleBCIt) {
//	   	std::cout << "            ->BasicCluster:" << std::endl;
	   	std::vector< std::pair<DetId,float> > eleBC_hitsAndFracs = (*eleBCIt)->hitsAndFractions();
	   	for(unsigned int recHitIdx = 0; recHitIdx<eleBC_hitsAndFracs.size(); recHitIdx++){
	   		double recHitEtMin = 0.0;
	   		double recHitEnergyMin = 0.0;
	   		bool recHit_isFromEB = true;
	   		DetId recHitDetId = eleBC_hitsAndFracs.at(recHitIdx).first;
	   		// Grab (an iterator for) this recHit for this collection; first look in the EB recHits collection
	   		CaloRecHitMetaCollectionV::const_iterator recHit_tmpIt = recHitsMeta_EB->find(recHitDetId);
   			recHitEtMin = 0.0; recHitEnergyMin = 0.08;
   			recHit_isFromEB = true;
	   		// ... then, if this recHit isn't in EB coll'n, look in the EE coll'n ...
	   		if(recHit_tmpIt==recHitsMeta_EB->end()){
	   			recHit_tmpIt = recHitsMeta_EE->find(recHitDetId);
	   			recHitEtMin = 0.1; recHitEnergyMin = 0.0; recHit_isFromEB = false;}
	   		// ... then if still haven't been able to find this recHit, skip it and move onto the next one! (This is exactly what is done in calculating GSF electron isolation variables)
	   		if(recHit_tmpIt==recHitsMeta_EE->end()){
	   			continue;
	   			std::cout << "  *** N.B: SC recHit not found in recHitsCollns ***" << std::endl;
	   		}

	   		// Now, add the (eta, phi) position, and the corresponding Et value for the recHit to the recHits vectors ...
	   		const GlobalPoint recHitPosn = handle_caloGeom.product()->getPosition(recHitDetId);
	   		float recHit_energy = recHit_tmpIt->energy();
	   		recHits_totEnergy += recHit_energy; recHits_totNum += 1;
	   		float recHit_Et = recHit_energy*recHitPosn.perp()/recHitPosn.mag();
//	   		std::cout << "                recHit #" << recHitIdx << ": (energy=" << recHit_energy << "; et=" << recHit_Et << "; eta=" << recHitPosn.eta() << "; phi=" << recHitPosn.phi() << ")" << std::endl;
	   		if(fabs(recHit_energy)>recHitEnergyMin && fabs(recHit_Et)>recHitEtMin){
	   			recHits_EtValues.push_back(  recHit_Et        );
	   			recHits_etaValues.push_back( recHitPosn.eta() );
	   			recHits_phiValues.push_back( recHitPosn.phi() );
	   			recHits_isFromEBFlags.push_back( recHit_isFromEB );
	   		}
	   	}// end of recHits for loop
	   }// end of BCs for loop

	   normHEEPEles_SCposn_eta_.push_back( eleSC_eta );
	   normHEEPEles_SCposn_phi_.push_back( eleSC_phi );
	   normHEEPEles_SC_rawEnergy_.push_back( eleSC_rawEnergy );
	   normHEEPEles_SC_recHits_Et_.push_back(  recHits_EtValues  );
	   normHEEPEles_SC_recHits_eta_.push_back( recHits_etaValues );
	   normHEEPEles_SC_recHits_phi_.push_back( recHits_phiValues );
	   normHEEPEles_SC_recHits_isFromEB_.push_back( recHits_isFromEBFlags );
	   normHEEPEles_SC_totEnergyRecHits_.push_back( recHits_totEnergy );
	   normHEEPEles_SC_totNumRecHits_.push_back( recHits_totNum );

	   // Variable storing eta, phi & vz values for GSF track associated with electron ...
	   const float ithEle_gsfTrk_eta = ithHEEPEle.gsfEle().gsfTrack()->eta();
	   normHEEPEles_gsfTrk_eta_.push_back( ithEle_gsfTrk_eta );
	   const float ithEle_gsfTrk_phi = ithHEEPEle.gsfEle().gsfTrack()->phi();
	   normHEEPEles_gsfTrk_phi_.push_back( ithEle_gsfTrk_phi );
	   normHEEPEles_gsfTrk_vz_.push_back( ithHEEPEle.gsfEle().gsfTrack()->vz() );

	   // Getting values of variable storing pt, eta, phi and vz values for each CTF track with pT>0.7GeV,
	   // within the electron's inner isolation cone.
	   std::vector<float> innerIsoConeTrks_pt;
	   std::vector<float> innerIsoConeTrks_eta;
	   std::vector<float> innerIsoConeTrks_phi;
	   std::vector<float> innerIsoConeTrks_vz;

	   for(reco::TrackCollection::const_iterator itrTrk = handle_ctfTracks.product()->begin();
	   		itrTrk != handle_ctfTracks.product()->end(); itrTrk++){
	   	// Get pT of track, and skip track if pT<0.07GeV
	   	float ithTrk_pt = (*itrTrk).pt();
	   	if(ithTrk_pt<0.7)
	   		continue;

	   	// Get eta, phi and vz values for track
	   	float ithTrk_eta = (*itrTrk).eta();
	   	float ithTrk_phi = (*itrTrk).phi();
	   	float ithTrk_vz  = (*itrTrk).vz();

	   	// Skip track if it is not in this electron's inner track isolation 'area' - i.e. if it is outside of the ...
	   	// ... inner strip and outside of the inner cone, OR if it is NOT within the outer track isol cone
	   	float dEta = ithTrk_eta - ithEle_gsfTrk_eta;
	   	float dPhi = ithTrk_phi - ithEle_gsfTrk_phi;
			while( dPhi<=-1.0*TMath::Pi() ){ dPhi = dPhi+2.0*TMath::Pi(); }
			while( dPhi>TMath::Pi()      ){ dPhi = dPhi-2.0*TMath::Pi(); }
			float dR_ithTrk = sqrt(dEta*dEta+dPhi*dPhi);
	   	if(fabs(dEta)>0.015 && dR_ithTrk>0.015)
	   		continue;
	   	if(dR_ithTrk>0.3)
	   		continue;

	   	// If have got this far, then the track has high enough pT to be included in isolation variable, and ...
	   	// ... it is within this electron's inner isolation 'area' => store its information in branch
	   	innerIsoConeTrks_pt.push_back(  ithTrk_pt  );
	   	innerIsoConeTrks_eta.push_back( ithTrk_eta );
	   	innerIsoConeTrks_phi.push_back( ithTrk_phi );
	   	innerIsoConeTrks_vz.push_back(  ithTrk_vz  );
	   }

	   normHEEPEles_innerIsoConeTrks_pt_.push_back(  innerIsoConeTrks_pt  );
	   normHEEPEles_innerIsoConeTrks_eta_.push_back( innerIsoConeTrks_eta );
	   normHEEPEles_innerIsoConeTrks_phi_.push_back( innerIsoConeTrks_phi );
	   normHEEPEles_innerIsoConeTrks_vz_.push_back(  innerIsoConeTrks_vz  );

		// Reading-in ele isolation values re-calculated with IsoDeposits ...
		reco::GsfElectronRef ithGsfEleRef(handle_normGsfEles, eleIdx);
		event_->AddStdEleInfo_isoDep_std( (*h_stdEleIsoFromDeps_Tk)[ithGsfEleRef], (*h_stdEleIsoFromDeps_Ecal)[ithGsfEleRef], (*h_stdEleIsoFromDeps_HcalD1)[ithGsfEleRef] );
		event_->AddStdEleInfo_isoDep_inrVeto( (*h_modEleIsoFromDeps_Tk)[ithGsfEleRef], (*h_modEleIsoFromDeps_Ecal)[ithGsfEleRef], (*h_modEleIsoFromDeps_HcalD1)[ithGsfEleRef] );
		event_->AddStdEleInfo_inrVetoModIso(tsw::Event::xSmallVeto, (*h_inrXSVetoEleIso_Tk)[ithGsfEleRef], (*h_inrXSVetoEleIso_Ecal)[ithGsfEleRef], (*h_inrXSVetoEleIso_HcalD1)[ithGsfEleRef] );
		event_->AddStdEleInfo_inrVetoModIso(tsw::Event::smallVeto,  (*h_inrSVetoEleIso_Tk )[ithGsfEleRef], (*h_inrSVetoEleIso_Ecal )[ithGsfEleRef], (*h_inrSVetoEleIso_HcalD1 )[ithGsfEleRef] );
		event_->AddStdEleInfo_inrVetoModIso(tsw::Event::mediumVeto, (*h_inrMVetoEleIso_Tk )[ithGsfEleRef], (*h_inrMVetoEleIso_Ecal )[ithGsfEleRef], (*h_inrMVetoEleIso_HcalD1 )[ithGsfEleRef] );
		event_->AddStdEleInfo_inrVetoModIso(tsw::Event::largeVeto,  (*h_inrLVetoEleIso_Tk )[ithGsfEleRef], (*h_inrLVetoEleIso_Ecal )[ithGsfEleRef], (*h_inrLVetoEleIso_HcalD1 )[ithGsfEleRef] );
		event_->AddStdEleInfo_inrVetoModIso(tsw::Event::xLargeVeto, (*h_inrXLVetoEleIso_Tk)[ithGsfEleRef], (*h_inrXLVetoEleIso_Ecal)[ithGsfEleRef], (*h_inrXLVetoEleIso_HcalD1)[ithGsfEleRef] );


		typedef std::vector<edm::Handle< edm::ValueMap<double> > >::const_iterator  DblValMapHandleVecIt ;
		DblValMapHandleVecIt hIt_inrVetoModIsosPhantomEle_dEta = hVec_inrVetoModIsosPhantomEle_dEta.begin();
		DblValMapHandleVecIt hIt_inrVetoModIsosPhantomEle_dPhi = hVec_inrVetoModIsosPhantomEle_dPhi.begin();
		DblValMapHandleVecIt hIt_inrVetoModTkIsosPhantomEle    = hVec_inrVetoModTkIsosPhantomEle.begin();
		DblValMapHandleVecIt hIt_inrVetoModEmIsosPhantomEle    = hVec_inrVetoModEmIsosPhantomEle.begin();
		DblValMapHandleVecIt hIt_inrVetoModH1IsosPhantomEle    = hVec_inrVetoModH1IsosPhantomEle.begin();
		std::vector<double> vec_inrVetoModIsosPhantomEle_dEta, vec_inrVetoModIsosPhantomEle_dPhi,
				vec_inrVetoModTkIsosPhantomEle, vec_inrVetoModEmIsosPhantomEle, vec_inrVetoModH1IsosPhantomEle;
		for( ; hIt_inrVetoModIsosPhantomEle_dEta != hVec_inrVetoModIsosPhantomEle_dEta.end(); ){
			vec_inrVetoModIsosPhantomEle_dEta.push_back( (**hIt_inrVetoModIsosPhantomEle_dEta)[ithGsfEleRef] );
			hIt_inrVetoModIsosPhantomEle_dEta++;
			vec_inrVetoModIsosPhantomEle_dPhi.push_back( (**hIt_inrVetoModIsosPhantomEle_dPhi)[ithGsfEleRef] );
			hIt_inrVetoModIsosPhantomEle_dPhi++;

			vec_inrVetoModTkIsosPhantomEle.push_back( (**hIt_inrVetoModTkIsosPhantomEle)[ithGsfEleRef] );
			hIt_inrVetoModTkIsosPhantomEle++;
			vec_inrVetoModEmIsosPhantomEle.push_back( (**hIt_inrVetoModEmIsosPhantomEle)[ithGsfEleRef] );
			hIt_inrVetoModEmIsosPhantomEle++;
			vec_inrVetoModH1IsosPhantomEle.push_back( (**hIt_inrVetoModH1IsosPhantomEle)[ithGsfEleRef] );
			hIt_inrVetoModH1IsosPhantomEle++;
		}
		event_->AddStdEleInfo_inrVetoModIsoWithPhantomEle(vec_inrVetoModIsosPhantomEle_dEta, vec_inrVetoModIsosPhantomEle_dPhi,
				vec_inrVetoModTkIsosPhantomEle, vec_inrVetoModEmIsosPhantomEle, vec_inrVetoModH1IsosPhantomEle );

		// Calculating the number of hadrons within a dR=0.4 cone around each electron (for separation of underlying event & other-electron-footprint effects from isolation cut efficiency curves)
		unsigned int nGenHadrons_dR04 = 0;
		double   ptSumGenHadrons_dR04 = 0.0;

		edm::Handle<reco::GenParticleCollection> h_genParticles;
		iEvent.getByLabel("genParticles", h_genParticles);

		if( !(h_genParticles.failedToGet()) ){
			const TLorentzVector gsfP4 = TLorentzVector(ithGsfEle.px(), ithGsfEle.py(), ithGsfEle.pz(), ithGsfEle.energy());
			for(size_t iGen=0; iGen<h_genParticles->size(); iGen++){
				const reco::GenParticle& genParticle = h_genParticles->at(iGen);

				if(genParticle.pt()<0.4) // Check that particle's pT is > 0.4GeV
					continue;

				if(genParticle.status()==3 || genParticle.status()==2) // Check that particle is not part of hard interaction
					continue;
				int pdgId = genParticle.pdgId();
				if(abs(pdgId)>=11 && abs(pdgId)<=18) // Check that particle is not a lepton
					continue;
				else if(abs(pdgId)==4000001 || abs(pdgId)==4000002 || abs(pdgId)==4000011 || abs(pdgId)==4000012 ) // Check that particle is not an excited fermion
					continue;
				else if(abs(pdgId)==9 || (abs(pdgId)>=21 && abs(pdgId)<=25) || (abs(pdgId)>=32 && abs(pdgId)<=37) ) // Check that particle is not a gauge/Higgs boson
					continue;
				// If have got this far then particle must be a hadron from the underlying event (at least for 2011 signal/DY+Jets MC)

				if( fabs(ithGsfEle.eta()-genParticle.eta())>0.4 )
					continue;
				const TLorentzVector genP4 = TLorentzVector(genParticle.px(),  genParticle.py(), genParticle.pz(), genParticle.energy());
				const double drVsGsfEle = genP4.DeltaR(gsfP4);
				if(drVsGsfEle<0.4){
					nGenHadrons_dR04++;
					ptSumGenHadrons_dR04 += genParticle.pt();
				}
			}
		} //End: if got genParticles from event
		event_->AddStdEleInfo_genHadronsDr04(nGenHadrons_dR04, ptSumGenHadrons_dR04);

	} //End: for loop over GSF eles
	
//	//Printing these values to screen...
//	for(unsigned int iEle = 0; iEle < normGsfEles_p4ptr_->size(); iEle++){
//		std::cout << "  * GSF ele no. " << iEle << ":" << std::endl
//				    << "       Et=" << normGsfEles_p4_.at(iEle).Pt() << " ; eta=" << normGsfEles_p4_.at(iEle).Eta() << " ; phi=" << normGsfEles_p4_.at(iEle).Phi() << std::endl;
//		std::cout << "       Isos: track=" << normHEEPEles_isolPtTrks_.at(iEle) << " ; ecal=" << normHEEPEles_isolEm_.at(iEle) << " ; hcalD1=" << normHEEPEles_isolHadDepth1_.at(iEle) << std::endl;
//		event_->PrintStdEleInfo_inrVetoModIsoWithPhantomEle(iEle);
//		std::cout << std::endl;
//	}

	if(beVerbose){
		for(unsigned int iEle = 0; iEle < normGsfEles_p4ptr_->size(); iEle++){
			std::cout << "     norm GSF ele no. " << iEle << ":";
			std::cout << " Charge = " << normGsfEles_charge_.at(iEle) << std::endl;

			std::cout << "       p4: Px=" << normGsfEles_p4_.at(iEle).Px() << "; Py=" << normGsfEles_p4_.at(iEle).Py() << "; Pz=" << normGsfEles_p4_.at(iEle).Pz() << std::endl;
			std::cout << "           Px=" << normGsfEles_p4ptr_->at(iEle).Px() << "; Py=" << normGsfEles_p4ptr_->at(iEle).Py() << "; Pz=" << normGsfEles_p4ptr_->at(iEle).Pz() << std::endl;
			std::cout << "           Pt=" << normGsfEles_p4ptr_->at(iEle).Pt() << "; Eta=" << normGsfEles_p4ptr_->at(iEle).Eta() << "; Phi=" << normGsfEles_p4ptr_->at(iEle).Phi() << std::endl;

			std::cout << "       Et=" << normGsfEles_Et_.at(iEle);
			std::cout << "; HEEP_Et=" << normGsfEles_HEEP_Et_.at(iEle) << std::endl;
			std::cout << "       Eta=" << normGsfEles_Eta_.at(iEle);
			std::cout << "; scEta=" << normGsfEles_scEta_.at(iEle) << std::endl;
			std::cout << "       ecalDriven=" << normGsfEles_ecalDriven_.at(iEle);
			std::cout << "; ecalDrivenSeed=" << normGsfEles_ecalDrivenSeed_.at(iEle) << std::endl;

			std::cout << "       HEEP_dEtaIn=" << normGsfEles_HEEP_dEtaIn_.at(iEle);
			std::cout << "; HEEP_dPhiIn=" << normGsfEles_HEEP_dPhiIn_.at(iEle) << std::endl;
			std::cout << "       HoverE=" << normGsfEles_HoverE_.at(iEle);
			std::cout << "; sigmaIetaIeta=" << normGsfEles_sigmaIetaIeta_.at(iEle);
			std::cout << "; scSigmaIetaIeta=" << normGsfEles_scSigmaIetaIeta_.at(iEle) << std::endl;
	
			std::cout << "       dr03EmIsoEt=" << normGsfEles_dr03EmIsoEt_.at(iEle);	
			std::cout << "; dr03HadDepth1IsoEt=" << normGsfEles_dr03HadDepth1IsoEt_.at(iEle) << std::endl;
			std::cout << "       dr03HadDepth2IsoEt=" << normGsfEles_dr03HadDepth2IsoEt_.at(iEle);
			std::cout << "; dr03TkIsoPt=" << normGsfEles_dr03TkIsoPt_.at(iEle) << std::endl;

			std::cout << "       e2x5Max=" << normGsfEles_e2x5Max_.at(iEle);
			std::cout << "; e5x5=" << normGsfEles_e5x5_.at(iEle) << std::endl;
			
			std::cout << "    ->HEEP method values ..." << std::endl;
			//Variables storing the heep::Ele method values ...
			std::cout << "       et=" << normHEEPEles_et_.at(iEle);
			std::cout << "; gsfEt=" << normHEEPEles_gsfEt_.at(iEle);
			std::cout << "; scEt=" << normHEEPEles_scEt_.at(iEle) << std::endl;
			std::cout << "       energy=" << normHEEPEles_energy_.at(iEle);
			std::cout << "; gsfEnergy=" << normHEEPEles_gsfEnergy_.at(iEle);
			std::cout << "; caloEnergy=" << normHEEPEles_caloEnergy_.at(iEle);
			std::cout << "; ecalEnergyError=" << normHEEPEles_ecalEnergyError_.at(iEle) << std::endl;
			std::cout << "       eta=" << normHEEPEles_eta_.at(iEle);
			std::cout << "; scEta=" << normHEEPEles_scEta_.at(iEle);
			std::cout << "; detEta=" << normHEEPEles_detEta_.at(iEle);
			std::cout << "; detEtaAbs=" << normHEEPEles_detEtaAbs_.at(iEle) << std::endl;
			std::cout << "       phi=" << normHEEPEles_phi_.at(iEle);
			std::cout << "; scPhi=" << normHEEPEles_scPhi_.at(iEle);
			std::cout << "; detPhi=" << normHEEPEles_detPhi_.at(iEle);
			std::cout << "; zVtx=" << normHEEPEles_zVtx_.at(iEle) << std::endl;
			std::cout << "       p4,E=" << normHEEPEles_p4_.at(iEle).E() << "; p4,p=(" << normHEEPEles_p4_.at(iEle).Px() << ", " << normHEEPEles_p4_.at(iEle).Py() << ", " << normHEEPEles_p4_.at(iEle).Pz() << ")" << std::endl;
			std::cout << "       gsfP4,E=" << normHEEPEles_gsfP4_.at(iEle).E() << "; gsfP4,p=(" << normHEEPEles_gsfP4_.at(iEle).Px() << ", " << normHEEPEles_gsfP4_.at(iEle).Py() << ", " << normHEEPEles_gsfP4_.at(iEle).Pz() << ")" << std::endl;

			//Variables storing the heep::Ele method values - Classification...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       classification=" << normHEEPEles_classification_.at(iEle);
			std::cout << "; isEcalDriven=" << normHEEPEles_isEcalDriven_.at(iEle);
			std::cout << "; isTrakerDriven=" << normHEEPEles_isTrackerDriven_.at(iEle) << std::endl;
			std::cout << "       isEB=" << normHEEPEles_isEB_.at(iEle);
			std::cout << "; isEE=" << normHEEPEles_isEE_.at(iEle) <<std::endl;

			//Variables storing the heep::Ele method values - track methods ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       charge=" << normHEEPEles_charge_.at(iEle);
			std::cout << "; trkCharge=" << normHEEPEles_trkCharge_.at(iEle) << std::endl;
			std::cout << "       pVtx=" << normHEEPEles_pVtx_.at(iEle);
			std::cout << "; pCalo=" << normHEEPEles_pCalo_.at(iEle);
			std::cout << "; ptVtx=" << normHEEPEles_ptVtx_.at(iEle);
			std::cout << "; ptCalo=" << normHEEPEles_ptCalo_.at(iEle) << std::endl;
			std::cout << "       Closest CTF Track... (pT, eta, phi)=(" << normHEEPEles_closestCtfTrk_pt_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_eta_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_phi_.at(iEle) << ")" << std::endl;
			std::cout << "         Start of track (inner): (pT, eta, phi)=(" << normHEEPEles_closestCtfTrk_innerPt_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_innerEta_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_innerPhi_.at(iEle) << ")" << std::endl;
			std::cout << "         End of track (outer):   (pT, eta, phi)=(" << normHEEPEles_closestCtfTrk_outerPt_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_outerEta_.at(iEle) << ", " << normHEEPEles_closestCtfTrk_outerPhi_.at(iEle) << ")" << std::endl;

			//Variables storing the heep::Ele method values ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       hOver=" << normHEEPEles_hOverE_.at(iEle);
			std::cout << "; dEtaIn=" << normHEEPEles_dEtaIn_.at(iEle);
			std::cout << "; dPhiIn=" << normHEEPEles_dPhiIn_.at(iEle);
			std::cout << "; dPhiOut=" << normHEEPEles_dPhiOut_.at(iEle) << std::endl;
			std::cout << "       epIn=" << normHEEPEles_epIn_.at(iEle);
			std::cout << "; epOut=" << normHEEPEles_epOut_.at(iEle);
			std::cout << "; fbrem=" << normHEEPEles_fbrem_.at(iEle);
			std::cout << "; bremFrac=" << normHEEPEles_bremFrac_.at(iEle);
			std::cout << "; invEOverInvP=" << normHEEPEles_invEOverInvP_.at(iEle) << std::endl;

		   //Variables storing the heep::Ele method values - shower shape variables ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       sigmaEtaEta=" << normHEEPEles_sigmaEtaEta_.at(iEle);
		   std::cout << "; sigmaEtaEtaUnCorr=" << normHEEPEles_sigmaEtaEtaUnCorr_.at(iEle);
		   std::cout << "; sigmaIEtaIEta=" << normHEEPEles_sigmaIEtaIEta_.at(iEle) << std::endl;
		   std::cout << "       e1x5=" << normHEEPEles_e1x5_.at(iEle);
		   std::cout << "; e2x5Max=" << normHEEPEles_e2x5Max_.at(iEle);
		   std::cout << "; e5x5=" << normHEEPEles_e5x5_.at(iEle);
		   std::cout << "; e1x5Over5x5=" << normHEEPEles_e1x5Over5x5_.at(iEle);
		   std::cout << "; e2x5MaxOver5x5=" << normHEEPEles_e2x5MaxOver5x5_.at(iEle) << std::endl;

		   //Variables storing the heep::Ele method values - isolation variables ...
		   std::cout << "         -=-=-" << std::endl;
		   std::cout << "       isolEm=" << normHEEPEles_isolEm_.at(iEle);
		   std::cout << "; isolHad=" << normHEEPEles_isolHad_.at(iEle);
		   std::cout << "; isolHadDepth1=" << normHEEPEles_isolHadDepth1_.at(iEle);
		   std::cout << "; isolHadDepth2=" << normHEEPEles_isolHadDepth2_.at(iEle) << std::endl;
		   std::cout << "       isolPtTrks=" << normHEEPEles_isolPtTrks_.at(iEle);
		   std::cout << "; isolEmHadDepth1=" << normHEEPEles_isolEmHadDepth1_.at(iEle) << std::endl;

		   std::cout << "       numMissInnerHits = " << normHEEPEles_numMissInnerHits_.at(iEle) << std::endl;

		   event_->PrintStdEleInfo_isoValues(iEle);

		   std::cout << "         -=-=-" << std::endl;
		   std::cout << "       recHits... (tot num =" << normHEEPEles_SC_totNumRecHits_.at(iEle) << ", totEnergy=" << normHEEPEles_SC_totEnergyRecHits_.at(iEle) << "; num stored=" << normHEEPEles_SC_recHits_Et_.at(iEle).size() << "=" << normHEEPEles_SC_recHits_eta_.at(iEle).size() << "=" << normHEEPEles_SC_recHits_phi_.at(iEle).size() << "=" << normHEEPEles_SC_recHits_isFromEB_.at(iEle).size() << "?? of them)" << std::endl;
		   std::cout << "           GlobalPoint pos'n: (eta,phi) = (" << normHEEPEles_SCposn_eta_.at(iEle) << "," << normHEEPEles_SCposn_phi_.at(iEle) << "); rawEnergy=" << normHEEPEles_SC_rawEnergy_.at(iEle) << std::endl;
		   for(unsigned int recHitIdx=0; recHitIdx<normHEEPEles_SC_recHits_Et_.at(iEle).size(); recHitIdx++){
		   	std::cout << "          #" << recHitIdx << ": (Et, eta, phi; isFromEB)=(" << normHEEPEles_SC_recHits_Et_.at(iEle).at(recHitIdx) << ", " << normHEEPEles_SC_recHits_eta_.at(iEle).at(recHitIdx) << ", " << normHEEPEles_SC_recHits_phi_.at(iEle).at(recHitIdx) << "; " << normHEEPEles_SC_recHits_isFromEB_.at(iEle).at(recHitIdx) << ")" << std::endl;
		   }

		   // Variables storing eta, phi & vz values for GSF track associated with electron ...
		   std::cout << "         -=-=-" << std::endl;
		   std::cout << "       gsfTrack... (eta, phi) = (" << normHEEPEles_gsfTrk_eta_.at(iEle) << ", " << normHEEPEles_gsfTrk_phi_.at(iEle) << "), vz=" <<	normHEEPEles_gsfTrk_vz_.at( iEle ) << std::endl;

			// Getting values of variable storing pt, eta, phi and vz values for each CTF track with pT>0.7GeV,
			// within the electron's inner isolation cone.
		   std::cout << "         -=-=-" << std::endl;
		   const std::vector<float> innerIsoConeTrks_pt = normHEEPEles_innerIsoConeTrks_pt_.at(iEle);
		   const std::vector<float> innerIsoConeTrks_eta = normHEEPEles_innerIsoConeTrks_eta_.at(iEle);
		   const std::vector<float> innerIsoConeTrks_phi = normHEEPEles_innerIsoConeTrks_phi_.at(iEle);
		   const std::vector<float> innerIsoConeTrks_vz = normHEEPEles_innerIsoConeTrks_vz_.at(iEle);
		   std::cout << "       CTF tracks w/i inner isol area (" << innerIsoConeTrks_pt.size() << "=" << innerIsoConeTrks_eta.size() << "=" << innerIsoConeTrks_phi.size() << "=" << innerIsoConeTrks_vz.size() << " of them)... " << std::endl;
		   for(unsigned int trkIdx=0; trkIdx<innerIsoConeTrks_pt.size(); trkIdx++)
		   	std::cout << "          #" << trkIdx << ": (pT, eta, phi) = (" << innerIsoConeTrks_pt.at(trkIdx) << ", " << innerIsoConeTrks_eta.at(trkIdx) << ", " << innerIsoConeTrks_phi.at(trkIdx) << "); vz=" << innerIsoConeTrks_vz.at(trkIdx) << std::endl;
		}
	}
	
}


//------------ method for reading in the values of the standard GSF electron variables ---------------
void BstdZeeNTupler::ReadInBstdGsfEles(bool beVerbose, const edm::Handle<reco::GsfElectronCollection>& handle_bstdGsfEles){
	
	reco::GsfElectron ithGsfEle;
	heep::Ele ithHEEPEle(ithGsfEle);

	//Setting the values of the standard reconstruction GSF electron variables...
	bstdGsfEles_number_ = handle_bstdGsfEles.product()->size();
	if(beVerbose){std::cout << " ->There are " << bstdGsfEles_number_ << " GSF from the special boosted Z(ee) reco'n in this event."<< std::endl;}
	
	for(unsigned int eleIdx = 0; eleIdx < bstdGsfEles_number_; eleIdx++){
		
		ithGsfEle = handle_bstdGsfEles.product()->at(eleIdx);
		ithHEEPEle = heep::Ele(ithGsfEle);

		bstdGsfEles_p4ptr_->push_back(            ithGsfEle.p4()     );
		bstdGsfEles_charge_.push_back(            ithGsfEle.charge() );
 	
		bstdGsfEles_Et_.push_back(                ithGsfEle.et() );
		bstdGsfEles_HEEP_Et_.push_back(           ithGsfEle.caloEnergy()*sin(ithGsfEle.p4().theta()) );
		bstdGsfEles_Eta_.push_back(               ithGsfEle.eta()              );
		bstdGsfEles_scEta_.push_back(             ithGsfEle.caloPosition().eta() );

		bstdGsfEles_ecalDriven_.push_back(        ithGsfEle.ecalDriven()       );
		bstdGsfEles_ecalDrivenSeed_.push_back(    ithGsfEle.ecalDrivenSeed()   );

		bstdGsfEles_HEEP_dEtaIn_.push_back(       ithGsfEle.deltaEtaSuperClusterTrackAtVtx() );
		bstdGsfEles_HEEP_dPhiIn_.push_back(       ithGsfEle.deltaPhiSuperClusterTrackAtVtx() );
		bstdGsfEles_HoverE_.push_back(            ithGsfEle.hadronicOverEm()  );
		bstdGsfEles_sigmaIetaIeta_.push_back(     ithGsfEle.sigmaIetaIeta()   );
		bstdGsfEles_scSigmaIetaIeta_.push_back(   ithGsfEle.scSigmaIEtaIEta() );

		bstdGsfEles_dr03EmIsoEt_.push_back(        ithGsfEle.dr03EcalRecHitSumEt()      );
		bstdGsfEles_dr03HadDepth1IsoEt_.push_back( ithGsfEle.dr03HcalDepth1TowerSumEt() );
		bstdGsfEles_dr03HadDepth2IsoEt_.push_back( ithGsfEle.dr03HcalDepth2TowerSumEt() );
		bstdGsfEles_dr03TkIsoPt_.push_back(        ithGsfEle.dr03TkSumPt()              );

		bstdGsfEles_e2x5Max_.push_back(            ithGsfEle.e2x5Max() );
		bstdGsfEles_e5x5_.push_back(               ithGsfEle.e5x5()    );

		//Variables storing the heep::Ele method values ...
		bstdHEEPEles_et_.push_back(        ithHEEPEle.et() );
		bstdHEEPEles_gsfEt_.push_back(     ithHEEPEle.gsfEt() );
		bstdHEEPEles_scEt_.push_back(      ithHEEPEle.scEt() );
		bstdHEEPEles_energy_.push_back(    ithHEEPEle.energy() );
		bstdHEEPEles_gsfEnergy_.push_back( ithHEEPEle.gsfEnergy() );
		bstdHEEPEles_caloEnergy_.push_back(ithHEEPEle.caloEnergy() );
		bstdHEEPEles_eta_.push_back(       ithHEEPEle.eta() );
		bstdHEEPEles_scEta_.push_back(     ithHEEPEle.scEta() );
		bstdHEEPEles_detEta_.push_back(    ithHEEPEle.detEta() );
		bstdHEEPEles_detEtaAbs_.push_back( ithHEEPEle.detEtaAbs() );
		bstdHEEPEles_phi_.push_back(       ithHEEPEle.phi() );
		bstdHEEPEles_scPhi_.push_back(     ithHEEPEle.scPhi() );
		bstdHEEPEles_detPhi_.push_back(    ithHEEPEle.detPhi() );
		bstdHEEPEles_zVtx_.push_back(      ithHEEPEle.zVtx() );
		bstdHEEPEles_p4_.push_back(        ithHEEPEle.p4() );
		bstdHEEPEles_gsfP4_.push_back(     ithHEEPEle.gsfP4() );

		//Variables storing the heep::Ele method values - Classification...
		bstdHEEPEles_classification_.push_back( ithHEEPEle.classification() );
		bstdHEEPEles_isEcalDriven_.push_back(   ithHEEPEle.isEcalDriven() );
		bstdHEEPEles_isTrackerDriven_.push_back(ithHEEPEle.isTrackerDriven() );
		bstdHEEPEles_isEB_.push_back(           ithHEEPEle.isEB() );
		bstdHEEPEles_isEE_.push_back(           ithHEEPEle.isEE() );

		//Variables storing the heep::Ele method values - track methods ...
		bstdHEEPEles_charge_.push_back(   ithHEEPEle.charge() );
		bstdHEEPEles_trkCharge_.push_back(ithHEEPEle.trkCharge() );
		bstdHEEPEles_pVtx_.push_back(     ithHEEPEle.pVtx() );
		bstdHEEPEles_pCalo_.push_back(    ithHEEPEle.pCalo() );
		bstdHEEPEles_ptVtx_.push_back(    ithHEEPEle.ptVtx() );
		bstdHEEPEles_ptCalo_.push_back(   ithHEEPEle.ptCalo() );

		//Variables storing the heep::Ele method values ...
		bstdHEEPEles_hOverE_.push_back(      ithHEEPEle.hOverE() );
		bstdHEEPEles_dEtaIn_.push_back(      ithHEEPEle.dEtaIn() );
		bstdHEEPEles_dPhiIn_.push_back(      ithHEEPEle.dPhiIn() );
		bstdHEEPEles_dPhiOut_.push_back(     ithHEEPEle.dPhiOut() );
		bstdHEEPEles_epIn_.push_back(        ithHEEPEle.epIn() );
		bstdHEEPEles_epOut_.push_back(       ithHEEPEle.epOut() );
		bstdHEEPEles_fbrem_.push_back(       ithHEEPEle.fbrem() );
		bstdHEEPEles_bremFrac_.push_back(    ithHEEPEle.bremFrac() );
		bstdHEEPEles_invEOverInvP_.push_back(ithHEEPEle.invEOverInvP() );

		//Variables storing the heep::Ele method values - shower shape variables ...
		bstdHEEPEles_sigmaEtaEta_.push_back(      ithHEEPEle.sigmaEtaEta() );
		bstdHEEPEles_sigmaEtaEtaUnCorr_.push_back(ithHEEPEle.sigmaEtaEtaUnCorr() );
		bstdHEEPEles_sigmaIEtaIEta_.push_back(    ithHEEPEle.sigmaIEtaIEta() );
		bstdHEEPEles_e1x5_.push_back(             ithHEEPEle.e1x5() );
		bstdHEEPEles_e2x5Max_.push_back(          ithHEEPEle.e2x5Max() );
		bstdHEEPEles_e5x5_.push_back(             ithHEEPEle.e5x5() );
		bstdHEEPEles_e1x5Over5x5_.push_back(      ithHEEPEle.e1x5Over5x5() );
		bstdHEEPEles_e2x5MaxOver5x5_.push_back(   ithHEEPEle.e2x5MaxOver5x5() );

		//Variables storing the heep::Ele method values - isolation variables ...
		bstdHEEPEles_isolEm_.push_back(         ithHEEPEle.isolEm() );
		bstdHEEPEles_isolHad_.push_back(        ithHEEPEle.isolHad() );
		bstdHEEPEles_isolHadDepth1_.push_back(  ithHEEPEle.isolHadDepth1() );
		bstdHEEPEles_isolHadDepth2_.push_back(  ithHEEPEle.isolHadDepth2() );
		bstdHEEPEles_isolPtTrks_.push_back(     ithHEEPEle.isolPtTrks() );
		bstdHEEPEles_isolEmHadDepth1_.push_back(ithHEEPEle.isolEmHadDepth1() );

	}
	
	//Printing these values to screen...
	if(beVerbose){
		for(unsigned int iEle = 0; iEle < bstdGsfEles_p4ptr_->size(); iEle++){
			std::cout << "     bstd GSF ele no. " << iEle << ":";
			std::cout << " Charge = " << bstdGsfEles_charge_.at(iEle) << std::endl;

			std::cout << "       p4: Px=" << bstdGsfEles_p4_.at(iEle).Px() << "; Py=" << bstdGsfEles_p4_.at(iEle).Py() << "; Pz=" << bstdGsfEles_p4_.at(iEle).Pz() << std::endl;
			std::cout << "           Px=" << bstdGsfEles_p4ptr_->at(iEle).Px() << "; Py=" << bstdGsfEles_p4ptr_->at(iEle).Py() << "; Pz=" << bstdGsfEles_p4ptr_->at(iEle).Pz() << std::endl;
			std::cout << "           Pt=" << bstdGsfEles_p4ptr_->at(iEle).Pt() << "; Eta=" << bstdGsfEles_p4ptr_->at(iEle).Eta() << "; Phi=" << bstdGsfEles_p4ptr_->at(iEle).Phi() << std::endl;

			std::cout << "       Et=" << bstdGsfEles_Et_.at(iEle);
			std::cout << "; HEEP_Et=" << bstdGsfEles_HEEP_Et_.at(iEle) << std::endl;
			std::cout << "       Eta=" << bstdGsfEles_Eta_.at(iEle);
			std::cout << "; scEta=" << bstdGsfEles_scEta_.at(iEle) << std::endl;
			std::cout << "       ecalDriven=" << bstdGsfEles_ecalDriven_.at(iEle);
			std::cout << "; ecalDrivenSeed=" << bstdGsfEles_ecalDrivenSeed_.at(iEle) << std::endl;

			std::cout << "       HEEP_dEtaIn=" << bstdGsfEles_HEEP_dEtaIn_.at(iEle);
			std::cout << "; HEEP_dPhiIn=" << bstdGsfEles_HEEP_dPhiIn_.at(iEle) << std::endl;
			std::cout << "       HoverE=" << bstdGsfEles_HoverE_.at(iEle);
			std::cout << "; sigmaIetaIeta=" << bstdGsfEles_sigmaIetaIeta_.at(iEle);
			std::cout << "; scSigmaIetaIeta=" << bstdGsfEles_scSigmaIetaIeta_.at(iEle) << std::endl;
	
			std::cout << "       dr03EmIsoEt=" << bstdGsfEles_dr03EmIsoEt_.at(iEle);
			std::cout << "; dr03HadDepth1IsoEt=" << bstdGsfEles_dr03HadDepth1IsoEt_.at(iEle) << std::endl;
			std::cout << "       dr03HadDepth2IsoEt=" << bstdGsfEles_dr03HadDepth2IsoEt_.at(iEle);
			std::cout << "; dr03TkIsoPt=" << bstdGsfEles_dr03TkIsoPt_.at(iEle) << std::endl;

			std::cout << "       e2x5Max=" << bstdGsfEles_e2x5Max_.at(iEle);
			std::cout << "; e5x5=" << bstdGsfEles_e5x5_.at(iEle) << std::endl;
			
			std::cout << "    ->HEEP method values ..." << std::endl;
			//Variables storing the heep::Ele method values ...
			std::cout << "       et=" << bstdHEEPEles_et_.at(iEle);
			std::cout << "; gsfEt=" << bstdHEEPEles_gsfEt_.at(iEle);
			std::cout << "; scEt=" << bstdHEEPEles_scEt_.at(iEle) << std::endl;
			std::cout << "       energy=" << bstdHEEPEles_energy_.at(iEle);
			std::cout << "; gsfEnergy=" << bstdHEEPEles_gsfEnergy_.at(iEle);
			std::cout << "; caloEnergy=" << bstdHEEPEles_caloEnergy_.at(iEle) << std::endl;
			std::cout << "       eta=" << bstdHEEPEles_eta_.at(iEle);
			std::cout << "; scEta=" << bstdHEEPEles_scEta_.at(iEle);
			std::cout << "; detEta=" << bstdHEEPEles_detEta_.at(iEle);
			std::cout << "; detEtaAbs=" << bstdHEEPEles_detEtaAbs_.at(iEle) << std::endl;
			std::cout << "       phi=" << bstdHEEPEles_phi_.at(iEle);
			std::cout << "; scPhi=" << bstdHEEPEles_scPhi_.at(iEle);
			std::cout << "; detPhi=" << bstdHEEPEles_detPhi_.at(iEle);
			std::cout << "; zVtx=" << bstdHEEPEles_zVtx_.at(iEle) << std::endl;
			std::cout << "       p4,E=" << bstdHEEPEles_p4_.at(iEle).E() << "; p4,p=(" << bstdHEEPEles_p4_.at(iEle).Px() << ", " << bstdHEEPEles_p4_.at(iEle).Py() << ", " << bstdHEEPEles_p4_.at(iEle).Pz() << ")" << std::endl;
			std::cout << "       gsfP4,E=" << bstdHEEPEles_gsfP4_.at(iEle).E() << "; gsfP4,p=(" << bstdHEEPEles_gsfP4_.at(iEle).Px() << ", " << bstdHEEPEles_gsfP4_.at(iEle).Py() << ", " << bstdHEEPEles_gsfP4_.at(iEle).Pz() << ")" << std::endl;

			//Variables storing the heep::Ele method values - Classification...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       classification=" << bstdHEEPEles_classification_.at(iEle);
			std::cout << "; isEcalDriven=" << bstdHEEPEles_isEcalDriven_.at(iEle);
			std::cout << "; isTrakerDriven=" << bstdHEEPEles_isTrackerDriven_.at(iEle) << std::endl;
			std::cout << "       isEB=" << bstdHEEPEles_isEB_.at(iEle);
			std::cout << "; isEE=" << bstdHEEPEles_isEE_.at(iEle) <<std::endl;

			//Variables storing the heep::Ele method values - track methods ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       charge=" << bstdHEEPEles_charge_.at(iEle);
			std::cout << "; trkCharge=" << bstdHEEPEles_trkCharge_.at(iEle) << std::endl;
			std::cout << "       pVtx=" << bstdHEEPEles_pVtx_.at(iEle);
			std::cout << "; pCalo=" << bstdHEEPEles_pCalo_.at(iEle);
			std::cout << "; ptVtx=" << bstdHEEPEles_ptVtx_.at(iEle);
			std::cout << "; ptCalo=" << bstdHEEPEles_ptCalo_.at(iEle) << std::endl;

			//Variables storing the heep::Ele method values ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       hOver=" << bstdHEEPEles_hOverE_.at(iEle);
			std::cout << "; dEtaIn=" << bstdHEEPEles_dEtaIn_.at(iEle);
			std::cout << "; dPhiIn=" << bstdHEEPEles_dPhiIn_.at(iEle);
			std::cout << "; dPhiOut=" << bstdHEEPEles_dPhiOut_.at(iEle) << std::endl;
			std::cout << "       epIn=" << bstdHEEPEles_epIn_.at(iEle);
			std::cout << "; epOut=" << bstdHEEPEles_epOut_.at(iEle);
			std::cout << "; fbrem=" << bstdHEEPEles_fbrem_.at(iEle);
			std::cout << "; bremFrac=" << bstdHEEPEles_bremFrac_.at(iEle);
			std::cout << "; invEOverInvP=" << bstdHEEPEles_invEOverInvP_.at(iEle) << std::endl;

			//Variables storing the heep::Ele method values - shower shape variables ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       sigmaEtaEta=" << bstdHEEPEles_sigmaEtaEta_.at(iEle);
			std::cout << "; sigmaEtaEtaUnCorr=" << bstdHEEPEles_sigmaEtaEtaUnCorr_.at(iEle);
			std::cout << "; sigmaIEtaIEta=" << bstdHEEPEles_sigmaIEtaIEta_.at(iEle) << std::endl;
			std::cout << "       e1x5=" << bstdHEEPEles_e1x5_.at(iEle);
			std::cout << "; e2x5Max=" << bstdHEEPEles_e2x5Max_.at(iEle);
			std::cout << "; e5x5=" << bstdHEEPEles_e5x5_.at(iEle);
			std::cout << "; e1x5Over5x5=" << bstdHEEPEles_e1x5Over5x5_.at(iEle);
			std::cout << "; e2x5MaxOver5x5=" << bstdHEEPEles_e2x5MaxOver5x5_.at(iEle) << std::endl;

			//Variables storing the heep::Ele method values - isolation variables ...
			std::cout << "         -=-=-" << std::endl;
			std::cout << "       isolEm=" << bstdHEEPEles_isolEm_.at(iEle);
			std::cout << "; isolHad=" << bstdHEEPEles_isolHad_.at(iEle);
			std::cout << "; isolHadDepth1=" << bstdHEEPEles_isolHadDepth1_.at(iEle);
			std::cout << "; isolHadDepth2=" << bstdHEEPEles_isolHadDepth2_.at(iEle) << std::endl;
			std::cout << "       isolPtTrks" << bstdHEEPEles_isolPtTrks_.at(iEle);
			std::cout << "; isolEmHadDepth1=" << bstdHEEPEles_isolEmHadDepth1_.at(iEle) << std::endl;


		}
	}
	
}

//----------------------------------------------------------------------------------------------------
//-------- Method for reading in the muon information, and dumping it into tsw::Event class ----------
void BstdZeeNTupler::ReadInMuons(bool beVerbose, const edm::Event& edmEvent){
	// Declarations ...
	edm::Handle <reco::TrackToTrackMap> tevMapH1, tevMapH2, tevMapH3;
	edm::Handle<reco::MuonCollection> MuCollection;

	//edm::Handle<reco::MuonCollection> refitMuonCollnH;
	std::string MuonTags_ = "muons"; // It is suggested at the start of the "Available information" section of the "WorkBookMuonAnalysis" TWiki page
												// that the vector<reco::Muon> with label "muons" should be used for the high-pT muons (along with other global muons, as well as the stand-alone and tracker muons).
	reco::MuonCollection::const_iterator imuon;
	//edm::View<reco::Muon>::const_iterator imuon;

	//Getting the TrackToTrackMap's that link the refitted tracks to the corresponding global muon
	/*edmEvent.getByLabel("tevMuons", "default", tevMapH1);
	const reco::TrackToTrackMap tevMap1 = *(tevMapH1.product());
	edmEvent.getByLabel("tevMuons", "firstHit", tevMapH2);
	const reco::TrackToTrackMap tevMap2 = *(tevMapH2.product());
	edmEvent.getByLabel("tevMuons", "picky", tevMapH3);
	const reco::TrackToTrackMap tevMap3 = *(tevMapH3.product());*/
	edmEvent.getByLabel(MuonTags_, MuCollection);
	const reco::MuonCollection muonC = *(MuCollection.product());

	// Grab beamspot and highest pt-sum vertex
	edm::Handle<reco::BeamSpot> h_beamSpot;
	edmEvent.getByLabel("offlineBeamSpot", h_beamSpot);
	if( ! h_beamSpot.isValid() )
		edm::LogWarning("BstdZeeNTupler") << "WARNING!! Error will occur soon since handle to beam spot is not valid!\n";

	edm::Handle<reco::VertexCollection> h_vertices;
	edmEvent.getByLabel("offlinePrimaryVertices", h_vertices);
	if( ! h_vertices.isValid() )
		edm::LogWarning("BstdZeeNTupler") << "WARNING!! Error will occur soon since hande to primary vertices is not valid!\n";
//	const reco::Vertex& firstPrimaryVertex = ;//( *h_vertices->begin() );
	std::vector<reco::Vertex>::const_iterator mainPrimaryVertexIt(h_vertices->begin());
	double mainPrimaryVtx_sumTrackPts = 0.0;

	for(std::vector<reco::Vertex>::const_iterator vertexIt = h_vertices->begin();
						vertexIt != h_vertices->end(); vertexIt++ )
	{
		// Apply basic quality cuts to vertex ...
		bool goodQualityVtx = (vertexIt->tracksSize()>=4) && ( vertexIt->z()<24. && sqrt( pow(vertexIt->x(),2.0)+pow(vertexIt->y(),2.0) )<2.0 );
		if( goodQualityVtx )
		{
			double ithVtx_sumTrackPts = 0.0;
			for(std::vector<reco::Track>::const_iterator trkIt=vertexIt->refittedTracks().begin(); trkIt!=vertexIt->refittedTracks().end(); trkIt++)
				ithVtx_sumTrackPts += sqrt(trkIt->innerMomentum().perp2());

			if(ithVtx_sumTrackPts > mainPrimaryVtx_sumTrackPts)
				mainPrimaryVertexIt = vertexIt;
		}
	}

	if( mainPrimaryVertexIt!=h_vertices->begin() )
		edm::LogWarning("BstdZeeNTupler") << "WARNING!! The highest sum-pT primary vertex used in " << __func__ << "is not the first from the collection!\n  This is unexpected ... \n";

	// Running over the muons in the collection ...
	if(beVerbose){std::cout << " ->There are " << muonC.size() << " muons in this event." << std::endl;}
	for(imuon = muonC.begin(); imuon != muonC.end(); ++imuon){
		tsw::MuStruct ithMuon;

		//General muon information ...
		ithMuon.p4                 = imuon->p4();
		ithMuon.charge             = imuon->charge();
		ithMuon.isGlobalMuon       = imuon->isGlobalMuon();
		ithMuon.isTrackerMuon      = imuon->isTrackerMuon();
		ithMuon.isStandAloneMuon   = imuon->isStandAloneMuon();
		ithMuon.numMatchedMuonStns = imuon->numberOfMatchedStations();
		ithMuon.isolR03_sumPt      = imuon->isolationR03().sumPt;
		// N.B. The track method is not used below since it just returns the inner track - i.e. it just returns the output of the innerTrack() method

		// Global track information ...
		if(imuon->globalTrack().get()!=0){
			ithMuon.globTrk_exists = true;
			ithMuon.globTrk_pT            = imuon->globalTrack()->pt();
			ithMuon.globTrk_eta           = imuon->globalTrack()->eta();
			ithMuon.globTrk_phi           = imuon->globalTrack()->phi();
			ithMuon.globTrk_charge        = imuon->globalTrack()->charge();
			ithMuon.globTrk_numberOfValidMuonHits = imuon->globalTrack()->hitPattern().numberOfValidMuonHits();
			ithMuon.globTrk_normalisedChi2 = imuon->globalTrack()->normalizedChi2();
		}
		else
			ithMuon.globTrk_exists = false;

		// Inner track information ...
		if(imuon->innerTrack().get()!=0){
			ithMuon.inTrk_exists = true;
			ithMuon.inTrk_pT               = imuon->innerTrack()->pt();
			ithMuon.inTrk_eta              = imuon->innerTrack()->eta();
			ithMuon.inTrk_phi              = imuon->innerTrack()->phi();
			ithMuon.inTrk_charge           = imuon->innerTrack()->charge();
			ithMuon.inTrk_numValidPixHits  = imuon->innerTrack()->hitPattern().numberOfValidPixelHits();
			ithMuon.inTrk_numValidTrkrHits = imuon->innerTrack()->hitPattern().numberOfValidTrackerHits();
			ithMuon.inTrk_dxyVsOrigin      = imuon->innerTrack()->dxy(); // Really, this should be calculated as innerTrack()->dxy(vertex->position()), but no vertex information is read in at the moment, and MuonRecoPerformance2010 TWiki page => Can calculate this approximately just relative to (0,0,0)
			ithMuon.trk_trkrLayersWHits    = imuon->track()->hitPattern().trackerLayersWithMeasurement();
		}
		else
			ithMuon.inTrk_exists = false;

		// Outer track information ...
		if(imuon->outerTrack().get()!=0){
			ithMuon.outTrk_exists = true;
			ithMuon.outTrk_pT     = imuon->outerTrack()->pt();
			ithMuon.outTrk_eta    = imuon->outerTrack()->eta();
			ithMuon.outTrk_phi    = imuon->outerTrack()->phi();
			ithMuon.outTrk_charge = imuon->outerTrack()->charge();
		}
		else
			ithMuon.outTrk_exists = false;

		if( imuon->muonBestTrack().get()!=0 ){
			ithMuon.bestTrk_exists = true;
			ithMuon.bestTrk_dxy_bspot   = imuon->muonBestTrack()->dxy( h_beamSpot->position() );
			ithMuon.bestTrk_dxy_vtx     = imuon->muonBestTrack()->dxy( mainPrimaryVertexIt->position() );
			ithMuon.bestTrk_dz_vtx      = imuon->muonBestTrack()->dz( mainPrimaryVertexIt->position() );
		}
		else
			ithMuon.bestTrk_exists = false;

		//Printing this information to screen (if desired) ...
		if(beVerbose){
			//General muon information ...
			std::cout << "     Muon ...: charge=" << ithMuon.charge << std::endl;
			std::cout << "         p=" << ithMuon.p4.P() << "; (pT, eta, phi) = " << ithMuon.p4.pt() << ", " << ithMuon.p4.eta() << ", " << ithMuon.p4.phi() << ")" << std::endl;
			std::cout << "         is{Global,Tracker,StandAlone}Muon = {" << ithMuon.isGlobalMuon << ", " << ithMuon.isTrackerMuon << ", " << ithMuon.isStandAloneMuon << "}" << std::endl;
			std::cout << "         numberOfMatchedMuonStations = " << ithMuon.numMatchedMuonStns << std::endl;
			std::cout << "         isolationR03().sumPt = " << ithMuon.isolR03_sumPt << std::endl;
			// N.B. No use of the track method since it just returns the inner track - i.e. it just returns the output of the innerTrack() method

			// Global track information ...
			std::cout << "        *GlobTrk:  " << std::endl;
			if(ithMuon.globTrk_exists){
				std::cout << "            (pT, eta, phi) = (" << ithMuon.globTrk_pT << ", " << ithMuon.globTrk_eta << ", " << ithMuon.globTrk_phi << "); charge=" << ithMuon.globTrk_charge << std::endl;
				std::cout << "            numberOfValidMuonHits=" << ithMuon.globTrk_numberOfValidMuonHits << "; normalizedChi2=" << ithMuon.globTrk_normalisedChi2 << std::endl;
			}
			else
				std::cout << "            *** globalTrack() method => NULL pointer. ***" << std::endl;

			// Inner track information ...
			std::cout << "        *InnerTrk: " << std::endl;
			if(ithMuon.inTrk_exists){
				std::cout << "            (pT, eta, phi) = (" << ithMuon.inTrk_pT << ", " << ithMuon.inTrk_eta << ", " << ithMuon.inTrk_phi  << "); charge=" << ithMuon.inTrk_charge << std::endl;
				std::cout << "            NumPixHits=" << ithMuon.inTrk_numValidPixHits << "; Num hits in trkr=" << ithMuon.inTrk_numValidTrkrHits << std::endl;
				std::cout << "            d_xy rel to origin = " << ithMuon.inTrk_dxyVsOrigin << std::endl; // Really, this should be calculated as innerTrack()->dxy(vertex->position()), but no vertex information is read in at the moment, and MuonRecoPerformance2010 TWiki page => Can calculate this approximately just relative to (0,0,0)
			}
			else
				std::cout << "            *** innerTrack() method => NULL pointer. ***" << std::endl;

			// Outer track information ...
			std::cout << "        *OuterTrk:" << std::endl;
			if(ithMuon.outTrk_exists){
				std::cout << "            (pT, eta, phi) = (" << ithMuon.outTrk_pT  << ", " << ithMuon.outTrk_eta << ", " << ithMuon.outTrk_phi << ")" << std::endl;
				std::cout << "            charge=" << ithMuon.outTrk_charge << std::endl;
			}
			else
				std::cout << "            *** outerTrack() method => NULL pointer. ***" << std::endl;
		}

		// And finally, adding this information to muon vectors in tsw::Event branch ...
		event_->AddNormMuon( &ithMuon );

	}

}

//----------------------------------------------------------------------------------------------------
//------------ method for setting up the values of the MC truth Z candidate variables... -------------
void
BstdZeeNTupler::SetMCZcandidateVariables(bool beVerbose)
{
	const Double_t value_pi = 3.14159265;

	//Setting the value of the Z candidate 4 momentum, and simple associated kinematic variables
	mcZcandidate_p4 = mcEles_HighestEt_p4 + mcEles_2ndHighestEt_p4;

	mcZcandidate_pt   = mcZcandidate_p4.Pt();
	mcZcandidate_eta  = mcZcandidate_p4.Eta();
	mcZcandidate_phi  = mcZcandidate_p4.Phi();
	mcZcandidate_mass = mcZcandidate_p4.mass();

	//Setting the variable values for dEta, dPhi and dR of the pair of eles that make up the Z candidate
   mcZcandidate_dEtaEles = std::abs(mcEles_HighestEt_p4.Eta() - mcEles_2ndHighestEt_p4.Eta());
      
   mcZcandidate_dPhiEles = mcEles_HighestEt_p4.Phi() - mcEles_2ndHighestEt_p4.Phi();  //Here dPhi can range from -2pi to +2pi
   mcZcandidate_dPhiEles = std::abs(mcZcandidate_dPhiEles);//Now dPhi can range from 0 to 2pi
   if (mcZcandidate_dPhiEles >= value_pi){
      mcZcandidate_dPhiEles = 2.0*value_pi - mcZcandidate_dPhiEles;        //Now dPhi can only range from 0 to pi
   }

	mcZcandidate_dREles = sqrt(mcZcandidate_dEtaEles*mcZcandidate_dEtaEles + mcZcandidate_dPhiEles*mcZcandidate_dPhiEles);

	//Similarly for the lab frame opening angle of the ele pair...
	mcZcandidate_openingAngle = CalcOpeningAngle(mcEles_HighestEt_p4,mcEles_2ndHighestEt_p4);	

	if(beVerbose==true)
	{
		std::cout << " ->For the Z candidate formed from the highest two Et final state electrons in this event..." << std::endl;
		std::cout << "       P=" << mcZcandidate_p4.P() << "; Pt=" << mcZcandidate_pt << "; inv mass=" << mcZcandidate_mass << std::endl;
		std::cout << "       eta=" << mcZcandidate_eta << "; phi=" << mcZcandidate_phi << std::endl;
		std::cout << "       dEta=" << mcZcandidate_dEtaEles << "; dPhi=" << mcZcandidate_dPhiEles << "; dR=" << mcZcandidate_dREles << std::endl;
		std::cout << "       lab frame opening angle = " << mcZcandidate_openingAngle << std::endl;
	}
}

//***--**--------------------------------------------------------------------------------------**--***

Double_t CalcOpeningAngle(const ROOT::Math::XYZTVector& p4a, const ROOT::Math::XYZTVector& p4b)
{

	Double_t momentaDotProduct = 0.0;
	Double_t openingAngle = -999.9;

	momentaDotProduct  = p4a.Px()*p4b.Px();
   momentaDotProduct += p4a.Py()*p4b.Py();
   momentaDotProduct += p4a.Pz()*p4b.Pz();
   openingAngle = acos( momentaDotProduct/(p4a.P()*p4b.P()) );
   
	return openingAngle;

}
