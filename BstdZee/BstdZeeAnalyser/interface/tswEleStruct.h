#ifndef tswEleStruct_h
#define tswEleStruct_h

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h"

namespace tsw {
	// Forward declaration
	struct ModEleIsoWithPhantom;

	// Main ele data struct
   struct EleStruct {
   	EleStruct();

	   // Kinematic and geometric variables
   	float et_;
   	float gsfEt_;
   	float scEt_;
   	float energy_;
   	float gsfEnergy_;
   	float caloEnergy_;
   	float ecalEnergyError_;
   	float eta_;
   	float scEta_;
   	float detEta_;
   	float detEtaAbs_;
   	float phi_;
   	float scPhi_;
   	float detPhi_;
   	float zVtx_;
   	ROOT::Math::XYZTVector p4_;
   	ROOT::Math::XYZTVector gsfP4_;

   	// 'Classification'
   	int classification_;
   	bool isEcalDriven_;   	bool isTrackerDriven_;
   	bool isEB_;
   	bool isEE_;

   	// Track variables ...
   	int charge_;
   	int trkCharge_;
   	float pVtx_;
	  	float pCalo_;
	  	float ptVtx_;
	  	float ptCalo_;
	  	float closestCtfTrk_pt_; //The closestCtfTrk variables are all assigned the value -999.9 if the pointer to this track was a NULL pointer
	  	float closestCtfTrk_eta_;
	  	float closestCtfTrk_phi_;
	  	float closestCtfTrk_innerPt_;
	  	float closestCtfTrk_innerEta_;
	  	float closestCtfTrk_innerPhi_;
	  	float closestCtfTrk_outerPt_;
	  	float closestCtfTrk_outerEta_;
	  	float closestCtfTrk_outerPhi_;

	  	// Various other variables ...
	  	float hOverE_;
	  	float dEtaIn_;
	  	float dPhiIn_;
	  	float dPhiOut_;
	  	float epIn_;
	  	float epOut_;
	  	float fbrem_;
	  	float bremFrac_;
	  	float invEOverInvP_;

	  	// Shower shape variables
	  	float sigmaEtaEta_;
	  	float sigmaEtaEtaUnCorr_;
	  	float sigmaIEtaIEta_;
	  	float e1x5_;
	  	float e2x5Max_;
	  	float e5x5_;
	  	float e1x5Over5x5_;
	  	float e2x5MaxOver5x5_;

	  	// Isolation variables ...
	  	float isolEm_;
	  	float isolHad_;
	  	float isolHadDepth1_;
	  	float isolHadDepth2_;
	  	float isolPtTrks_;
	  	float isolEmHadDepth1_;

	  	// Alternative isolation variables, calculated by IsoDeps modules ...
	  	double isol_isoDep_stdTrk_;
	  	double isol_isoDep_stdEm_;
	  	double isol_isoDep_stdHadD1_;
	  	double isol_isoDep_inrVetoModTrk_;
	  	double isol_isoDep_inrVetoModEm_;
	  	double isol_isoDep_inrVetoModHadD1_;

	  	// Alternative isolation variables, calculated by BstdZeeModIsolProducer ...
	  	double isol_inrXSVetoModTrk_;
	  	double isol_inrSVetoModTrk_;
	  	double isol_inrMVetoModTrk_;
	  	double isol_inrLVetoModTrk_;
	  	double isol_inrXLVetoModTrk_;

	  	double isol_inrXSVetoModEm_;
	  	double isol_inrSVetoModEm_;
	  	double isol_inrMVetoModEm_;
	  	double isol_inrLVetoModEm_;
	  	double isol_inrXLVetoModEm_;

	  	double isol_inrXSVetoModHadD1_;
	  	double isol_inrSVetoModHadD1_;
	  	double isol_inrMVetoModHadD1_;
	  	double isol_inrLVetoModHadD1_;
	  	double isol_inrXLVetoModHadD1_;

	  	double isol_inrVetoModTrk_otherEleAreaForSelf_;
	  	double isol_inrVetoModEcal_otherEleAreaForSelf_;
	  	double isol_inrVetoModHcalD1_otherEleAreaForSelf_;

	  	std::vector<tsw::ModEleIsoWithPhantom> isol_inrVetoModIsosWithPhantomEle_;

	  	// Number & pT sum of final-state gen-level hadrons within dR cone of 0.4
	  	unsigned int isol_nGenHadronsDr04_;
	  	double isol_ptSumGenHadronsDr04_;

	  	// Number of missing hits
	  	unsigned int numMissInnerHits_;
	   double dxy_;

	  	// Information about this ele's SC, and it's recHits ...
	  	float SC_posn_eta_;
	  	float SC_posn_phi_;
	  	float SC_rawEnergy_;
	  	std::vector<float> SC_recHits_Et_;
	  	std::vector<float> SC_recHits_eta_;
	  	std::vector<float> SC_recHits_phi_;
	  	std::vector<bool>  SC_recHits_isFromEB_;
	  	float SC_totEnergyRecHits_;
	  	unsigned int SC_totNumRecHits_;

	  	// Electron's GSF track eta, phi and vz values
	  	float gsfTrk_eta_;
	  	float gsfTrk_phi_;
	  	float gsfTrk_vz_;

	  	// Information about CTF tracks which fall into the inner isolation cone area of this electron
	  	std::vector<float> innerIsoConeTrks_pt_;
	  	std::vector<float> innerIsoConeTrks_eta_;
	  	std::vector<float> innerIsoConeTrks_phi_;
	  	std::vector<float> innerIsoConeTrks_vz_;
   };

   EleStruct::EleStruct() :
   		et_(-999.9), gsfEt_(-999.9), scEt_(-999.9), energy_(-999.9), gsfEnergy_(-999.9), caloEnergy_(-999.9), ecalEnergyError_(-999.9),
   		eta_(-999.9), scEta_(-999.9), detEta_(-999.9), detEtaAbs_(-999.9), phi_(-999.9), scPhi_(-999.9), detPhi_(-999.9), zVtx_(-999.9),
			classification_(false), isEcalDriven_(false), isTrackerDriven_(false), isEB_(false), isEE_(false),
			charge_(0), trkCharge_(0), pVtx_(-999.9), pCalo_(-999.9), ptVtx_(-999.9), ptCalo_(-999.9),
			closestCtfTrk_pt_(-999.9), closestCtfTrk_eta_(-999.9), closestCtfTrk_phi_(-999.9), closestCtfTrk_innerPt_(-999.9), closestCtfTrk_innerEta_(-999.9), closestCtfTrk_innerPhi_(-999.9), closestCtfTrk_outerPt_(-999.9), closestCtfTrk_outerEta_(-999.9), closestCtfTrk_outerPhi_(-999.9),
			hOverE_(-999.9),  dEtaIn_(-999.9), dPhiIn_(-999.9), dPhiOut_(-999.9), epIn_(-999.9), epOut_(-999.9), fbrem_(-999.9), bremFrac_(-999.9), invEOverInvP_(-999.9),
			sigmaEtaEta_(-999.9), sigmaEtaEtaUnCorr_(-999.9), sigmaIEtaIEta_(-999.9), e1x5_(-999.9), e2x5Max_(-999.9), e5x5_(-999.9), e1x5Over5x5_(-999.9), e2x5MaxOver5x5_(-999.9),
			isolEm_(-999.9), isolHad_(-999.9), isolHadDepth1_(-999.9), isolHadDepth2_(-999.9), isolPtTrks_(-999.9), isolEmHadDepth1_(-999.9),
			isol_isoDep_stdTrk_(-999.9), isol_isoDep_stdEm_(-999.9), isol_isoDep_stdHadD1_(-999.9),
			isol_isoDep_inrVetoModTrk_(-999.9), isol_isoDep_inrVetoModEm_(-999.9), isol_isoDep_inrVetoModHadD1_(-999.9),
			isol_inrXSVetoModTrk_(-999.9),isol_inrSVetoModTrk_(-999.9),isol_inrMVetoModTrk_(-999.9),isol_inrLVetoModTrk_(-999.9),isol_inrXLVetoModTrk_(-999.9),
			isol_inrXSVetoModEm_(-999.9),isol_inrSVetoModEm_(-999.9),isol_inrMVetoModEm_(-999.9),isol_inrLVetoModEm_(-999.9),isol_inrXLVetoModEm_(-999.9),
			isol_inrXSVetoModHadD1_(-999.9),isol_inrSVetoModHadD1_(-999.9),isol_inrMVetoModHadD1_(-999.9),isol_inrLVetoModHadD1_(-999.9),isol_inrXLVetoModHadD1_(-999.9),
			isol_nGenHadronsDr04_(9999), isol_ptSumGenHadronsDr04_(9999.9),
			numMissInnerHits_(9999), dxy_(-999.9),
			SC_posn_eta_(-999.9), SC_posn_phi_(-999.9), SC_rawEnergy_(-999.9), SC_recHits_Et_(), SC_recHits_eta_(), SC_recHits_phi_(), SC_recHits_isFromEB_(), SC_totEnergyRecHits_(-999.9), SC_totNumRecHits_(0),
			gsfTrk_eta_(-999.9), gsfTrk_phi_(-999.9), gsfTrk_vz_(-999.9),
			innerIsoConeTrks_pt_(), innerIsoConeTrks_eta_(), innerIsoConeTrks_phi_(), innerIsoConeTrks_vz_()
	{
   	ROOT::Math::XYZTVector defaultP4(0.0, 0.0, 0.0, 0.0);
		p4_    = defaultP4;
		gsfP4_ = defaultP4;
	}

   struct ModEleIsoWithPhantom{
   public:
   	ModEleIsoWithPhantom() :
   		dEta(-999.9), dPhi(-999.9),
   		trk(-999.9), ecal(-999.9), hcalD1(-999.9)
   	{  }
   	~ModEleIsoWithPhantom() { }

   	double dR() const { return sqrt(dEta*dEta + dPhi*dPhi); }

   	double dEta, dPhi;
   	double trk, ecal, hcalD1;
   };

}

#endif
