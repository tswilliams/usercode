#ifndef tswEleStruct_h
#define tswEleStruct_h

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h"

namespace tsw {
   struct EleStruct {
   	EleStruct() : et_(-999.9), gsfEt_(-999.9), scEt_(-999.9), energy_(-999.9), gsfEnergy_(-999.9), caloEnergy_(-999.9), eta_(-999.9), scEta_(-999.9), detEta_(-999.9), detEtaAbs_(-999.9), phi_(-999.9), scPhi_(-999.9), detPhi_(-999.9), zVtx_(-999.9), p4_(ROOT::Math::XYZTVector(0.0,0.0,0.0,0.0)), gsfP4_(ROOT::Math::XYZTVector(0.0,0.0,0.0,0.0)),
   			classification_(false), isEcalDriven_(false), isTrackerDriven_(false), isEB_(false), isEE_(false),
   	   	charge_(0), trkCharge_(0), pVtx_(-999.9), pCalo_(-999.9), ptVtx_(-999.9), ptCalo_(-999.9),
   	   	hOverE_(-999.9),  dEtaIn_(-999.9), dPhiIn_(-999.9), dPhiOut_(-999.9), epIn_(-999.9), epOut_(-999.9), fbrem_(-999.9), bremFrac_(-999.9), invEOverInvP_(-999.9),
   	   	sigmaEtaEta_(-999.9), sigmaEtaEtaUnCorr_(-999.9), sigmaIEtaIEta_(-999.9), e1x5_(-999.9), e2x5Max_(-999.9), e5x5_(-999.9), e1x5Over5x5_(-999.9), e2x5MaxOver5x5_(-999.9),
   	   	isolEm_(-999.9), isolHad_(-999.9), isolHadDepth1_(-999.9), isolHadDepth2_(-999.9), isolPtTrks_(-999.9), isolEmHadDepth1_(-999.9)
   		  	{}
	   // Kinematic and geometric variables
   	float et_;
   	float gsfEt_;
   	float scEt_;
   	float energy_;
   	float gsfEnergy_;
   	float caloEnergy_;
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
   	bool isEcalDriven_;
   	bool isTrackerDriven_;
   	bool isEB_;
   	bool isEE_;

   	// Track variables ...
   	int charge_;
   	int trkCharge_;
   	float pVtx_;
	  	float pCalo_;
	  	float ptVtx_;
	  	float ptCalo_;

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
   };

   void SetDefaultValuesInEleStruct(tsw::EleStruct& eleStruct){

		ROOT::Math::XYZTVector defaultP4(0.0, 0.0, 0.0, 0.0);

		// Kinematic and geometric variables
		eleStruct.et_         = -999.9;
		eleStruct.gsfEt_      = -999.9;
		eleStruct.scEt_       = -999.9;
		eleStruct.energy_     = -999.9;
		eleStruct.gsfEnergy_  = -999.9;
		eleStruct.caloEnergy_ = -999.9;
		eleStruct.eta_        = -999.9;
		eleStruct.scEta_      = -999.9;
		eleStruct.detEta_     = -999.9;
		eleStruct.detEtaAbs_  = -999.9;
		eleStruct.phi_        = -999.9;
		eleStruct.scPhi_      = -999.9;
		eleStruct.detPhi_     = -999.9;
		eleStruct.zVtx_       = -999.9;
		eleStruct.p4_    = defaultP4;
		eleStruct.gsfP4_ = defaultP4;

		// 'Classification'
		eleStruct.classification_  = false;
		eleStruct.isEcalDriven_    = false;
		eleStruct.isTrackerDriven_ = false;
		eleStruct.isEB_            = false;
		eleStruct.isEE_            = false;

		// Track variables ...
		eleStruct.charge_    = 0;
		eleStruct.trkCharge_ = 0;
		eleStruct.pVtx_   = -999.9;
		eleStruct.pCalo_  = -999.9;
		eleStruct.ptVtx_  = -999.9;
		eleStruct.ptCalo_ = -999.9;

		// Various other variables ...
		eleStruct.hOverE_       = -999.9;
		eleStruct.dEtaIn_       = -999.9;
		eleStruct.dPhiIn_       = -999.9;
		eleStruct.dPhiOut_      = -999.9;
		eleStruct.epIn_         = -999.9;
		eleStruct.epOut_        = -999.9;
		eleStruct.fbrem_        = -999.9;
		eleStruct.bremFrac_     = -999.9;
		eleStruct.invEOverInvP_ = -999.9;

		// Shower shape variables
		eleStruct.sigmaEtaEta_       = -999.9;
		eleStruct.sigmaEtaEtaUnCorr_ = -999.9;
		eleStruct.sigmaIEtaIEta_     = -999.9;
		eleStruct.e1x5_              = -999.9;
		eleStruct.e2x5Max_           = -999.9;
		eleStruct.e5x5_              = -999.9;
		eleStruct.e1x5Over5x5_       = -999.9;
		eleStruct.e2x5MaxOver5x5_    = -999.9;

		// Isolation variables ...
		eleStruct.isolEm_          = -999.9;
		eleStruct.isolHad_         = -999.9;
		eleStruct.isolHadDepth1_   = -999.9;
		eleStruct.isolHadDepth2_   = -999.9;
		eleStruct.isolPtTrks_      = -999.9;
		eleStruct.isolEmHadDepth1_ = -999.9;
	}
}

#endif
