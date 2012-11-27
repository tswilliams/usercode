#ifndef tswMuStruct_h
#define tswMuStruct_h

#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h"

namespace tsw {
   struct MuStruct {
   	MuStruct() : p4(-999.9, 0.0, 0.0, 0.0), charge(-999), isGlobalMuon(false), isTrackerMuon(false), isStandAloneMuon(false), numMatchedMuonStns(-999), isolR03_sumPt(-999.9),
   			globTrk_exists(false), globTrk_pT(-999.9), globTrk_eta(-999.9), globTrk_phi(-999.9), globTrk_charge(-999), globTrk_numberOfValidMuonHits(-999), globTrk_normalisedChi2(-999.9),
   			inTrk_exists(false), inTrk_pT(-999.9), inTrk_eta(-999.9), inTrk_phi(-999.9), inTrk_charge(-999), inTrk_numValidPixHits(-999), inTrk_numValidTrkrHits(-999), inTrk_dxyVsOrigin(-999.9), trk_trkrLayersWHits(-999),
   			outTrk_exists(false), outTrk_pT(-999.9), outTrk_eta(-999.9), outTrk_phi(-999.9), outTrk_charge(-999),
   			bestTrk_exists(false), bestTrk_pT(-999.9), bestTrk_ptError(-999.9), bestTrk_dxy_bspot(-999.9), bestTrk_dxy_vtx(-999.9), bestTrk_dz_vtx(-999.9)
   		{}
   	// General variables ...
   	ROOT::Math::XYZTVector p4;
   	int    charge;
   	bool   isGlobalMuon;
   	bool   isTrackerMuon;
   	bool   isStandAloneMuon;
   	int    numMatchedMuonStns;
   	double isolR03_sumPt;

   	// Kinematic variables from the global track ...
   	bool   globTrk_exists;
   	double globTrk_pT;
   	double globTrk_eta;
   	double globTrk_phi;
   	int    globTrk_charge;
   	int    globTrk_numberOfValidMuonHits;
   	double globTrk_normalisedChi2;
   	// ... and from the inner track ...
   	bool   inTrk_exists;
   	double inTrk_pT;
   	double inTrk_eta;
   	double inTrk_phi;
   	int    inTrk_charge;
   	int    inTrk_numValidPixHits;
   	int    inTrk_numValidTrkrHits;
   	double inTrk_dxyVsOrigin;
   	int    trk_trkrLayersWHits;

   	// ... and from the outer track ...
   	bool   outTrk_exists;
   	double outTrk_pT;
   	double outTrk_eta;
   	double outTrk_phi;
   	int    outTrk_charge;

   	// ... and from 'best' track ( = tuneP track developed for high pT muon ID)...
   	bool   bestTrk_exists;
   	double bestTrk_pT;
   	double bestTrk_ptError;
   	double bestTrk_dxy_bspot;
   	double bestTrk_dxy_vtx;
   	double bestTrk_dz_vtx;
   };
}

#endif
