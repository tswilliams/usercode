
#include "TSWilliams/BstdZeeNTupler/interface/tswEvent.h"
#include <iostream>

namespace tsw{

	Event::Event() :
		runNum_(0), lumiSec_(0), evtNum_(0),
		mc_genWeight_(-999.9), mcLHE_ZbosonP4_(),
		mc_bx0_nPUVtxPoissonMean_(-999.9), mc_bx0_nPUVtx_(9999), mc_bx0_PUVtxZPosns_(), mc_puWeight_1D_(-999.9),
		recoVtx_totalNum_(9999), recoVtx_numGoodVtxs_(9999),
		pu_rho_(0.0),
		normMuons_charge_(), normMuons_isGlobalMuon_(), normMuons_isTrackerMuon_(), normMuons_isStandAloneMuon_(), normMuons_numMatchedMuonStns_(), normMuons_isolR03_sumPt_(),
		normMuons_globTrk_exists_(), normMuons_globTrk_pT_(), normMuons_globTrk_eta_(), normMuons_globTrk_phi_(), normMuons_globTrk_charge_(), normMuons_globTrk_numberOfValidMuonHits_(), normMuons_globTrk_normalisedChi2_(),
		normMuons_inTrk_exists_(), normMuons_inTrk_pT_(), normMuons_inTrk_eta_(), normMuons_inTrk_phi_(), normMuons_inTrk_charge_(), normMuons_inTrk_numValidPixHits_(), normMuons_inTrk_numValidTrkrHits_(), normMuons_inTrk_dxyVsOrigin_(),
		normMuons_outTrk_exists_(), normMuons_outTrk_pT_(), normMuons_outTrk_eta_(), normMuons_outTrk_phi_(), normMuons_outTrk_charge_()
	{

	}
	Event::~Event() {}

	void Event::SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber)
	{
		runNum_  = runNumber;
		lumiSec_ = lumiSection;
		evtNum_  = eventNumber;
	}
	void Event::SetMCGenWeight(float genWeight){
		mc_genWeight_ = genWeight;
	}
	void Event::SetLHEZbosonInfo(ROOT::Math::XYZTVector p4_Zboson){
		mcLHE_ZbosonP4_ = p4_Zboson;
	}
	void Event::SetMCPUInfo(float nVtxPoissonMean, unsigned int nVertices, std::vector<float> vtxZposns, double puWeight_1D)
	{
		mc_bx0_nPUVtxPoissonMean_ = nVtxPoissonMean;
		mc_bx0_nPUVtx_ = nVertices;
		mc_bx0_PUVtxZPosns_ = vtxZposns;
		mc_puWeight_1D_ = puWeight_1D;
	}
	void Event::SetRecoVtxInfo(unsigned int nRecoVtxs, unsigned int nGoodRecoVtxs)
	{
		recoVtx_totalNum_ = nRecoVtxs;
		recoVtx_numGoodVtxs_ = nGoodRecoVtxs;
	}
	void Event::SetPURho(const double rho){
		pu_rho_ = rho;	}

	void Event::PrintBasicEventInformation(){
		std::cout << "Run " << runNum_ << ", LumiSec " << lumiSec_ << ", event " << evtNum_ << std::endl;
	}
	void Event::PrintMCGenWeight()
	{
		std::cout << "  * MC gen weight = " << mc_genWeight_ << std::endl;
	}
	void Event::PrintLHEZbosonInfo(){
		std::cout << "  * LHE info:" << std::endl;
		std::cout << "       Z boson p4(px,py,pz,E) = " << "(" << mcLHE_ZbosonP4_.Px() << ", " <<  mcLHE_ZbosonP4_.Py() << ", " <<  mcLHE_ZbosonP4_.Pz() << ", " <<  mcLHE_ZbosonP4_.E() << "); pT=" << mcLHE_ZbosonP4_.Pt() << std::endl;
	}
	void Event::PrintPUVtxInfo(){
		std::cout << "  * MC PU Info ... " << std::endl;
		std::cout << "       bx0: true num PU vtxs = " << mc_bx0_nPUVtxPoissonMean_ << "; num PU vtxs = " << mc_bx0_nPUVtx_ << "; puWeight=" << mc_puWeight_1D_ << std::endl;
		std::cout << "            z positions: ";
		for(unsigned int iVtx=0; iVtx<mc_bx0_PUVtxZPosns_.size(); iVtx++)
			std::cout << mc_bx0_PUVtxZPosns_.at(iVtx) << ", ";
		std::cout << std::endl;
		std::cout << "  * Reconstructed vertices:" << std::endl;
		std::cout << "       totNum=" << recoVtx_totalNum_ << "; numGoodVtxs=" << recoVtx_numGoodVtxs_ << std::endl;
	}
	void Event::PrintPURho(){
		std::cout << "  * PU jet rho = " << pu_rho_ << std::endl;
	}

	void Event::SetEMuTriggerInfo(const std::string eMuPathName, const bool eMuPathDecision)
	{
		trg_emuPath_name_ = eMuPathName;
		trg_emuPath_decision_ = eMuPathDecision;
	}

	void Event::AddStdEleInfo_isoDep_std(const double stdTrkIso, const double stdEcalIso, const double stdHcalD1Iso)
	{
		stdEles_isoDeps_stdTrkIso_.push_back( stdTrkIso );
		stdEles_isoDeps_stdEcalIso_.push_back( stdEcalIso );
		stdEles_isoDeps_stdHcalD1Iso_.push_back( stdHcalD1Iso );
	}
	void Event::AddStdEleInfo_isoDep_inrVeto(const double modTrkIso, const double modEcalIso, const double modHcalD1Iso)
	{
		stdEles_isoDeps_inrVetoModTrkIso_.push_back( modTrkIso );
		stdEles_isoDeps_inrVetoModEcalIso_.push_back( modEcalIso );
		stdEles_isoDeps_inrVetoModHcalD1Iso_.push_back( modHcalD1Iso );
	}

	void Event::AddStdEleInfo_inrVetoModIso(const tsw::Event::InnerVetoSize vetoSize, const double modTrkIso, const double modEcalIso, const double modHcalD1Iso)
	{
		switch (vetoSize)
		{
			case tsw::Event::xSmallVeto:
				stdEles_inrVetoXSModIso_Trk_.push_back(modTrkIso);
				stdEles_inrVetoXSModIso_Ecal_.push_back(modEcalIso);
				stdEles_inrVetoXSModIso_HcalD1_.push_back(modHcalD1Iso);
				break;
			case tsw::Event::smallVeto:
				stdEles_inrVetoSModIso_Trk_.push_back(modTrkIso);
				stdEles_inrVetoSModIso_Ecal_.push_back(modEcalIso);
				stdEles_inrVetoSModIso_HcalD1_.push_back(modHcalD1Iso);
				break;
			case tsw::Event::mediumVeto:
				stdEles_inrVetoModIso_Trk_.push_back(modTrkIso);
				stdEles_inrVetoModIso_Ecal_.push_back(modEcalIso);
				stdEles_inrVetoModIso_HcalD1_.push_back(modHcalD1Iso);
				break;
			case tsw::Event::largeVeto:
				stdEles_inrVetoLModIso_Trk_.push_back(modTrkIso);
				stdEles_inrVetoLModIso_Ecal_.push_back(modEcalIso);
				stdEles_inrVetoLModIso_HcalD1_.push_back(modHcalD1Iso);
				break;
			case tsw::Event::xLargeVeto:
				stdEles_inrVetoXLModIso_Trk_.push_back(modTrkIso);
				stdEles_inrVetoXLModIso_Ecal_.push_back(modEcalIso);
				stdEles_inrVetoXLModIso_HcalD1_.push_back(modHcalD1Iso);
				break;
			default:
				std::cout << std::endl;
				std::cout << " *** ERROR : Unknown InnerVetoSize enum value passed to tsw::Event::AddStdEleInfo_inrVetoModIso method! ***" << std::endl;
				std::cout << std::endl;
				break;
		}
	}

	void Event::AddStdEleInfo_inrVetoModIsoWithPhantomEle(const std::vector<double>& dEtaPhantomEles, const std::vector<double>& dPhiPhantomEles, const std::vector<double>& trkIsos, const std::vector<double>& ecalIsos, const std::vector<double>& hcalD1Isos)
	{
		stdEles_inrVetoModIsoPhantomEles_dEta_.push_back( dEtaPhantomEles );
		stdEles_inrVetoModIsoPhantomEles_dPhi_.push_back( dPhiPhantomEles );
		stdEles_inrVetoModIsoPhantomEles_Trk_.push_back( trkIsos );
		stdEles_inrVetoModIsoPhantomEles_Ecal_.push_back( ecalIsos );
		stdEles_inrVetoModIsoPhantomEles_HcalD1_.push_back( hcalD1Isos );
	}

	void Event::AddStdEleInfo_genHadronsDr04(unsigned int nGenHadrons_dR04, double ptSumGenHadrons_dR04){
		stdEles_nGenHadronsDr04_.push_back(nGenHadrons_dR04);
		stdEles_ptSumGenHadronsDr04_.push_back(ptSumGenHadrons_dR04);
	}

	void Event::PrintStdEleInfo_isoValues(unsigned int iEle){
		std::cout << "    Iso values from isoDeposits... " << std::endl;
		std::cout << "          Std values:    Trk=" << stdEles_isoDeps_stdTrkIso_.at(iEle) << ", Ecal=" << stdEles_isoDeps_stdEcalIso_.at(iEle) << ", HcalDepth1=" << stdEles_isoDeps_stdHcalD1Iso_.at(iEle) << std::endl;
		std::cout << "          Mod (inrVeto): Trk=" << stdEles_isoDeps_inrVetoModTrkIso_.at(iEle) << ", Ecal=" << stdEles_isoDeps_inrVetoModEcalIso_.at(iEle) << ", HcalDepth1=" << stdEles_isoDeps_inrVetoModHcalD1Iso_.at(iEle) << std::endl;
		std::cout << "     InnerVetoMod iso values from bstdZee code ..." << std::endl;
		std::cout << "          Trk=" << stdEles_inrVetoModIso_Trk_.at(iEle) << ", Ecal=" << stdEles_inrVetoModIso_Ecal_.at(iEle) << ", HcalDepth1=" << stdEles_inrVetoModIso_HcalD1_.at(iEle) << std::endl;
	}

	void Event::PrintStdEleInfo_inrVetoModIsoWithPhantomEle(unsigned int iEle)
	{
		std::vector<double>::const_iterator dEtaIt = stdEles_inrVetoModIsoPhantomEles_dEta_.at(iEle).begin();
		std::vector<double>::const_iterator dPhiIt = stdEles_inrVetoModIsoPhantomEles_dPhi_.at(iEle).begin();
		std::vector<double>::const_iterator trkIsoIt    = stdEles_inrVetoModIsoPhantomEles_Trk_.at(iEle).begin();
		std::vector<double>::const_iterator ecalIsoIt   = stdEles_inrVetoModIsoPhantomEles_Ecal_.at(iEle).begin();
		std::vector<double>::const_iterator hcalD1IsoIt = stdEles_inrVetoModIsoPhantomEles_HcalD1_.at(iEle).begin();

		for( ; dEtaIt != stdEles_inrVetoModIsoPhantomEles_dEta_.at(iEle).end(); ){
			std::cout << "     Phantom with d(Eta,Phi)=(" << *dEtaIt << ", " << *dPhiIt << "): " << std::endl
						 << "       Isos: track=" << *trkIsoIt << " ; ecal=" << *ecalIsoIt << " ; hcalD1=" << *hcalD1IsoIt << std::endl;
			dEtaIt++;  dPhiIt++;
			trkIsoIt++;  ecalIsoIt++;  hcalD1IsoIt++;
		}
	}


	void Event::AddNormMuon(tsw::MuStruct* theMuon){
		normMuons_p4_                .push_back(theMuon->p4);
		normMuons_charge_            .push_back(theMuon->charge);
		normMuons_isGlobalMuon_      .push_back(theMuon->isGlobalMuon);
		normMuons_isTrackerMuon_     .push_back(theMuon->isTrackerMuon);
		normMuons_isStandAloneMuon_  .push_back(theMuon->isStandAloneMuon);
		normMuons_numMatchedMuonStns_.push_back(theMuon->numMatchedMuonStns);
		normMuons_isolR03_sumPt_     .push_back(theMuon->isolR03_sumPt);

		normMuons_globTrk_exists_.push_back(theMuon->globTrk_exists);
		normMuons_globTrk_pT_    .push_back(theMuon->globTrk_pT);
		normMuons_globTrk_eta_   .push_back(theMuon->globTrk_eta);
		normMuons_globTrk_phi_   .push_back(theMuon->globTrk_phi);
		normMuons_globTrk_charge_.push_back(theMuon->globTrk_charge);
		normMuons_globTrk_numberOfValidMuonHits_.push_back(theMuon->globTrk_numberOfValidMuonHits);
		normMuons_globTrk_normalisedChi2_.push_back(theMuon->globTrk_normalisedChi2);

		normMuons_inTrk_exists_         .push_back(theMuon->inTrk_exists);
		normMuons_inTrk_pT_              .push_back(theMuon->inTrk_pT);
		normMuons_inTrk_eta_             .push_back(theMuon->inTrk_eta);
		normMuons_inTrk_phi_             .push_back(theMuon->inTrk_phi);
		normMuons_inTrk_charge_          .push_back(theMuon->inTrk_charge);
		normMuons_inTrk_numValidPixHits_ .push_back(theMuon->inTrk_numValidPixHits);
		normMuons_inTrk_numValidTrkrHits_.push_back(theMuon->inTrk_numValidTrkrHits);
		normMuons_inTrk_dxyVsOrigin_     .push_back(theMuon->inTrk_dxyVsOrigin);
   	normMuons_trk_trkrLayersWHits.push_back(theMuon->trk_trkrLayersWHits);

		normMuons_outTrk_exists_.push_back(theMuon->outTrk_exists);
		normMuons_outTrk_pT_    .push_back(theMuon->outTrk_pT);
		normMuons_outTrk_eta_   .push_back(theMuon->outTrk_eta);
		normMuons_outTrk_phi_   .push_back(theMuon->outTrk_phi);
		normMuons_outTrk_charge_.push_back(theMuon->outTrk_charge);

		normMuons_bestTrk_exists_   .push_back(theMuon->bestTrk_exists);
   	normMuons_bestTrk_dxy_bspot_.push_back(theMuon->bestTrk_dxy_bspot);
   	normMuons_bestTrk_dxy_vtx_  .push_back(theMuon->bestTrk_dxy_vtx);
   	normMuons_bestTrk_dz_vtx_   .push_back(theMuon->bestTrk_dz_vtx);
	}
}


#if !defined(__CINT__)
  ClassImp(tsw::Event)
#endif
