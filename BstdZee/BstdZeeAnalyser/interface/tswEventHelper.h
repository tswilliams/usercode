#ifndef tswEventHelper_h
#define tswEventHelper_h

// BstdZee includes
#include "TSWilliams/BstdZeeNTupler/interface/tswEvent.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswEleStruct.h"
#include "TSWilliams/BstdZeeAnalyser/interface/tswMuonCollection.h"

namespace tsw{
	class EventHelper{
		private:
			tsw::Event* theEvent_;
			float xsecWeight_;

		public:
			//CTORs and DTOR ...
			EventHelper() : theEvent_(0) {}
			EventHelper(tsw::Event* ptr_event) : theEvent_(ptr_event){}

			// Methods for getting information out of tsw::Event ...
			unsigned int runNum()   const  { return theEvent_->runNum_; }
			unsigned int lumiSec()  const  { return theEvent_->lumiSec_; }
			unsigned int eventNum() const  { return theEvent_->evtNum_; }

			void setXsecWeight(float xsecWeight) { xsecWeight_ = xsecWeight; }
			float xsecWeight() const { return xsecWeight_; }
			float totWeight() const {
				return (xsecWeight()*genWeight()*puWeight());  }

			float genWeight() const { return theEvent_->mc_genWeight_;	}
			TLorentzVector GetLHEZbosonP4(){
				TLorentzVector p4_Z = ConvertToTLorentzVector( &(theEvent_->mcLHE_ZbosonP4_) );
				return p4_Z; }

			float GetMCPU_TrueNumVtx()             const  { return theEvent_->mc_bx0_nPUVtxPoissonMean_;	}
			float GetMCPU_nVtx()                   const  { return theEvent_->mc_bx0_nPUVtx_; }
			std::vector<float> GetMCPU_VtxZPosns() const  { return theEvent_->mc_bx0_PUVtxZPosns_;	}
			double mcPU_1Dweight()                 const  { return theEvent_->mc_puWeight_1D_; }
			double puWeight()  const  { return mcPU_1Dweight(); }

			float GetRecoVtx_nVtxs()     const  { return theEvent_->recoVtx_totalNum_; }
			float GetRecoVtx_nGoodVtxs() const  { return theEvent_->recoVtx_numGoodVtxs_; }

			double GetPURho() const {
					return ( (theEvent_->pu_rho_)>0.0 ? theEvent_->pu_rho_ : 0.0 ); }

			std::string GetTrigInfo_eMuPath_name() const { return theEvent_->trg_emuPath_name_; }
			bool GetTrigInfo_eMuPath_decision()    const { return theEvent_->trg_emuPath_decision_; }

			const double GetNormEle_IsoDep_stdTrkIso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_stdTrkIso_.at(iEle); }
			const double GetNormEle_IsoDep_stdEcalIso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_stdEcalIso_.at(iEle); }
			const double GetNormEle_IsoDep_stdHcalD1Iso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_stdHcalD1Iso_.at(iEle); }
			const double GetNormEle_IsoDep_inrVetoModTrkIso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_inrVetoModTrkIso_.at(iEle); }
			const double GetNormEle_IsoDep_inrVetoModEcalIso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_inrVetoModEcalIso_.at(iEle); }
			const double GetNormEle_IsoDep_inrVetoModHcalD1Iso(const unsigned int iEle) const {return theEvent_->stdEles_isoDeps_inrVetoModHcalD1Iso_.at(iEle); }

			const double GetNormEle_inrVetoModTrkIso(const unsigned int iEle, const tsw::Event::InnerVetoSize vetoSize) const
			{
				double tkIso = 9999.9;

				switch (vetoSize)
				{
				case tsw::Event::xSmallVeto:
					tkIso = theEvent_->stdEles_inrVetoXSModIso_Trk_.at(iEle); break;
				case tsw::Event::smallVeto:
					tkIso = theEvent_->stdEles_inrVetoSModIso_Trk_.at(iEle); break;
				case tsw::Event::mediumVeto:
					tkIso = theEvent_->stdEles_inrVetoModIso_Trk_.at(iEle); break;
				case tsw::Event::largeVeto:
					tkIso = theEvent_->stdEles_inrVetoLModIso_Trk_.at(iEle); break;
				case tsw::Event::xLargeVeto:
					tkIso = theEvent_->stdEles_inrVetoXLModIso_Trk_.at(iEle); break;
				default:
					std::cout << " *** ERROR : Unknown InnerVetoSize enum value passed to tsw::EventHelper::GetNormEle_inrVetoModTrkIso ! ***" << std::endl;
					break;
				}

				return tkIso;
			}

			const double GetNormEle_inrVetoModEmIso(const unsigned int iEle, const tsw::Event::InnerVetoSize vetoSize) const
			{
				double ecalIso = 9999.9;

				switch (vetoSize)
				{
				case tsw::Event::xSmallVeto:
					ecalIso = theEvent_->stdEles_inrVetoXSModIso_Ecal_.at(iEle); break;
				case tsw::Event::smallVeto:
					ecalIso = theEvent_->stdEles_inrVetoSModIso_Ecal_.at(iEle); break;
				case tsw::Event::mediumVeto:
					ecalIso = theEvent_->stdEles_inrVetoModIso_Ecal_.at(iEle); break;
				case tsw::Event::largeVeto:
					ecalIso = theEvent_->stdEles_inrVetoLModIso_Ecal_.at(iEle); break;
				case tsw::Event::xLargeVeto:
					ecalIso = theEvent_->stdEles_inrVetoXLModIso_Ecal_.at(iEle); break;
				default:
					std::cout << " *** ERROR : Unknown InnerVetoSize enum value passed to tsw::EventHelper::GetNormEle_inrVetoModEcalIso ! ***" << std::endl;
					break;
				}

				return ecalIso;
			}

			const double GetNormEle_inrVetoModHadD1Iso(const unsigned int iEle, const tsw::Event::InnerVetoSize vetoSize) const
			{
				double hcalD1Iso = 9999.9;

				switch (vetoSize)
				{
				case tsw::Event::xSmallVeto:
					hcalD1Iso = theEvent_->stdEles_inrVetoXSModIso_HcalD1_.at(iEle); break;
				case tsw::Event::smallVeto:
					hcalD1Iso = theEvent_->stdEles_inrVetoSModIso_HcalD1_.at(iEle); break;
				case tsw::Event::mediumVeto:
					hcalD1Iso = theEvent_->stdEles_inrVetoModIso_HcalD1_.at(iEle); break;
				case tsw::Event::largeVeto:
					hcalD1Iso = theEvent_->stdEles_inrVetoLModIso_HcalD1_.at(iEle); break;
				case tsw::Event::xLargeVeto:
					hcalD1Iso = theEvent_->stdEles_inrVetoXLModIso_HcalD1_.at(iEle); break;
				default:
					std::cout << " *** ERROR : Unknown InnerVetoSize enum value passed to tsw::EventHelper::GetNormEle_inrVetoModHcalD1Iso ! ***" << std::endl;
					break;
				}

				return hcalD1Iso;
			}

			unsigned int GetNormEle_nGenHadronsDr04(const unsigned int iEle) const {
				return theEvent_->stdEles_nGenHadronsDr04_.at(iEle); }
			double GetNormEle_ptSumGenHadronsDr04(const unsigned int iEle) const {
				return theEvent_->stdEles_ptSumGenHadronsDr04_.at(iEle); }

			void SetEleStructValues(const unsigned int eleIdx, tsw::EleStruct& ) const;

			tsw::MuonCollection GetNormMuons();
	};
}

void tsw::EventHelper::SetEleStructValues(const unsigned int eleIdx, tsw::EleStruct& theStruct) const
{
	theStruct.dxy_ = theEvent_->stdEles_dxy_.at(eleIdx);

	theStruct.isol_inrVetoModTrk_otherEleAreaForSelf_    = theEvent_->stdEles_inrVetoModIso_otherEleAreaForSelf_Trk_.at(eleIdx);
	theStruct.isol_inrVetoModEcal_otherEleAreaForSelf_   = theEvent_->stdEles_inrVetoModIso_otherEleAreaForSelf_Ecal_.at(eleIdx);
	theStruct.isol_inrVetoModHcalD1_otherEleAreaForSelf_ = theEvent_->stdEles_inrVetoModIso_otherEleAreaForSelf_HcalD1_.at(eleIdx);

	for(unsigned int i=0; i<theEvent_->stdEles_inrVetoModIsoPhantomEles_dEta_.at(eleIdx).size(); i++){
		tsw::ModEleIsoWithPhantom phantomEleModIsoStruct;
		phantomEleModIsoStruct.dEta   = theEvent_->stdEles_inrVetoModIsoPhantomEles_dEta_.at(eleIdx).at(i);
		phantomEleModIsoStruct.dPhi   = theEvent_->stdEles_inrVetoModIsoPhantomEles_dPhi_.at(eleIdx).at(i);
		phantomEleModIsoStruct.trk    = theEvent_->stdEles_inrVetoModIsoPhantomEles_Trk_.at(eleIdx).at(i);
		phantomEleModIsoStruct.ecal   = theEvent_->stdEles_inrVetoModIsoPhantomEles_Ecal_.at(eleIdx).at(i);
		phantomEleModIsoStruct.hcalD1 = theEvent_->stdEles_inrVetoModIsoPhantomEles_HcalD1_.at(eleIdx).at(i);

		theStruct.isol_inrVetoModIsosWithPhantomEle_.push_back( phantomEleModIsoStruct );
	}
}

	tsw::MuonCollection tsw::EventHelper::GetNormMuons(){
		std::vector<tsw::Muon> vecOfMuons; vecOfMuons.clear();
		unsigned int numOfMuons = theEvent_->normMuons_charge_.size();

		//Loop over the normMuon vectors in the event, storing the information for each muon in a tsw::MuStruct ...
		//... and using the muStruct to initialise instances of tsw::Muon ...
		for(unsigned int muIdx=0; muIdx<numOfMuons; muIdx++){
			MuStruct ithMuStruct;
			// General variables ...
			ithMuStruct.p4                 = theEvent_->normMuons_p4_.at(muIdx);
			ithMuStruct.charge             = theEvent_->normMuons_charge_.at(muIdx);
			ithMuStruct.isGlobalMuon       = theEvent_->normMuons_isGlobalMuon_.at(muIdx);
			ithMuStruct.isTrackerMuon      = theEvent_->normMuons_isTrackerMuon_.at(muIdx);
			ithMuStruct.isStandAloneMuon   = theEvent_->normMuons_isStandAloneMuon_.at(muIdx);
			ithMuStruct.numMatchedMuonStns = theEvent_->normMuons_numMatchedMuonStns_.at(muIdx);
			ithMuStruct.isolR03_sumPt      = theEvent_->normMuons_isolR03_sumPt_.at(muIdx);
			// Kinematic variables from the global track ...
			ithMuStruct.globTrk_exists = theEvent_->normMuons_globTrk_exists_.at(muIdx);
			ithMuStruct.globTrk_pT     = theEvent_->normMuons_globTrk_pT_.at(muIdx);
			ithMuStruct.globTrk_eta    = theEvent_->normMuons_globTrk_eta_.at(muIdx);
			ithMuStruct.globTrk_phi    = theEvent_->normMuons_globTrk_phi_.at(muIdx);
			ithMuStruct.globTrk_charge = theEvent_->normMuons_globTrk_charge_.at(muIdx);
			ithMuStruct.globTrk_numberOfValidMuonHits = theEvent_->normMuons_globTrk_numberOfValidMuonHits_.at(muIdx);
			ithMuStruct.globTrk_normalisedChi2 = theEvent_->normMuons_globTrk_normalisedChi2_.at(muIdx);
			// ... and from the inner track ...
			ithMuStruct.inTrk_exists           = theEvent_->normMuons_inTrk_exists_.at(muIdx);
			ithMuStruct.inTrk_pT               = theEvent_->normMuons_inTrk_pT_.at(muIdx);
			ithMuStruct.inTrk_eta              = theEvent_->normMuons_inTrk_eta_.at(muIdx);
			ithMuStruct.inTrk_phi              = theEvent_->normMuons_inTrk_phi_.at(muIdx);
			ithMuStruct.inTrk_charge           = theEvent_->normMuons_inTrk_charge_.at(muIdx);
			ithMuStruct.inTrk_numValidPixHits  = theEvent_->normMuons_inTrk_numValidPixHits_.at(muIdx);
			ithMuStruct.inTrk_numValidTrkrHits = theEvent_->normMuons_inTrk_numValidTrkrHits_.at(muIdx);
			ithMuStruct.inTrk_dxyVsOrigin      = theEvent_->normMuons_inTrk_dxyVsOrigin_.at(muIdx);
			ithMuStruct.trk_trkrLayersWHits = theEvent_->normMuons_trk_trkrLayersWHits.at(muIdx);
	   	// ... and from the outer track ...
			ithMuStruct.outTrk_exists = theEvent_->normMuons_outTrk_exists_.at(muIdx);
			ithMuStruct.outTrk_pT     = theEvent_->normMuons_outTrk_pT_.at(muIdx);
			ithMuStruct.outTrk_eta    = theEvent_->normMuons_outTrk_eta_.at(muIdx);
			ithMuStruct.outTrk_phi    = theEvent_->normMuons_outTrk_phi_.at(muIdx);
			ithMuStruct.outTrk_charge = theEvent_->normMuons_outTrk_charge_.at(muIdx);

			ithMuStruct.bestTrk_exists    = theEvent_->normMuons_bestTrk_exists_.at(muIdx);
			ithMuStruct.bestTrk_dxy_bspot = theEvent_->normMuons_bestTrk_dxy_bspot_.at(muIdx);
			ithMuStruct.bestTrk_dxy_vtx   = theEvent_->normMuons_bestTrk_dxy_vtx_.at(muIdx);
			ithMuStruct.bestTrk_dz_vtx    = theEvent_->normMuons_bestTrk_dz_vtx_.at(muIdx);

			Muon ithMuon(ithMuStruct);
			vecOfMuons.push_back(ithMuon);
		}
		MuonCollection theMuColln(&vecOfMuons, "normal muons");
		return theMuColln;
	}

#endif
