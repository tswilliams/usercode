#ifndef tswHEEPEle_h
#define tswHEEPEle_h

// C++ includes
#include <vector>
#include <functional>

// ROOT includes
#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h"

// BstdZee includes
#include "TSWilliams/BstdZeeAnalyser/interface/tswEleStruct.h"
#include "TSWilliams/BstdZeeNTupler/interface/tswUsefulFunctions.h"


namespace tsw{
	class HEEPEle{
		friend class HEEPDiEle;
		private:
			EleStruct eleStr_;
		public:
			HEEPEle();
			HEEPEle(tsw::EleStruct eleStruct);
			~HEEPEle();

			void PrintOutVariables();

			bool ApplySimpleCuts();
			bool isFiducialAndEcalDriven();
			const bool isHEEPEB() const;
			const bool isHEEPEE() const;
		private:
			const bool ApplyHEEPIsoCut_Trk(const double trkIsoValue) const {
				bool tmpCutsFlag = false;
				if( isHEEPEB() || isHEEPEE() )
					tmpCutsFlag = (trkIsoValue<5.0);
				return tmpCutsFlag;
			}
			const bool ApplyHEEPIsoCut_EmHad1(float ) const;
		public:

			// PUBLIC STRUCTS //
			enum ModIsoPhantomEleDrRange{
				PHANTOM_DR_005_010=0,
				PHANTOM_DR_010_015=1,
				PHANTOM_DR_015_020=2,
				PHANTOM_DR_020_025=3,
				PHANTOM_DR_025_030=4,
				PHANTOM_DR_030_035=5,
				PHANTOM_DR_035_040=6
			};

			/////////////////////////////////////////////////
			// Methods for reading out all of the HEEP variables ...
			// Kinematic and geometric variables
			const float et()  const {return eleStr_.et_;}
			float gsfEt()     const {return eleStr_.gsfEt_;}
			float scEt()      const {return eleStr_.scEt_;}
			float energy()    const {return eleStr_.energy_;}
//			float gsfEnergy(){return eleStr_.gsfEnergy_;}
			float caloEnergy() const {return eleStr_.caloEnergy_;}
			float ecalEnergyError(){return eleStr_.ecalEnergyError_;}
			float eta()         const {return eleStr_.eta_;}
			const float scEta() const {return eleStr_.scEta_;}
//			float detEta(){return eleStr_.detEta_;}
//			float detEtaAbs(){return eleStr_.detEtaAbs_;}
			float phi()    const  {return eleStr_.phi_;}
			float scPhi()  const  {return eleStr_.scPhi_;}
//			float detPhi(){return eleStr_.detPhi_;}
//			float zVtx(){return eleStr_.zVtx_;}
			TLorentzVector p4()    const  {return ConvertToTLorentzVector( &(eleStr_.p4_) );}
			TLorentzVector gsfP4() const  {return ConvertToTLorentzVector( &(eleStr_.gsfP4_) );}

			// 'Classification'
//			int classification(){return eleStr_.classification_;}
			bool isEcalDriven() const {return eleStr_.isEcalDriven_;}
//			bool isTrackerDriven(){return eleStr_.isTrackerDriven_;}
			bool isEB() const {return eleStr_.isEB_;}
			bool isEE() const {return eleStr_.isEE_;}

			// Track variables ...
			int charge() const {return eleStr_.charge_;}
//			int trkCharge(){return eleStr_.trkCharge_;}
//			float pVtx(){return eleStr_.pVtx_;}
//			float pCalo(){return eleStr_.pCalo_;}
			float ptVtx(){return eleStr_.ptVtx_;}
			float ptCalo(){return eleStr_.ptCalo_;}
			bool  closestCtfTrk_exists() const {
				if(closestCtfTrk_pt()<-0.1){
					//std::cout << std::endl << "    ***** closestCtfTrackRef is NULL *****" << std::endl << std::endl;
					return false;}
				else
					return true;
			}
		  	float closestCtfTrk_pt()  const {return eleStr_.closestCtfTrk_pt_;}
		  	float closestCtfTrk_eta() const {return eleStr_.closestCtfTrk_eta_;}
		  	float closestCtfTrk_phi() const {return eleStr_.closestCtfTrk_phi_;}
//		  	float closestCtfTrk_innerPt(){return eleStr_.closestCtfTrk_innerPt_;}
//		  	float closestCtfTrk_innerEta(){return eleStr_.closestCtfTrk_innerEta_;}
//		  	float closestCtfTrk_innerPhi(){return eleStr_.closestCtfTrk_innerPhi_;}
//		  	float closestCtfTrk_outerPt(){return eleStr_.closestCtfTrk_outerPt_;}
//		  	float closestCtfTrk_outerEta(){return eleStr_.closestCtfTrk_outerEta_;}
//		  	float closestCtfTrk_outerPhi(){return eleStr_.closestCtfTrk_outerPhi_;}

			// Various other variables ...
			float hOverE() const {return eleStr_.hOverE_;}
			float dEtaIn() const {return eleStr_.dEtaIn_;}
			float dPhiIn() const {return eleStr_.dPhiIn_;}
//			float dPhiOut(){return eleStr_.dPhiOut_;}
			float epIn(){return eleStr_.epIn_;}
			float epOut(){return eleStr_.epOut_;}
//			float fbrem(){return eleStr_.fbrem_;}
//			float bremFrac(){return eleStr_.bremFrac_;}
			float invEOverInvP(){return eleStr_.invEOverInvP_;}

			// Shower shape variables
//			float sigmaEtaEta(){return eleStr_.sigmaEtaEta_;}
//			float sigmaEtaEtaUnCorr(){return eleStr_.sigmaEtaEtaUnCorr_;}
			float sigmaIEtaIEta() const { return eleStr_.sigmaIEtaIEta_; }
//			float e1x5(){return eleStr_.e1x5_;}
//			float e2x5Max(){return eleStr_.e2x5Max_;}
//			float e5x5(){return eleStr_.e5x5_;}
			float e1x5Over5x5() const {return eleStr_.e1x5Over5x5_;}
			float e2x5MaxOver5x5() const {return eleStr_.e2x5MaxOver5x5_;}

			// Isolation variables ...
			float isolEm() const {return eleStr_.isolEm_;}
//			float isolHad(){return eleStr_.isolHad_;}
			float isolHadDepth1() const {return eleStr_.isolHadDepth1_;}
			float isolHadDepth2() const {return eleStr_.isolHadDepth2_;}
			float isolPtTrks() const {return eleStr_.isolPtTrks_;}
			float isolEmHadDepth1() const {return eleStr_.isolEmHadDepth1_;}

			// Alternative isolation variables from isoDeps module ...
			const double isol_isoDep_stdTrk() const { return eleStr_.isol_isoDep_stdTrk_; }
			const double isol_isoDep_stdEm()  const { return eleStr_.isol_isoDep_stdEm_; }
			const double isol_isoDep_stdHadD1() const { return eleStr_.isol_isoDep_stdHadD1_; }
			const double isol_isoDep_stdEmHadD1() const { return (eleStr_.isol_isoDep_stdEm_ + eleStr_.isol_isoDep_stdHadD1_); }
			const bool isolCut_isoDep_stdEmHadD1() const {
				return ApplyHEEPIsoCut_EmHad1( isol_isoDep_stdEmHadD1() ); }

			const double isol_isoDep_inrVetoModTrk() const { return eleStr_.isol_isoDep_inrVetoModTrk_; }
			const double isol_isoDep_inrVetoModEm()  const { return eleStr_.isol_isoDep_inrVetoModEm_; }
			const double isol_isoDep_inrVetoModHadD1() const { return eleStr_.isol_isoDep_inrVetoModHadD1_; }
			const double isol_isoDep_inrVetoModEmHadD1() const { return (eleStr_.isol_isoDep_inrVetoModEm_ + eleStr_.isol_isoDep_inrVetoModHadD1_); }
			const bool isolCut_isoDep_inrVetoModEmHadD1() const {
				return ApplyHEEPIsoCut_EmHad1(isol_isoDep_inrVetoModEmHadD1()); }

			double isol_inrVetoModTrk(const tsw::Event::InnerVetoSize) const;
			double isol_inrVetoModEm(const tsw::Event::InnerVetoSize) const;
			double isol_inrVetoModHadD1(const tsw::Event::InnerVetoSize) const;
			double isol_inrVetoModEmHadD1(const tsw::Event::InnerVetoSize) const;
			const bool isolCut_inrVetoModEmHadD1(const tsw::Event::InnerVetoSize) const;

			const tsw::ModEleIsoWithPhantom& modVetoIsoWithPhantomEle(const tsw::HEEPEle::ModIsoPhantomEleDrRange drRange) const { return eleStr_.isol_inrVetoModIsosWithPhantomEle_.at(drRange); }

			double modIsoTk_otherEleVetoForSelf()     const { return eleStr_.isol_inrVetoModTrk_otherEleAreaForSelf_; }    ///< 'Inner veto' mod track iso using 'other ele' inner veto area for this ele as well
			double modIsoEcal_otherEleVetoForSelf()   const { return eleStr_.isol_inrVetoModEcal_otherEleAreaForSelf_; }   ///< 'Inner veto' mod ECAL iso using 'other ele' inner veto area for this ele as well
			double modIsoHcalD1_otherEleVetoForSelf() const { return eleStr_.isol_inrVetoModHcalD1_otherEleAreaForSelf_; } ///< 'Inner veto' mod HCAL depth 1 iso using 'other ele' inner veto area for this ele as well

			unsigned int isol_nGenHadronsDr04() const { return eleStr_.isol_nGenHadronsDr04_; }
			double isol_ptSumGenHadronsDr04()   const { return eleStr_.isol_ptSumGenHadronsDr04_; }

			double isol_rhoCorrnEmH1(const tsw::EventHelper& evtHelper) const {
				return 0.28*(evtHelper.GetPURho()); }
			unsigned int numMissInnerHits()  const  { return eleStr_.numMissInnerHits_; }
			double dxy() const { return eleStr_.dxy_; }

			// Methods for applying HEEP cuts
			bool heepFidCut_v40() const;
			bool heepFidCut() const { return heepFidCut_v40();  }

			int heepIdNoIsoCutCode_v40() const;
			int heepIdNoIsoCutCode_v41() const;
			int heepIdNoIsoCutCode() const { return heepIdNoIsoCutCode_v41(); }
			bool heepIdNoIsoCut() const { return (heepIdNoIsoCutCode()==0); }

			int heepIdModIsoCutCode_v41_v00(const tsw::EventHelper& evtHelper) const;
			int heepIdModIsoCutCode(const tsw::EventHelper& evtHelper) const { return heepIdModIsoCutCode_v41_v00(evtHelper);}
			bool heepIdModIsoCut(const tsw::EventHelper& evtHelper) const { return (heepIdModIsoCutCode(evtHelper)==0); }

			int heepIdModIsoStdThrCutCode_v41(const tsw::EventHelper& evtHelper) const;
			int heepIdModIsoStdThrCutCode(const tsw::EventHelper& evtHelper) const { return heepIdModIsoStdThrCutCode_v41(evtHelper); }
			bool heepIdModIsoStdThrCut(const tsw::EventHelper& evtHelper) const { return (heepIdModIsoStdThrCutCode(evtHelper)==0); }

			int heepIdStdIsoCutCode_v41(const tsw::EventHelper& evtHelper) const;
			int heepIdStdIsoCutCode(const tsw::EventHelper& evtHelper) const { return heepIdStdIsoCutCode_v41(evtHelper); }
			bool heepIdStdIsoCut(const tsw::EventHelper& evtHelper) const { return (heepIdStdIsoCutCode(evtHelper)==0); }

			/// Fake rate pre-selection cut code method
			int fakeRatePreSelnCutCode_ichep2012() const;
			int fakeRatePreSelnCutCode() const { return fakeRatePreSelnCutCode_ichep2012(); }
			bool fakeRatePreSelnCut() const { return (fakeRatePreSelnCutCode()==0); }


			//////////////////////////////////////////////////////////////////////////////
			// METHODS FOR CALCULATING OLD MODIFIED ISOLATION CUTS

			//SC information (incl. individual recHits) ...
			float SC_eta() const { return eleStr_.SC_posn_eta_; }
			float SC_phi() const { return eleStr_.SC_posn_phi_; }
			float SC_rawEnergy() const { return eleStr_.SC_rawEnergy_; }
			unsigned int numStoredRecHits()  const {return eleStr_.SC_recHits_Et_.size();}
			float recHits_Et(int recHitIdx)  const {return eleStr_.SC_recHits_Et_.at(recHitIdx);}
			float recHits_eta(int recHitIdx) const {return eleStr_.SC_recHits_eta_.at(recHitIdx);}
			float recHits_phi(int recHitIdx) const {return eleStr_.SC_recHits_phi_.at(recHitIdx);}
			bool recHits_isFromEB(int recHitIdx) const {return eleStr_.SC_recHits_isFromEB_.at(recHitIdx);}
			float SC_totEnergyRecHits() const { return eleStr_.SC_totEnergyRecHits_; }
			float ratio_CaloToRecHitEnergy() const { return this->caloEnergy()/this->SC_totEnergyRecHits(); }
			unsigned int SC_totNumRecHits()  const { return eleStr_.SC_totNumRecHits_; }

			// Access to GSF track variables
			float gsfTrk_eta() const {return eleStr_.gsfTrk_eta_;}
			float gsfTrk_phi() const {return eleStr_.gsfTrk_phi_;}
			float gsfTrk_vz()  const {return eleStr_.gsfTrk_vz_;}

			// Accesss to information about CTF tracks that are within this electron's inner isol. cone.
			unsigned int numInnerIsoConeTrks() const {return eleStr_.innerIsoConeTrks_pt_.size();}
			float innerIsoConeTrks_pt(const int trkIdx) const
					{return eleStr_.innerIsoConeTrks_pt_.at(trkIdx);}
			float innerIsoConeTrks_eta(const int trkIdx) const
					{return eleStr_.innerIsoConeTrks_eta_.at(trkIdx);}
			float innerIsoConeTrks_phi(const int trkIdx) const
					{return eleStr_.innerIsoConeTrks_phi_.at(trkIdx);}
			float innerIsoConeTrks_vz(const int trkIdx) const
					{return eleStr_.innerIsoConeTrks_vz_.at(trkIdx);}

			/////////////////////////////////////////////////////

			// Method that modifies the standard tracker isolation value
			float modTrkIso(const tsw::HEEPEle* theOtherEle) const;

			//Methods used in applying the modified EmHad1 isolation HEEP cut ...
			float dR_SCs(tsw::HEEPEle* theOtherEle){
				// Calculate deltaEta ...
				float dEta = this->scEta() - theOtherEle->scEta();
				// Calculate deltaPhi ...
				double value_pi = 3.14159265;
				float dPhi = this->scPhi() - theOtherEle->scPhi();
				while(dPhi >= value_pi)
					dPhi -= 2.0*value_pi;
				while(dPhi < -value_pi)
					dPhi += 2.0*value_pi;

				// Now calculate deltaR ...
				float dR = sqrt( pow(dEta, 2.0) + pow(dPhi,2.0) );
				return dR;
			}
			float modEmHad1Iso_v0(tsw::HEEPEle* theOtherEle){
				bool coutDebugTxt = false;
				double isolEmHad1 = this->isolEmHadDepth1();
				// If deltaR for the di-electron SCs is < 0.3, then the centres of each electron's SC lies in the other ele's isolation cone, ...
				// and so remove the other electron's SC Et from that electron ...
				if(dR_SCs(theOtherEle)<0.3){
					isolEmHad1 -= theOtherEle->scEt();
					if(coutDebugTxt){std::cout << "                       << An isolEmHadDepth1 value has been modified!! >>" << std::endl;}
				}
				return isolEmHad1;
			}
			float modEmIso(const tsw::HEEPEle* theOtherEle) const ;
			float modHad1Iso(const tsw::HEEPEle* theOtherEle) const ;
			float modEmHad1Iso(const tsw::HEEPEle* anotherEle) const {
				return (this->modEmIso(anotherEle)+this->modHad1Iso(anotherEle));}
			float modEmHad1Iso_v1(tsw::HEEPEle* theOtherEle);


		public:
			// Static public members storing the index for each cut in HEEP cut codes
			static const int cutCode_Et_       = 0x0001;
			static const int cutCode_eta_      = 0x0002;
			static const int cutCode_ecalDrvn_ = 0x0004;
			static const int cutCode_dEtaIn_   = 0x0008;
			static const int cutCode_dPhiIn_   = 0x0010;
			static const int cutCode_hOverE_   = 0x0020;
			static const int cutCode_sieie_    = 0x0040;
			static const int cutCode_eXx5_     = 0x0080;
			static const int cutCode_isoEmH1_  = 0x0100;
			static const int cutCode_isoHad2_  = 0x0200;
			static const int cutCode_isoTrk_   = 0x0400;
			static const int cutCode_missHits_ = 0x0800;
			static const int cutCode_dxy_      = 0x1000;

			// Static public members storing the index for each cut in fake rate pre-selection cut codes
			static const int cutCode_frPre_Et_       = 0x0001;
			static const int cutCode_frPre_eta_      = 0x0002;
			static const int cutCode_frPre_ecalDrvn_ = 0x0004;
			static const int cutCode_frPre_sieie_    = 0x0008;
			static const int cutCode_frPre_hOverE_   = 0x0010;
			static const int cutCode_frPre_missHits_ = 0x0020;

	};

	HEEPEle::HEEPEle(tsw::EleStruct eleStruct){
		eleStr_ = eleStruct;
		if( (eleStr_.SC_recHits_Et_.size()!=eleStr_.SC_recHits_eta_.size()) || (eleStr_.SC_recHits_Et_.size()!=eleStr_.SC_recHits_phi_.size()) || (eleStr_.SC_recHits_Et_.size()!=eleStr_.SC_recHits_isFromEB_.size()) )
			std::cout << std::endl << "   *** WARNING (in tsw::HEEPEle CTOR): The SC_recHits_* vectors in the electron struct are NOT the same size!!! ***" << std::endl << std::endl;
	}
	HEEPEle::HEEPEle(){
		EleStruct eleStruct;
		SetDefaultValuesInEleStruct(eleStruct);
		eleStr_ = eleStruct;
	}
	HEEPEle::~HEEPEle(){
		//Placeholder ...
	}

	void HEEPEle::PrintOutVariables(){
		// Kinematic and geometric variables
		std::cout << "       et=" << et();
		std::cout << "; gsfEt=" << gsfEt();
		std::cout << "; scEt=" << scEt() << std::endl;
		std::cout << "       energy=" << energy();
//		std::cout << "; gsfEnergy=" << gsfEnergy();
		std::cout << "; caloEnergy=" << caloEnergy();
		std::cout << "; ecalEnergyError=" << ecalEnergyError() << std::endl;
		std::cout << "       eta=" << eta();
		std::cout << "; scEta=" << scEta();
//		std::cout << "; detEta=" << detEta();
		std::cout << std::endl;//		std::cout << "; detEtaAbs=" << detEtaAbs() << std::endl;
		std::cout << "       phi=" << phi();
		std::cout << "; scPhi=" << scPhi();
//		std::cout << "; detPhi=" << detPhi();
		std::cout <<  std::endl;//		std::cout << "; zVtx=" << zVtx() << std::endl;
		std::cout << "       p4,E=" << p4().E() << "; p4,p=(" << p4().Px() << ", " << p4().Py() << ", " << p4().Pz() << ")" << std::endl;
		std::cout << "       gsfP4,E=" << gsfP4().E() << "; gsfP4,p=(" << gsfP4().Px() << ", " << gsfP4().Py() << ", " << gsfP4().Pz() << ")" << std::endl;

		// Classification variables ...
		std::cout << "         -=-=-" << std::endl;
//		std::cout << "       classification=" << classification();
		std::cout << "       isEcalDriven=" << isEcalDriven();//		std::cout << "; isEcalDriven=" << isEcalDriven();
//		std::cout << "; isTrakerDriven=" << isTrackerDriven() << std::endl;
		std::cout << "       isEB=" << isEB();
		std::cout << "; isEE=" << isEE() <<std::endl;

		// Track methods ...
		std::cout << "         -=-=-" << std::endl;
		std::cout << "       charge=" << charge();
		std::cout << std::endl;//		std::cout << "; trkCharge=" << trkCharge() << std::endl;
//		std::cout << "       pVtx=" << pVtx();
//		std::cout << "; pCalo=" << pCalo();
		std::cout << "       ptVtx=" << ptVtx();
		std::cout << "; ptCalo=" << ptCalo() << std::endl;
		if(!closestCtfTrk_exists())
			std::cout << std::endl << "    ***** closestCtfTrackRef is NULL *****" << std::endl << std::endl;
		std::cout << "       closestCtfTrk... (pt,eta,phi)=(" << closestCtfTrk_pt() << ", " << closestCtfTrk_eta() << ", " << closestCtfTrk_phi() << ")" << std::endl;
//		std::cout << "          ... inner:    (pt,eta,phi)=(" << closestCtfTrk_innerPt() << ", " << closestCtfTrk_innerEta() << ", " << closestCtfTrk_innerPhi() << ")" << std::endl;
//		std::cout << "          ... outer:    (pt,eta,phi)=(" << closestCtfTrk_outerPt() << ", " << closestCtfTrk_outerEta() << ", " << closestCtfTrk_outerPhi() << ")" << std::endl;

		// Various other methods ...
		std::cout << "         -=-=-" << std::endl;
		std::cout << "       hOver=" << hOverE();
		std::cout << "; dEtaIn=" << dEtaIn();
		std::cout << "; dPhiIn=" << dPhiIn();
		std::cout << std::endl;//		std::cout << "; dPhiOut=" << dPhiOut() << std::endl;
		std::cout << "       epIn=" << epIn();
		std::cout << "; epOut=" << epOut();
//		std::cout << "; fbrem=" << fbrem();
//		std::cout << "; bremFrac=" << bremFrac();
		std::cout << "; invEOverInvP=" << invEOverInvP() << std::endl;

		// Shower shape variables ...
		std::cout << "         -=-=-" << std::endl;
//		std::cout << "       sigmaEtaEta=" << sigmaEtaEta();
//		std::cout << "; sigmaEtaEtaUnCorr=" << sigmaEtaEtaUnCorr();
		std::cout << "       sigmaIEtaIEta=" << sigmaIEtaIEta() << std::endl;//		std::cout << "; sigmaIEtaIEta=" << sigmaIEtaIEta() << std::endl;
//		std::cout << "       e1x5=" << e1x5();
//		std::cout << "; e2x5Max=" << e2x5Max();
//		std::cout << "; e5x5=" << e5x5();
		std::cout << "       e1x5Over5x5=" << e1x5Over5x5();//		std::cout << "; e1x5Over5x5=" << e1x5Over5x5();
		std::cout << "; e2x5MaxOver5x5=" << e2x5MaxOver5x5() << std::endl;

		// Isolation variables ...
		std::cout << "         -=-=-" << std::endl;
		std::cout << "       isolEm=" << isolEm();
//		std::cout << "; isolHad=" << isolHad();
		std::cout << "; isolHadDepth1=" << isolHadDepth1();
		std::cout << "; isolHadDepth2=" << isolHadDepth2() << std::endl;
		std::cout << "       isolPtTrks=" << isolPtTrks();
		std::cout << "; isolEmHadDepth1=" << isolEmHadDepth1() << std::endl;
		std::cout << "       numMissInnerHits = " << numMissInnerHits() << std::endl;

		// SC information ...
		std::cout << "         -=-=-" << std::endl;
		std::cout << "       SC: (eta,phi) = (" << SC_eta() << "," << SC_phi() << "); rawEnergy=" << SC_rawEnergy() << std::endl;
		std::cout << "       SC recHits (" << SC_totNumRecHits() << " of them, total energy = " << SC_totEnergyRecHits() << "; " << numStoredRecHits() << " of them stored in NTuple)" << std::endl;
		for(unsigned int recHitIdx=0; recHitIdx<numStoredRecHits(); recHitIdx++)
			std::cout << "          #" << recHitIdx << ": (Et, eta, phi, isFromEB)=(" << recHits_Et(recHitIdx) << ", " << recHits_eta(recHitIdx) << ", " << recHits_phi(recHitIdx) << "; " << recHits_isFromEB(recHitIdx) << ")" << std::endl;

		std::cout << "         -=-=-" << std::endl;
		std::cout << "       gsfTrk: (eta, phi) = (" << gsfTrk_eta() << ", " << gsfTrk_phi() << "); vz=" << gsfTrk_vz() << std::endl;

		std::cout << "         -=-=-" << std::endl;
		std::cout << "       innerIsoCone CTF tracks (There are " << numInnerIsoConeTrks() << " of them) ..." << std::endl;
		for(unsigned int trkIdx=0; trkIdx<numInnerIsoConeTrks(); trkIdx++)
			std::cout << "          #" << trkIdx << ": (pT, eta, phi) = (" << innerIsoConeTrks_pt(trkIdx) << ", " << innerIsoConeTrks_eta(trkIdx) << ", " << innerIsoConeTrks_phi(trkIdx) << "); vz=" << innerIsoConeTrks_vz(trkIdx) << std::endl;
	}

	bool HEEPEle::ApplySimpleCuts(){
		bool tmpCutsFlag = false;

		if( fabs(scEta()) < 1.442 ){ tmpCutsFlag=true; }
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) ){ tmpCutsFlag=true; }
		else{ tmpCutsFlag = false; }

		tmpCutsFlag = tmpCutsFlag && ( et()>12.0 );

		return tmpCutsFlag;
	}

	bool HEEPEle::isFiducialAndEcalDriven()
	{
		bool tmpCutsFlag = false;
		// First apply fiducial cuts (i.e. HEEP eta and Et requirements)
		if( isHEEPEB() )
			tmpCutsFlag = ( et()>35.0 );
		else if( isHEEPEE() )
			tmpCutsFlag = ( et()>40.0 );
		// Then make sure that electron is ECAL-driven
		tmpCutsFlag = tmpCutsFlag && isEcalDriven();

		return tmpCutsFlag;
	}

	const bool HEEPEle::isHEEPEB() const {
		return ( fabs(scEta()) < 1.442 );
	}
	const bool HEEPEle::isHEEPEE() const {
		return ( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) );
	}


	const bool HEEPEle::ApplyHEEPIsoCut_EmHad1(float emHad1IsolValue) const {
		bool tmpCutsFlag = false;
		if( isHEEPEB() ){
			tmpCutsFlag = ( emHad1IsolValue<(2.0+0.03*et()) );
		}
		else if( isHEEPEE() ){
			if( et()<50.0 ){
				tmpCutsFlag = ( emHad1IsolValue<2.5 );}
			else{
				tmpCutsFlag = ( emHad1IsolValue<(2.5+0.03*(et()-50.0)) );}
		}
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}


}// End: tsw namespace

///////////////////////////////////////////
// Modified isol value accessor methods ...

double tsw::HEEPEle::isol_inrVetoModTrk(const tsw::Event::InnerVetoSize vetoSize) const
{
	double trkIso = 9999.9;
	switch (vetoSize)
	{
		case tsw::Event::xSmallVeto:
			trkIso = eleStr_.isol_inrXSVetoModTrk_; break;
		case tsw::Event::smallVeto:
			trkIso = eleStr_.isol_inrSVetoModTrk_; break;
		case tsw::Event::mediumVeto:
			trkIso = eleStr_.isol_inrMVetoModTrk_; break;
		case tsw::Event::largeVeto:
			trkIso = eleStr_.isol_inrLVetoModTrk_; break;
		case tsw::Event::xLargeVeto:
			trkIso = eleStr_.isol_inrXLVetoModTrk_; break;
		default:
			std::cout << " *** ERROR : Unknown VetoSize enum value passed to tsw::HEEPEle::isol_inrVetoModTrk ! ***" << std::endl;
			break;
	}
	return trkIso;
}

double tsw::HEEPEle::isol_inrVetoModEm(const tsw::Event::InnerVetoSize vetoSize) const
{
	double ecalIso = 9999.9;
	switch (vetoSize)
	{
		case tsw::Event::xSmallVeto:
			ecalIso = eleStr_.isol_inrXSVetoModEm_; break;
		case tsw::Event::smallVeto:
			ecalIso = eleStr_.isol_inrSVetoModEm_; break;
		case tsw::Event::mediumVeto:
			ecalIso = eleStr_.isol_inrMVetoModEm_; break;
		case tsw::Event::largeVeto:
			ecalIso = eleStr_.isol_inrLVetoModEm_; break;
		case tsw::Event::xLargeVeto:
			ecalIso = eleStr_.isol_inrXLVetoModEm_; break;
		default:
			std::cout << " *** ERROR : Unknown VetoSize enum value passed to tsw::HEEPEle::isol_inrVetoModEm ! ***" << std::endl;
			break;
	}
	return ecalIso;
}

double tsw::HEEPEle::isol_inrVetoModHadD1(const tsw::Event::InnerVetoSize vetoSize) const
{
	double hcalD1Iso = 9999.9;
	switch (vetoSize)
	{
		case tsw::Event::xSmallVeto:
			hcalD1Iso = eleStr_.isol_inrXSVetoModHadD1_; break;
		case tsw::Event::smallVeto:
			hcalD1Iso = eleStr_.isol_inrSVetoModHadD1_; break;
		case tsw::Event::mediumVeto:
			hcalD1Iso = eleStr_.isol_inrMVetoModHadD1_; break;
		case tsw::Event::largeVeto:
			hcalD1Iso = eleStr_.isol_inrLVetoModHadD1_; break;
		case tsw::Event::xLargeVeto:
			hcalD1Iso = eleStr_.isol_inrXLVetoModHadD1_; break;
		default:
			std::cout << " *** ERROR : Unknown VetoSize enum value passed to tsw::HEEPEle::isol_inrVetoModHadD1 ! ***" << std::endl;
			break;
	}
	return hcalD1Iso;
}

double tsw::HEEPEle::isol_inrVetoModEmHadD1(const tsw::Event::InnerVetoSize vetoSize) const
{
	return ( isol_inrVetoModEm(vetoSize) + isol_inrVetoModHadD1(vetoSize) );
}
const bool tsw::HEEPEle::isolCut_inrVetoModEmHadD1(const tsw::Event::InnerVetoSize vetoSize) const
{
	return ApplyHEEPIsoCut_EmHad1( isol_inrVetoModEmHadD1(vetoSize) );
}


bool tsw::HEEPEle::heepFidCut_v40() const
{
	bool isFid = (isHEEPEB() || isHEEPEE());
	isFid = isFid && et()>35.0;
	isFid = isFid && isEcalDriven();
	return isFid;
}

int tsw::HEEPEle::heepIdNoIsoCutCode_v40() const
{
	const bool isInEB = (fabs(scEta()) < 1.442);
	const bool isInEE = ( fabs(scEta())>1.56 && fabs(scEta())<2.5 );
	const float thr_dEtaIn = isInEB ? 0.005 : 0.007 ;

	int cutCode = 0;

	// Fiducial cuts
	if( !(et()>35.0) )
		cutCode |= cutCode_Et_;
	if( !(isInEB || isInEE) )
		cutCode |= cutCode_eta_;

	if( !(isEcalDriven()) )
		cutCode |= cutCode_ecalDrvn_;

	// Track-SC matching cuts
	if( !(fabs(dEtaIn())<thr_dEtaIn) )
		cutCode |= cutCode_dEtaIn_;
	if( !(fabs(dPhiIn())<0.06) )
		cutCode |= cutCode_dPhiIn_;

	// Calo-based cuts
	if( !(hOverE()<0.05) )
		cutCode |= cutCode_hOverE_;
	if(  isInEE && !(sigmaIEtaIEta()<0.03)  )
		cutCode |= cutCode_sieie_;
	if(  isInEB && !((e2x5MaxOver5x5()>0.94) || (e1x5Over5x5()>0.83))  )
		cutCode |= cutCode_eXx5_;

	// Missing hits cuts
	if( numMissInnerHits()!=0 )
		cutCode |= cutCode_missHits_;

	return cutCode;
}


int tsw::HEEPEle::heepIdNoIsoCutCode_v41() const
{
	const bool isInEB = (fabs(scEta()) < 1.442);
	const float thr_dxy = isInEB ? 0.02 : 0.05 ;

	// Grab V4.0 cut code & reset missing hit cut bit
	int cutCode = heepIdNoIsoCutCode_v40();
	cutCode = cutCode & (~cutCode_missHits_);

	// Missing hits cut
	if( numMissInnerHits()>1 )
		cutCode |= cutCode_missHits_;
	if( !(abs(dxy())<thr_dxy) )
		cutCode |= cutCode_dxy_;

	return cutCode;
}


int tsw::HEEPEle::heepIdModIsoCutCode_v41_v00(const tsw::EventHelper& evtHelper) const
{
	const bool isInEB = (fabs(scEta()) < 1.442);

	// Grab 'heepNoIso' cut code
	int cutCode = heepIdNoIsoCutCode_v41();

	// ... add in iso bits with mod isolation cuts
	const float modEmH1isol = isol_inrVetoModEmHadD1(tsw::Event::mediumVeto) - isol_rhoCorrnEmH1(evtHelper);
	const float modTrkIsol  = isol_inrVetoModTrk(tsw::Event::xSmallVeto);
	if(  !( isInEB && modEmH1isol<(10.0+0.04*et()) )  )
		cutCode |= cutCode_isoEmH1_;
	if(  !( isInEB && modTrkIsol<(8.0+0.06*et()) )  )
		cutCode |= cutCode_isoTrk_;

	return cutCode;
}

int tsw::HEEPEle::heepIdModIsoStdThrCutCode_v41(const tsw::EventHelper& evtHelper) const
{
	const bool isInEB = (fabs(scEta()) < 1.442);
	const float thr_isoEmH1 = isInEB ? (2.0+0.03*et()) : ( et()<50.0 ? 2.5 : (2.5+0.03*(et()-50.0)) ) ;

	// Grab 'heepNoIso' cut code
	int cutCode = heepIdNoIsoCutCode_v41();

	// ... add in iso bits with mod isolation cuts
	const float modEmH1isol = isol_inrVetoModEmHadD1(tsw::Event::mediumVeto) - isol_rhoCorrnEmH1(evtHelper);
	const float modTrkIsol  = isol_inrVetoModTrk(tsw::Event::xSmallVeto);
	if(  !( modEmH1isol<thr_isoEmH1 )  )
		cutCode |= cutCode_isoEmH1_;
	if(  !( modTrkIsol<5.0 )  )
		cutCode |= cutCode_isoTrk_;

	return cutCode;

}

int tsw::HEEPEle::heepIdStdIsoCutCode_v41(const tsw::EventHelper& evtHelper) const
{
	const bool isInEB = (fabs(scEta()) < 1.442);
	const float thr_isoEmH1 = isInEB ? (2.0+0.03*et()) : ( et()<50.0 ? 2.5 : (2.5+0.03*(et()-50.0)) ) ;

	// Grab 'HEEPIdNoIso' cut code ...
	int cutCode = heepIdNoIsoCutCode_v41();

	// ... then add in iso bits using STANDARD iso cut
	const float stdEmH1isol = isolEmHadDepth1() - isol_rhoCorrnEmH1(evtHelper);
	if(  !( stdEmH1isol<thr_isoEmH1 )  )
		cutCode |= cutCode_isoEmH1_;
	if(  !( isolPtTrks()<5.0)  )
		cutCode |= cutCode_isoTrk_;

	return cutCode;
}


int tsw::HEEPEle::fakeRatePreSelnCutCode_ichep2012() const
{
	const bool isInEB = (fabs(scEta()) < 1.442);
	const bool isInEE = ( fabs(scEta())>1.56 && fabs(scEta())<2.5 );
	const float thr_Et  = isInEB ? 35.0 : 40.0 ;
	const float thr_sieie  = isInEB ? 0.013 : 0.034 ;
	const float thr_hOverE = isInEB ? 0.15 : 0.10 ;

	int cutCode = 0;

	// Fiducial cuts
	if( !(et()>thr_Et) )
		cutCode |= cutCode_frPre_Et_;
	if( !(isInEB || isInEE) )
		cutCode |= cutCode_frPre_eta_;

	if( !(isEcalDriven()) )
		cutCode |= cutCode_frPre_ecalDrvn_;

	if( !(sigmaIEtaIEta()<thr_sieie) )
		cutCode |= cutCode_frPre_sieie_;

	if( !(hOverE()<thr_hOverE) )
		cutCode |= cutCode_frPre_hOverE_;

	if( numMissInnerHits()!=0 )
		cutCode |= cutCode_frPre_missHits_;

	return cutCode;
}

//////////////////////////////////////////////////////////////////////////////////////
// METHODS FOR CALCULATING SC-BASED MODIFIED ISOLATION VALUES (written Nov/Dec 2011)
//////////////////////////////////////////////////////////////////////////////////////

float tsw::HEEPEle::modTrkIso(const tsw::HEEPEle* theOtherEle) const
{
	const bool coutDebugTxt = false;
	float isolPtTrkValue = this->isolPtTrks();
	const float thisEle_gsfEta = this->gsfTrk_eta();
	const float thisEle_gsfPhi = this->gsfTrk_phi();
	const float thisEle_gsfVz  = this->gsfTrk_vz();

	const float theOtherEle_gsfEta = theOtherEle->gsfTrk_eta();
	const float theOtherEle_gsfPhi = theOtherEle->gsfTrk_phi();

	// Now running over the CTF tracks that are within the other electron's inner isolation area
	for(unsigned int trkIdx=0; trkIdx<theOtherEle->numInnerIsoConeTrks(); trkIdx++){
		const float trk_pt  = theOtherEle->innerIsoConeTrks_pt(trkIdx);
		const float trk_eta = theOtherEle->innerIsoConeTrks_eta(trkIdx);
		const float trk_phi = theOtherEle->innerIsoConeTrks_phi(trkIdx);
		const float trk_vz  = theOtherEle->innerIsoConeTrks_vz(trkIdx);
		if(coutDebugTxt){std::cout << "                 innerIsoConeTrk: (pT, eta, phi) = (" << trk_pt << ", " << trk_eta << ", " << trk_phi << "); vz=" << trk_vz << std::endl;}

		// Double-check that this track is within the other electron's inner isolation area (skip track if its not)
		const float dEta_vsOtherEle = trk_eta - theOtherEle_gsfEta;
		const float dPhi_vsOtherEle = tsw::deltaPhi(trk_phi, theOtherEle_gsfPhi);
		const float dR_vsOtherEle   = sqrt(dEta_vsOtherEle*dEta_vsOtherEle + dPhi_vsOtherEle*dPhi_vsOtherEle);
		if(coutDebugTxt){std::cout << "                      rel to other ele: d(eta, phi) = (" << dEta_vsOtherEle << ", " << dPhi_vsOtherEle << "); dR=" << dR_vsOtherEle << std::endl;}
		if(fabs(dEta_vsOtherEle)>0.015 && dR_vsOtherEle>0.015)
			continue;
		if(dR_vsOtherEle>0.3)
			continue;

		// Check that this track was included in the tracker isolation sum for this electron (skip track if it wasn't)
		const float dEta_vsThisEle = trk_eta - thisEle_gsfEta;
		const float dPhi_vsThisEle = tsw::deltaPhi(trk_phi, thisEle_gsfPhi);
		const float dR_vsThisEle   = sqrt(dEta_vsThisEle*dEta_vsThisEle + dPhi_vsThisEle*dPhi_vsThisEle);
		const float dVz_vsThisEle  = trk_vz - thisEle_gsfVz;
		if(coutDebugTxt){std::cout << "                      rel to this ele: d(eta, phi) = (" << dEta_vsThisEle << ", " << dPhi_vsThisEle << "); dR=" << dR_vsThisEle << "; dVz=" << dVz_vsThisEle << std::endl;}
		if(trk_pt<0.7) // pT threshold
			continue;
		if(dVz_vsThisEle>0.2) //dz cut
			continue;
		if(fabs(dEta_vsThisEle)>0.015 && (dR_vsThisEle>0.015 && dR_vsThisEle<0.3) ){
			// A track that reaches here is within the other electron's inner isolation area, but was included in the isolation sum for this electron
			//  => Take it's pT away from this electron's tracker isolation variable
			isolPtTrkValue -= trk_pt;
			if(coutDebugTxt){std::cout << "                        *** Trk is in other ele's inner isol area, and in this ele's trk isol sum ***" << std::endl;
				std::cout << "                                     => isolPtTrkValue modified to " << isolPtTrkValue << std::endl;}
		}
	}

	return isolPtTrkValue;
}

float tsw::HEEPEle::modEmIso(const tsw::HEEPEle* theOtherEle) const
{
	bool coutDebugTxt = false;
	double isolEm = this->isolEm();

	float otherEle_CaloToRecHitEnergyRatio = theOtherEle->ratio_CaloToRecHitEnergy();
//	float otherEle_CaloToRecHitEnergyRatio = 1.005;
	if (otherEle_CaloToRecHitEnergyRatio<1.0)
		otherEle_CaloToRecHitEnergyRatio = 1.0;
	const float otherEle_leakedEt = ( theOtherEle->caloEnergy() - theOtherEle->SC_totEnergyRecHits() )*( theOtherEle->et()/theOtherEle->caloEnergy() );
	bool otherEle_isClose = false;
	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numStoredRecHits(); recHitIdx++){
		double recHit_Et = theOtherEle->recHits_Et(recHitIdx);
		double recHit_eta = theOtherEle->recHits_eta(recHitIdx);
		double recHit_phi = theOtherEle->recHits_phi(recHitIdx);
		bool recHit_isFromEB = theOtherEle->recHits_isFromEB(recHitIdx);
		TVector3 tmp3Vector; tmp3Vector.SetPtEtaPhi(recHit_Et,recHit_eta,recHit_phi);
		double recHit_energy = tmp3Vector.Mag();
		if(coutDebugTxt){std::cout << "                 recHit: (energy, Et, eta, phi) = (" << recHit_energy << ", " << recHit_Et << ", " << recHit_eta << ", " << recHit_phi << "), isFromEB=" << recHit_isFromEB << std::endl;}

		//Calculate delta eta and phi of this recHit relative to cluster
		double etaDiff = recHit_eta - scEta();
		double phiDiff = tsw::deltaPhi(recHit_phi, scPhi());
		double recHitSC_dR = sqrt(etaDiff*etaDiff+phiDiff*phiDiff);
		if(coutDebugTxt){std::cout << "                         (eta,phi)Diff = (" << etaDiff << ", " << phiDiff << "), recHitSC_dR=" << recHitSC_dR << std::endl;}

		// Skip to the next recHit if this one is not in electron's isolation cone ...
		double etaStripHalfWidth = 1.5; double isoCone_intRadius = 3.0;
		if( (  ( (recHitSC_dR>0.3 && recHitSC_dR<(0.3+0.035)) || (recHitSC_dR<0.0174*isoCone_intRadius && recHitSC_dR>0.0174*(isoCone_intRadius-2.0)) ) || (recHitSC_dR<0.3 && recHitSC_dR>0.0174*isoCone_intRadius && fabs(etaDiff)<0.0174*etaStripHalfWidth) ) ){
//			isolEm -= (otherEle_CaloToRecHitEnergyRatio-1.0)*recHit_Et; /// <--- MODFICATION HERE!!!
			otherEle_isClose = true;
		}
		if( recHitSC_dR>0.3 ) continue;
		if( fabs(scEta())<1.479 ){ //EB: Crystal width = 0.0174 in eta & phi
			if( fabs(etaDiff)<0.0174*etaStripHalfWidth ) continue;
			if(recHitSC_dR<isoCone_intRadius*0.0174) continue;
		}
		else{ //EE: Crystal width = 0.00864*fabs(sinh(recHit_eta))
			double xtalWidth_eta = 0.00864*fabs(sinh(recHit_eta));
			if( fabs(etaDiff)<etaStripHalfWidth*xtalWidth_eta ) continue;
			if( recHitSC_dR< isoCone_intRadius*xtalWidth_eta ) continue;
		}

		// However, if have got this far, and IF the recHit is above threshold Et/energy, then it contributed to the isolEm value ...
		// ... and so must take this recHit's Et away from the isolEm value of the electron.
		double thrEt_recHit = 0.0; double thrEnergy_recHit = 0.0;
		if(recHit_isFromEB){
			thrEt_recHit = 0.0; thrEnergy_recHit = 0.08;}
		else{
			thrEt_recHit = 0.1; thrEnergy_recHit = 0.0;}
		//if( (fabs(recHit_Et)>thrEt_recHit) && (fabs(recHit_energy)>thrEnergy_recHit) ){
		if( (recHit_Et>thrEt_recHit) && (recHit_energy>thrEnergy_recHit) ){
			isolEm -= recHit_Et; /// <--- MODFICATION HERE!!!
			otherEle_isClose = true;
			if(coutDebugTxt){std::cout << "                       << isolEm value has been modified!! (to " << isolEm << ") >>" << std::endl;}
		}
	}// End of for loop over recHits
//	if(otherEle_isClose && otherEle_leakedEt>0.0)
//		isolEm -= otherEle_leakedEt;

	return isolEm;
}

float tsw::HEEPEle::modHad1Iso(const tsw::HEEPEle* theOtherEle) const
{
	bool coutDebugTxt = false;
	double isolHad1 = this->isolHadDepth1();

	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numStoredRecHits(); recHitIdx++){
		double recHit_Et = theOtherEle->recHits_Et(recHitIdx);
		double recHit_eta = theOtherEle->recHits_eta(recHitIdx);
		double recHit_phi = theOtherEle->recHits_phi(recHitIdx);
		bool recHit_isFromEB = theOtherEle->recHits_isFromEB(recHitIdx);
		TVector3 tmp3Vector; tmp3Vector.SetPtEtaPhi(recHit_Et,recHit_eta,recHit_phi);
		double recHit_energy = tmp3Vector.Mag();
		if(coutDebugTxt){std::cout << "                 recHit: (energy, Et, eta, phi) = (" << recHit_energy << ", " << recHit_Et << ", " << recHit_eta << ", " << recHit_phi << "), isFromEB=" << recHit_isFromEB << std::endl;}

		//Calculate delta eta and phi of this recHit relative to cluster
		double etaDiff = recHit_eta - scEta();
		double phiDiff = tsw::deltaPhi(recHit_phi, scPhi());
		double recHitSC_dR = sqrt(etaDiff*etaDiff+phiDiff*phiDiff);
		if(coutDebugTxt){std::cout << "                         (eta,phi)Diff = (" << etaDiff << ", " << phiDiff << "), recHitSC_dR=" << recHitSC_dR << std::endl;}

		// Skip to the next recHit if this one is not in electron's isolation cone ...
		if( recHitSC_dR>0.3 )
			continue;
		if(recHitSC_dR<0.15)
			continue;

		// However, if have got this far, and IF the recHit is above threshold Et/energy, then it contributed to the isolHad1 value ...
		// ... and so must take this recHit's Et away from the isolHad1 value of the electron.
		double thrEt_recHit = 0.0; double thrEnergy_recHit = 0.0;
		if(recHit_isFromEB){
			thrEt_recHit = 0.0; thrEnergy_recHit = 0.08;}
		else{
			thrEt_recHit = 0.1; thrEnergy_recHit = 0.0;}
		//if( (fabs(recHit_Et)>thrEt_recHit) && (fabs(recHit_energy)>thrEnergy_recHit) ){
		if( (recHit_Et>thrEt_recHit) && (recHit_energy>thrEnergy_recHit) ){
			isolHad1 -= 1.00*(theOtherEle->hOverE())*recHit_Et; /// <--- MODFICATION HERE!!!
			if(coutDebugTxt){std::cout << "                       << isolHadDepth1 value has been modified!! (to " << isolHad1 << ") >>" << std::endl;}
		}
	}// End of for loop over recHits

	return isolHad1;
}

float tsw::HEEPEle::modEmHad1Iso_v1(tsw::HEEPEle* theOtherEle){
	bool coutDebugTxt = false;
	double isolEmHad1 = this->isolEmHadDepth1();

	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numStoredRecHits(); recHitIdx++){
		double recHit_Et = theOtherEle->recHits_Et(recHitIdx);
		double recHit_eta = theOtherEle->recHits_eta(recHitIdx);
		double recHit_phi = theOtherEle->recHits_phi(recHitIdx);
		bool recHit_isFromEB = theOtherEle->recHits_isFromEB(recHitIdx);
		TVector3 tmp3Vector; tmp3Vector.SetPtEtaPhi(recHit_Et,recHit_eta,recHit_phi);
		double recHit_energy = tmp3Vector.Mag();
		if(coutDebugTxt){std::cout << "                 recHit: (energy, Et, eta, phi) = (" << recHit_energy << ", " << recHit_Et << ", " << recHit_eta << ", " << recHit_phi << "), isFromEB=" << recHit_isFromEB << std::endl;}

		//Calculate delta eta and phi of this recHit relative to cluster
		double etaDiff = recHit_eta - scEta();
		double phiDiff = tsw::deltaPhi(recHit_phi, scPhi());
		double recHitSC_dR = sqrt(etaDiff*etaDiff+phiDiff*phiDiff);
		if(coutDebugTxt){std::cout << "                         (eta,phi)Diff = (" << etaDiff << ", " << phiDiff << "), recHitSC_dR=" << recHitSC_dR << std::endl;}

		// Skip to the next recHit if this one is not in electron's isolation cone ...
		double etaStripHalfWidth = 1.5; double isoCone_intRadius = 3.0;
		if( recHitSC_dR>0.3 ) continue;
		if( fabs(scEta())<1.479 ){ //EB: Crystal width = 0.0174 in eta & phi
			if( fabs(etaDiff)<0.0174*etaStripHalfWidth ) continue;
			if(recHitSC_dR<isoCone_intRadius*0.0174) continue;
		}
		else{ //EE: Crystal width = 0.00864*fabs(sinh(recHit_eta))
			double xtalWidth_eta = 0.00864*fabs(sinh(recHit_eta));
			if( fabs(etaDiff)<etaStripHalfWidth*xtalWidth_eta ) continue;
			if( recHitSC_dR< isoCone_intRadius*xtalWidth_eta ) continue;
		}

		// However, if have got this far, and IF the recHit is above threshold Et/energy, then it contributed to the isolEmHad1 value ...
		// ... and so must take this recHit's Et away from the isolEmHad1 value of the electron.
		double thrEt_recHit = 0.0; double thrEnergy_recHit = 0.0;
		if(recHit_isFromEB){
			thrEt_recHit = 0.0; thrEnergy_recHit = 0.08;}
		else{
			thrEt_recHit = 0.1; thrEnergy_recHit = 0.0;}
		//if( (fabs(recHit_Et)>thrEt_recHit) && (fabs(recHit_energy)>thrEnergy_recHit) ){
		if( (recHit_Et>thrEt_recHit) && (recHit_energy>thrEnergy_recHit) ){
			isolEmHad1 -= recHit_Et;
			if(coutDebugTxt){std::cout << "                       << isolEmHadDepth1 value has been modified!! (to " << isolEmHad1 << ") >>" << std::endl;}
		}
	}// End of for loop over recHits

	return isolEmHad1;
}

// *** HELPER FUNCTIONS *** //

namespace tsw{
	std::vector<unsigned int> IndicesOfElesPassingCuts(std::vector<tsw::HEEPEle> theElectrons, const std::vector<bool>& passCutsFlags, const int ecalRegionFlag);
	std::vector<tsw::HEEPEle> OrderHEEPElesByEt(std::vector<tsw::HEEPEle> vecOfEles);
} //end of namespace tsw


std::vector<unsigned int> tsw::IndicesOfElesPassingCuts(std::vector<tsw::HEEPEle> theElectrons, const std::vector<bool>& passCutsFlags, const int ecalRegionFlag)
{
	//ecalRegionFlag:
	//			=0 => Consider electrons in all regions
	//			=1 => Consider only EB electrons
	//			=2 => Consider only EE electrons
	std::vector<unsigned int> indicesVector;
	indicesVector.clear();

	if(ecalRegionFlag==0){
		for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
			if(passCutsFlags.at(idx))
				indicesVector.push_back(idx);
		}
	}
	else if(ecalRegionFlag==1){
		for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
			if( passCutsFlags.at(idx) && theElectrons.at(idx).isEB() )
				indicesVector.push_back(idx);
		}
	}
	else if(ecalRegionFlag==2){
		for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
			if( passCutsFlags.at(idx) && theElectrons.at(idx).isEE() )
				indicesVector.push_back(idx);
		}
	}

	return indicesVector;
}

	std::vector<tsw::HEEPEle> tsw::OrderHEEPElesByEt(std::vector<tsw::HEEPEle> vecOfEles)
	{
		// Declare variables used to temporarily store i'th and (i+1)'th electron
		tsw::HEEPEle ithEle;
		tsw::HEEPEle iPlusOnethEle;
		// Declare boolean used as flag for whether any electrons needed to be swapped in latest iteration ...
		bool swappedEles = true;

		bool beVerbose = false;
		if(beVerbose){
			std::cout << "   ***** Electrons being ordered:" << std::endl;
			std::cout << "   * Current ordering:" << std::endl;
			std::cout << "   *   et= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).et();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   *   eta= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).eta();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   *   phi= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).phi();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
		}

		// Can only compare two adjacent electrons' E_T's if there is more than one electron ...
		if(vecOfEles.size()>1){
			while(swappedEles==true){
				//Reset electron-swap boolean to false ...
				swappedEles = false;

				// Go through all elements of the vector, looking at the pairs of electrons that are adjacent under the current order ...
				for(unsigned int iEle=0; iEle < vecOfEles.size()-1 ; iEle++){
					ithEle = vecOfEles.at(iEle);
					iPlusOnethEle = vecOfEles.at(iEle+1);

					// If the i'th electron has less E_T than the (i+1)'th electron then swap these two electrons so that they are in descending E_T order ...
					if(ithEle.et()<iPlusOnethEle.et()){
						swappedEles = true;
						vecOfEles.at(iEle)   = iPlusOnethEle;
						vecOfEles.at(iEle+1) = ithEle;
					}
				}

				if(beVerbose){
					std::cout << "   * Reordering =>" << std::endl;
					std::cout << "   *   et= ";
					for(unsigned int i=0; i<vecOfEles.size(); i++){
						std::cout << vecOfEles.at(i).et();
						if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
					}
					std::cout << std::endl;
				}

			}
		}

		if(beVerbose){
			std::cout << "   * Final ordering:" << std::endl;
			std::cout << "   *   et= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).et();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   *   eta= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).eta();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   *   phi= ";
			for(unsigned int i=0; i<vecOfEles.size(); i++){
				std::cout << vecOfEles.at(i).phi();
				if(i!=( vecOfEles.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   *****Electrons should now be ordered in descending Et" << std::endl;
		}

		return vecOfEles;
	}




#endif
