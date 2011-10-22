#ifndef tswHEEPEle_h
#define tswHEEPEle_h

//#include "NTupler/BstdZeeNTupler/interface/tswEleStruct.h"
#include "tswEleStruct.h"
#include "TObject.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include "BstdZeeFirst/Analyser/interface/tswUsefulFunctions.h"

#include <vector>
#include <functional>

namespace tsw{
	//class HEEPEle : public TObject {
	class HEEPEle{
		private:
			EleStruct eleStr_;
		public:
			HEEPEle();
			HEEPEle(tsw::EleStruct eleStruct);
			~HEEPEle();

			void PrintOutVariables();

			bool ApplySimpleCuts();
			bool isFiducialAndEcalDriven();
			bool isFidEcalDrAndPassesFRPre();
			bool isHEEPEB();
			bool isHEEPEE();
			bool ApplyHEEPCutsNoIso();
			bool ApplyHEEPIsoCut_Trk();
			bool ApplyHEEPIsoCut_EmHad1();
			bool ApplyHEEPIsoCut_EmHad1(float );
			bool ApplyHEEPIsoCut_Had2();
			bool ApplyIsoVarHEEPCuts();
			bool ApplyAllHEEPCuts();
			/////////////////////////////////////////////////
			//Methods for reading out all of the HEEP variables ...
			// Kinematic and geometric variables
			float et(){return eleStr_.et_;}
			float gsfEt(){return eleStr_.gsfEt_;}
			float scEt(){return eleStr_.scEt_;}
			float energy(){return eleStr_.energy_;}
//			float gsfEnergy(){return eleStr_.gsfEnergy_;}
			float caloEnergy(){return eleStr_.caloEnergy_;}
			float eta(){return eleStr_.eta_;}
			float scEta(){return eleStr_.scEta_;}
//			float detEta(){return eleStr_.detEta_;}
//			float detEtaAbs(){return eleStr_.detEtaAbs_;}
			float phi(){return eleStr_.phi_;}
			float scPhi(){return eleStr_.scPhi_;}
//			float detPhi(){return eleStr_.detPhi_;}
//			float zVtx(){return eleStr_.zVtx_;}
			TLorentzVector p4(){return ConvertToTLorentzVector( &(eleStr_.p4_) );}
			TLorentzVector gsfP4(){return ConvertToTLorentzVector( &(eleStr_.gsfP4_) );}

			// 'Classification'
//			int classification(){return eleStr_.classification_;}
			bool isEcalDriven(){return eleStr_.isEcalDriven_;}
//			bool isTrackerDriven(){return eleStr_.isTrackerDriven_;}
			bool isEB(){return eleStr_.isEB_;}
			bool isEE(){return eleStr_.isEE_;}

			// Track variables ...
			int charge(){return eleStr_.charge_;}
//			int trkCharge(){return eleStr_.trkCharge_;}
//			float pVtx(){return eleStr_.pVtx_;}
//			float pCalo(){return eleStr_.pCalo_;}
			float ptVtx(){return eleStr_.ptVtx_;}
			float ptCalo(){return eleStr_.ptCalo_;}
			bool  closestCtfTrk_exists(){
				if(closestCtfTrk_pt()<-0.1){
					//std::cout << std::endl << "    ***** closestCtfTrackRef is NULL *****" << std::endl << std::endl;
					return false;}
				else
					return true;
			}
		  	float closestCtfTrk_pt(){return eleStr_.closestCtfTrk_pt_;}
		  	float closestCtfTrk_eta(){return eleStr_.closestCtfTrk_eta_;}
		  	float closestCtfTrk_phi(){return eleStr_.closestCtfTrk_phi_;}
//		  	float closestCtfTrk_innerPt(){return eleStr_.closestCtfTrk_innerPt_;}
//		  	float closestCtfTrk_innerEta(){return eleStr_.closestCtfTrk_innerEta_;}
//		  	float closestCtfTrk_innerPhi(){return eleStr_.closestCtfTrk_innerPhi_;}
//		  	float closestCtfTrk_outerPt(){return eleStr_.closestCtfTrk_outerPt_;}
//		  	float closestCtfTrk_outerEta(){return eleStr_.closestCtfTrk_outerEta_;}
//		  	float closestCtfTrk_outerPhi(){return eleStr_.closestCtfTrk_outerPhi_;}

			// Various other variables ...
			float hOverE(){return eleStr_.hOverE_;}
			float dEtaIn(){return eleStr_.dEtaIn_;}
			float dPhiIn(){return eleStr_.dPhiIn_;}
//			float dPhiOut(){return eleStr_.dPhiOut_;}
			float epIn(){return eleStr_.epIn_;}
			float epOut(){return eleStr_.epOut_;}
//			float fbrem(){return eleStr_.fbrem_;}
//			float bremFrac(){return eleStr_.bremFrac_;}
			float invEOverInvP(){return eleStr_.invEOverInvP_;}

			// Shower shape variables
//			float sigmaEtaEta(){return eleStr_.sigmaEtaEta_;}
//			float sigmaEtaEtaUnCorr(){return eleStr_.sigmaEtaEtaUnCorr_;}
			float sigmaIEtaIEta(){return eleStr_.sigmaIEtaIEta_;}
//			float e1x5(){return eleStr_.e1x5_;}
//			float e2x5Max(){return eleStr_.e2x5Max_;}
//			float e5x5(){return eleStr_.e5x5_;}
			float e1x5Over5x5(){return eleStr_.e1x5Over5x5_;}
			float e2x5MaxOver5x5(){return eleStr_.e2x5MaxOver5x5_;}

			// Isolation variables ...
			float isolEm(){return eleStr_.isolEm_;}
//			float isolHad(){return eleStr_.isolHad_;}
			float isolHadDepth1(){return eleStr_.isolHadDepth1_;}
			float isolHadDepth2(){return eleStr_.isolHadDepth2_;}
			float isolPtTrks(){return eleStr_.isolPtTrks_;}
			float isolEmHadDepth1(){return eleStr_.isolEmHadDepth1_;}

			unsigned int numMissInnerHits(){ return eleStr_.numMissInnerHits_; }

			//SC recHits information ...
			unsigned int numRecHits(){return eleStr_.SC_recHits_Et_.size();}
			float recHits_Et(int recHitIdx){return eleStr_.SC_recHits_Et_.at(recHitIdx);}
			float recHits_eta(int recHitIdx){return eleStr_.SC_recHits_eta_.at(recHitIdx);}
			float recHits_phi(int recHitIdx){return eleStr_.SC_recHits_phi_.at(recHitIdx);}
			bool recHits_isFromEB(int recHitIdx){return eleStr_.SC_recHits_isFromEB_.at(recHitIdx);}

			// Access to GSF track variables
			float gsfTrk_eta(){return eleStr_.gsfTrk_eta_;}
			float gsfTrk_phi(){return eleStr_.gsfTrk_phi_;}
			float gsfTrk_vz(){return eleStr_.gsfTrk_vz_;}

			// Accesss to information about CTF tracks that are within this electron's inner isol. cone.
			unsigned int numInnerIsoConeTrks(){return eleStr_.innerIsoConeTrks_pt_.size();}
			float innerIsoConeTrks_pt(const int trkIdx)
					{return eleStr_.innerIsoConeTrks_pt_.at(trkIdx);}
			float innerIsoConeTrks_eta(const int trkIdx)
					{return eleStr_.innerIsoConeTrks_eta_.at(trkIdx);}
			float innerIsoConeTrks_phi(const int trkIdx)
					{return eleStr_.innerIsoConeTrks_phi_.at(trkIdx);}
			float innerIsoConeTrks_vz(const int trkIdx)
					{return eleStr_.innerIsoConeTrks_vz_.at(trkIdx);}

			/////////////////////////////////////////////////////

			// Method that modifies the standard tracker isolation value
			float modTrkIso(tsw::HEEPEle* theOtherEle);

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
			float modEmIso(tsw::HEEPEle* theOtherEle);
			float modHad1Iso(tsw::HEEPEle* theOtherEle);
			float modEmHad1Iso(tsw::HEEPEle* anotherEle){
				return (this->modEmIso(anotherEle)+this->modHad1Iso(anotherEle));}
			float modEmHad1Iso_v1(tsw::HEEPEle* theOtherEle);

		//private:
			//ClassDef(tsw::HEEPEle,1);
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
		std::cout << "; caloEnergy=" << caloEnergy() << std::endl;
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

		std::cout << "         -=-=-" << std::endl;
		std::cout << "       recHits (This ele is made up from " << numRecHits() << " of them)" << std::endl;
		for(unsigned int recHitIdx=0; recHitIdx<numRecHits(); recHitIdx++)
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

	bool HEEPEle::isFidEcalDrAndPassesFRPre()
	{
		// First apply the fiducial and ECAL-driven requirements ...
		bool tmpCutsFlag = isFiducialAndEcalDriven();
		// Then apply fake rate pre-selection from v6 of AN2011/159 (mentioned in trigger selection section, and detailed in table 10) ...
		if( isHEEPEB() ){
			tmpCutsFlag = tmpCutsFlag && ( sigmaIEtaIEta()<0.013 );
			tmpCutsFlag = tmpCutsFlag && ( hOverE()<0.15 );
		}
		else if ( isHEEPEE() ){
			tmpCutsFlag = tmpCutsFlag && ( sigmaIEtaIEta()<0.034 );
			tmpCutsFlag = tmpCutsFlag && ( hOverE()<0.10 );
		}
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}

	bool HEEPEle::isHEEPEB(){
		return ( fabs(scEta()) < 1.442 );
	}
	bool HEEPEle::isHEEPEE(){
		return ( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) );
	}

	bool HEEPEle::ApplyHEEPCutsNoIso(){
		bool tmpCutsFlag = false;

		if( fabs(scEta()) < 1.442 ){
			// E_T cut ...
			tmpCutsFlag = ( et()>35.0 );
			// isEcalDriven cut ...
			tmpCutsFlag = tmpCutsFlag && isEcalDriven();
			// dEtaIn cut ...
			tmpCutsFlag = tmpCutsFlag && ( fabs(dEtaIn())<0.005 );
			// dPhiIn cut ...
			tmpCutsFlag = tmpCutsFlag && ( fabs(dPhiIn())<0.09 );
			// H/E cut ...
			tmpCutsFlag = tmpCutsFlag && ( hOverE()<0.05 );
			// sigmaIEtaIEta cut ...
			//tmpCutsFlag = tmpCutsFlag && ( sigmaIEtaIEta()< 0.01);
			// E2x5/E5x5 cut ...
			tmpCutsFlag = tmpCutsFlag && ( (e2x5MaxOver5x5()>0.94) || (e1x5Over5x5()>0.83) );
		}
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) ){
			// E_T cut ...
			tmpCutsFlag = ( et()>40.0 );
			// isEcalDriven cut ...
			tmpCutsFlag = tmpCutsFlag && isEcalDriven();
			// dEtaIn cut ...
			tmpCutsFlag = tmpCutsFlag && ( fabs(dEtaIn())<0.007 );
			// dPhiIn cut ...
			tmpCutsFlag = tmpCutsFlag && ( fabs(dPhiIn())<0.09 );
			// H/E cut ...
			tmpCutsFlag = tmpCutsFlag && ( hOverE()<0.05 );
			// sigmaIEtaIEta cut ...
			tmpCutsFlag = tmpCutsFlag && ( sigmaIEtaIEta()< 0.03);
			// E2x5/E5x5 cut ...
			// ---> N/A
		}
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}

	bool HEEPEle::ApplyHEEPIsoCut_Trk(){
		bool tmpCutsFlag = false;
		if( fabs(scEta()) < 1.442 ){
			tmpCutsFlag = ( isolPtTrks()<7.5 );
		}
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) )
			tmpCutsFlag = ( isolPtTrks()<15.0 );
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}

	bool HEEPEle::ApplyHEEPIsoCut_EmHad1(){
		bool tmpCutsFlag = false;
		if( fabs(scEta()) < 1.442 ){
			tmpCutsFlag = ( isolEmHadDepth1()<(2.0+0.03*et()) );
		}
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) ){
			if( et()<50.0 ){
				tmpCutsFlag = ( isolEmHadDepth1()<2.5 );}
			else{
				tmpCutsFlag = ( isolEmHadDepth1()<(2.5+0.03*(et()-50.0)) );}
		}
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}
	bool HEEPEle::ApplyHEEPIsoCut_EmHad1(float emHad1IsolValue){
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

	bool HEEPEle::ApplyHEEPIsoCut_Had2(){
		bool tmpCutsFlag = false;

		if( fabs(scEta()) < 1.442 ){
			tmpCutsFlag = true; // (No had depth 2 isol cut)
		}
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) )
			tmpCutsFlag = ( isolHadDepth2()<0.5 );
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}

	bool HEEPEle::ApplyIsoVarHEEPCuts(){
		bool tmpCutsFlag = ( ApplyHEEPIsoCut_Trk() && ApplyHEEPIsoCut_EmHad1() );
		tmpCutsFlag = tmpCutsFlag && ApplyHEEPIsoCut_Had2();

		return tmpCutsFlag;
	}
	/*bool HEEPEle::ApplyIsoVarHEEPCuts(){
		bool tmpCutsFlag = false;

		if( fabs(scEta()) < 1.442 ){
			// EM + Had depth 1 isol cut ...
			tmpCutsFlag = ( isolEmHadDepth1()<(2.0+0.03*et()) );
			// (No had depth 2 isol cut)
			// Track isol cut ...
			tmpCutsFlag = tmpCutsFlag && ( isolPtTrks()<7.5 );
		}
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) ){
			// EM + Had depth 1 isol cut ...
			if( et()<50.0 ){
				tmpCutsFlag = ( isolEmHadDepth1()<2.5 );}
			else{
				tmpCutsFlag = ( isolEmHadDepth1()<(2.5+0.03*(et()-50.0)) );}
			// Had depth 2 isol cut ...
			tmpCutsFlag = tmpCutsFlag && ( isolHadDepth2()<0.5 );
			// Track isol cut ...
			tmpCutsFlag = tmpCutsFlag && ( isolPtTrks()<15.0 );
		}
		else
			tmpCutsFlag = false;

		return tmpCutsFlag;
	}*/
	bool HEEPEle::ApplyAllHEEPCuts(){
		return ( ApplyHEEPCutsNoIso() && ApplyIsoVarHEEPCuts() );
	}

	//ClassImp(tsw::HEEPEle)

	/*void HEEPEle::PrintOutVariables(){
		//Placeholder
	}*/

}

float tsw::HEEPEle::modTrkIso(tsw::HEEPEle* theOtherEle)
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

float tsw::HEEPEle::modEmIso(tsw::HEEPEle* theOtherEle){
	bool coutDebugTxt = false;
	double isolEm = this->isolEm();

	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numRecHits(); recHitIdx++){
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

		// However, if have got this far, and IF the recHit is above threshold Et/energy, then it contributed to the isolEm value ...
		// ... and so must take this recHit's Et away from the isolEm value of the electron.
		double thrEt_recHit = 0.0; double thrEnergy_recHit = 0.0;
		if(recHit_isFromEB){
			thrEt_recHit = 0.0; thrEnergy_recHit = 0.08;}
		else{
			thrEt_recHit = 0.1; thrEnergy_recHit = 0.0;}
		//if( (fabs(recHit_Et)>thrEt_recHit) && (fabs(recHit_energy)>thrEnergy_recHit) ){
		if( (recHit_Et>thrEt_recHit) && (recHit_energy>thrEnergy_recHit) ){
			isolEm -= recHit_Et;
			if(coutDebugTxt){std::cout << "                       << isolEm value has been modified!! (to " << isolEm << ") >>" << std::endl;}
		}
	}// End of for loop over recHits

	return isolEm;
}

float tsw::HEEPEle::modHad1Iso(tsw::HEEPEle* theOtherEle){
	bool coutDebugTxt = false;
	double isolHad1 = this->isolHadDepth1();

	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numRecHits(); recHitIdx++){
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
			isolHad1 -= (theOtherEle->hOverE())*recHit_Et;
			if(coutDebugTxt){std::cout << "                       << isolHadDepth1 value has been modified!! (to " << isolHad1 << ") >>" << std::endl;}
		}
	}// End of for loop over recHits

	return isolHad1;
}

float tsw::HEEPEle::modEmHad1Iso_v1(tsw::HEEPEle* theOtherEle){
	bool coutDebugTxt = false;
	double isolEmHad1 = this->isolEmHadDepth1();

	// Now running over the recHits from 'the other electron' ...
	for(unsigned int recHitIdx=0; recHitIdx<theOtherEle->numRecHits(); recHitIdx++){
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



namespace tsw{
	std::vector<tsw::HEEPEle> OrderHEEPElesByEt(std::vector<tsw::HEEPEle> vecOfEles){
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

} //end of namespace tsw

#endif
