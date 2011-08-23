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
				if( (closestCtfTrk_pt()<-0.1 && closestCtfTrk_innerPt()<-0.1) && closestCtfTrk_outerPt()<-0.1){
					//std::cout << std::endl << "    ***** closestCtfTrackRef is NULL *****" << std::endl << std::endl;
					return false;}
				else
					return true;
			}
		  	float closestCtfTrk_pt(){return eleStr_.closestCtfTrk_pt_;}
		  	float closestCtfTrk_eta(){return eleStr_.closestCtfTrk_eta_;}
		  	float closestCtfTrk_phi(){return eleStr_.closestCtfTrk_phi_;}
		  	float closestCtfTrk_innerPt(){return eleStr_.closestCtfTrk_innerPt_;}
		  	float closestCtfTrk_innerEta(){return eleStr_.closestCtfTrk_innerEta_;}
		  	float closestCtfTrk_innerPhi(){return eleStr_.closestCtfTrk_innerPhi_;}
		  	float closestCtfTrk_outerPt(){return eleStr_.closestCtfTrk_outerPt_;}
		  	float closestCtfTrk_outerEta(){return eleStr_.closestCtfTrk_outerEta_;}
		  	float closestCtfTrk_outerPhi(){return eleStr_.closestCtfTrk_outerPhi_;}

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
			/////////////////////////////////////////////////////

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
			float modEmHad1Iso(tsw::HEEPEle* theOtherEle){
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

		//private:
			//ClassDef(tsw::HEEPEle,1);
	};

	HEEPEle::HEEPEle(tsw::EleStruct eleStruct){
		eleStr_ = eleStruct;
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
		std::cout << "          ... inner:    (pt,eta,phi)=(" << closestCtfTrk_innerPt() << ", " << closestCtfTrk_innerEta() << ", " << closestCtfTrk_innerPhi() << ")" << std::endl;
		std::cout << "          ... outer:    (pt,eta,phi)=(" << closestCtfTrk_outerPt() << ", " << closestCtfTrk_outerEta() << ", " << closestCtfTrk_outerPhi() << ")" << std::endl;

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
	}

	bool HEEPEle::ApplySimpleCuts(){
		bool tmpCutsFlag = false;

		if( fabs(scEta()) < 1.442 ){ tmpCutsFlag=true; }
		else if( (fabs(scEta())>1.56) && (fabs(scEta())<2.5) ){ tmpCutsFlag=true; }
		else{ tmpCutsFlag = false; }

		tmpCutsFlag = tmpCutsFlag && ( et()>12.0 );

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
			// ---> N/A
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
