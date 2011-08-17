#ifndef tswHEEPDiEle_h
#define tswHEEPDiEle_h

namespace tsw{
	class HEEPDiEle{
		public:
			// Constructors and destructors
			HEEPDiEle():
				eleA_(),
				eleB_()
			{}
			HEEPDiEle(const std::vector<tsw::HEEPEle>& etOrderedEles, const std::vector<bool> cutsPassFailFlags);
			~HEEPDiEle(){}

			//Methods that return properties of the di-electron pair ...
			TLorentzVector totalP4(){return (eleA_.p4() + eleB_.p4());}
			float pT(){return totalP4().Pt();}
			float invMass(){return totalP4().M();}
			float openingAngle(){return eleA_.p4().Angle(eleB_.p4().Vect());};
			float deltaEta(){return eleB_.eta() - eleA_.eta();};
			float deltaPhi(){return eleB_.phi() - eleA_.phi();};
			float scDeltaEta(){return eleB_.scEta() - eleA_.scEta();};
			float scDeltaPhi(){return eleB_.scPhi() - eleA_.scPhi();};
			float deltaR(){return sqrt( pow(deltaEta(), 2.0) + pow(deltaPhi(),2.0) ); };

			bool isEBEB(){
				if( eleA_.isEB() && eleB_.isEB() )
					return true;
				else
					return false;
			}
			bool isEBEE(){
				if( eleA_.isEB() && eleB_.isEE() )
					return true;
				else if( eleA_.isEE() && eleB_.isEB() )
					return true;
				else
					return false;
			}
			bool isEEEE(){
				if( eleA_.isEE() && eleB_.isEE() )
					return true;
				else
					return false;
			}

			bool isInZMassRange(){
				if( invMass()<120 && invMass()>60)
					return true;
				else
					return false;
			}

			//Methods that apply the modified track & EmHad1 isolation HEEP cuts ...
			bool ApplyDiEleTrkIsolCut();
//			bool ApplyDiEleEmHad1IsolCut();
//			bool ApplyDiEleHEEPIsolCuts(){
//				return ( ApplyDiEleTrkIsolCut() && ApplyDiEleEmHad1IsolCut() );
//			}

			//Methods that return the electrons ...
			tsw::HEEPEle eleA(){return eleA_;}
			tsw::HEEPEle eleB(){return eleB_;}

			//Method to write out general quantities of interest about the di-ele
			void PrintOutInfo();

		private:
			tsw::HEEPEle eleA_;
			tsw::HEEPEle eleB_;


	};

	HEEPDiEle::HEEPDiEle(const std::vector<tsw::HEEPEle>& etOrderedEles, const std::vector<bool> cutsPassFailFlags){
		// General declarations ...
		// Get vector of the indices of the electrons that passed the cuts ....
		std::vector<unsigned int> idxs_elesPassedCuts   = IndicesOfElesPassingCuts(etOrderedEles, cutsPassFailFlags, 0); //Flag of 0 => Include both EB and EE eles

		// Warn user if they are trying to form a di-electron pair from less than two cut-passing electrons ...
		if( NumPassingCuts(cutsPassFailFlags)<2 ){
			std::cout << std::endl;
			std::cout << "  ************************                              ************************ " << std::endl;
			std::cout << "  *** ERROR: You are trying to form a di-ele pair from < 2 cut-passing eles!!" << std::endl;
			std::cout << "  *** An out-of-range vector exception will now be thrown ..." << std::endl;
			std::cout << "  ************************                              ************************ " << std::endl;
			std::cout << std::endl;
		}

		// Assign the 1st electron (i.e. the highest Et ele) from this vector to be eleA_ ...
		eleA_ = etOrderedEles.at( idxs_elesPassedCuts.at(0) );
		// Assign the 2nd electron (i.e. the 2nd highest Et ele) from this vector to be eleB_ ...
		eleB_ = etOrderedEles.at( idxs_elesPassedCuts.at(1) );

		if(tsw::NumPassingCuts(cutsPassFailFlags)>2 && ( fabs(eleB_.eta()-eleA_.eta())>0.001 || fabs(eleB_.phi()-eleA_.phi())>0.001 ) ){
			eleB_ = etOrderedEles.at( idxs_elesPassedCuts.at(2) );
		}

		// Print out info about the di-electron ...
		//PrintOutInfo();
	}

	//-------------------------------------------------------------------***
	//Methods that apply the modified track & EmHad1 isolation HEEP cuts ...
	bool HEEPDiEle::ApplyDiEleTrkIsolCut(){
//		std::cout << "      **Applying the di-electron modified HEEP tracker isolation cut ...:" << std::endl;
		double eleA_isolPtTrks = eleA().isolPtTrks();
		bool eleA_cutFlag = false;
		double eleB_isolPtTrks = eleB().isolPtTrks();
		bool eleB_cutFlag = false;

//		std::cout << "      ** isolPtTrks..." << std::endl;
//		std::cout << "      **     Initial values...  eleA: " << eleA_isolPtTrks << ", eleB: " << eleB_isolPtTrks << std::endl;
		// If deltaR for the di-electron is < 0.3, then the two electrons lie in each other's isolation cone, ...
		// and so remove the other electron's track pT from that electron ...
		if(deltaR()<0.3){
			eleA_isolPtTrks -= eleB().ptCalo(); //eleA_isolPtTrks -= eleB().gsfP4().Pt();
			eleB_isolPtTrks -= eleA().ptCalo(); //eleB_isolPtTrks -= eleA().gsfP4().Pt();
//			std::cout << "                       << isolPtTrks values have been modified!! >>" << std::endl;
		}
//		std::cout << "      **     Modified values... eleA: " << eleA_isolPtTrks << ", eleB: " << eleB_isolPtTrks << std::endl;

		// Now, apply the HEEP track pT isolation cuts to the modified isolPtTrks values ...
		if( eleA().isHEEPEB() && (eleA_isolPtTrks<7.5) )
			eleA_cutFlag = true;
		else if( eleA().isHEEPEE() && (eleA_isolPtTrks<15.0) )
			eleA_cutFlag = true;
		else
			eleA_cutFlag = false;
//		std::cout << "      ** Applying HEEP cut to electron A => " << eleA_cutFlag << std::endl;

		if( eleB().isHEEPEB() && (eleB_isolPtTrks<7.5) )
			eleB_cutFlag = true;
		else if( eleB().isHEEPEE() && (eleB_isolPtTrks<15.0) )
			eleB_cutFlag = true;
		else
			eleB_cutFlag = false;
//		std::cout << "      ** Applying HEEP cut to electron B => " << eleB_cutFlag << std::endl;

//		std::cout << "      ** Value returned = " << (eleA_cutFlag && eleB_cutFlag) << std::endl;
		return (eleA_cutFlag && eleB_cutFlag);
	}
	//bool HEEPDiEle::ApplyDiEleEmHad1IsolCut();

	void HEEPDiEle::PrintOutInfo(){
		std::cout << "  *** Di-electron information: " << std::endl;
		std::cout << "          p4,E=" << totalP4().E() << "; p4,p=(" << totalP4().Px() << ", " << totalP4().Py() << ", " <<totalP4().Pz() << ")" << std::endl;
		std::cout << "          pT=" << pT() << "; mass=" << invMass() << "; opening angle = " << openingAngle() << std::endl;
		std::cout << "          deltaEta=" << deltaEta() << "; deltaPhi=" << deltaPhi() << "; deltaR=" << deltaR() << std::endl;
		std::cout << "        EleA:" << std::endl;
		std::cout << "          p4,E=" << eleA_.p4().E() << "; p4,p=(" << eleA_.p4().Px() << ", " << eleA_.p4().Py() << ", " << eleA_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleA_.et() << "; energy=" << eleA_.energy() << "; eta=" << eleA_.eta() << "; phi=" << eleA_.phi() << std::endl;
		std::cout << "          ptVtx=" << eleA_.ptVtx() << "; ptCalo=" << eleA_.ptCalo() << std::endl;
		std::cout << "        EleB:" << std::endl;
		std::cout << "          p4,E=" << eleB_.p4().E() << "; p4,p=(" << eleB_.p4().Px() << ", " << eleB_.p4().Py() << ", " << eleB_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleB_.et() << "; energy=" << eleB_.energy() << "; eta=" << eleB_.eta() << "; phi=" << eleB_.phi() << std::endl;
		std::cout << "          ptVtx=" << eleB_.ptVtx() << "; ptCalo=" << eleB_.ptCalo() << std::endl;

	}

}// End of tsw namespace.

#endif
