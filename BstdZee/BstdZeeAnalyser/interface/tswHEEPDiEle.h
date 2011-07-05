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
			float et(){return totalP4().Pt();}
			float invMass(){return totalP4().M();}
			float openingAngle(){return eleA_.p4().Angle(eleB_.p4().Vect());};
			float deltaEta(){return eleB_.eta() - eleA_.eta();};
			float deltaPhi(){return eleB_.phi() - eleA_.phi();};
			float scDeltaEta(){return eleB_.scEta() - eleA_.scEta();};
			float scDeltaPhi(){return eleB_.scPhi() - eleA_.scPhi();};
			float deltaR(){return sqrt( pow(deltaEta(), 2.0) + pow(deltaPhi(),2.0) ); };

			bool isInZMassRange(){
				if( invMass()<120 && invMass()>60)
					return true;
				else
					return false;
			}

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

	void HEEPDiEle::PrintOutInfo(){
		std::cout << "  *** Di-electron information: " << std::endl;
		std::cout << "          p4,E=" << totalP4().E() << "; p4,p=(" << totalP4().Px() << ", " << totalP4().Py() << ", " <<totalP4().Pz() << ")" << std::endl;
		std::cout << "          et=" << et() << "; mass=" << invMass() << "; opening angle = " << openingAngle() << std::endl;
		std::cout << "          deltaEta=" << deltaEta() << "; deltaPhi=" << deltaPhi() << "; deltaR=" << deltaR() << std::endl;
		std::cout << "        EleA:" << std::endl;
		std::cout << "          p4,E=" << eleA_.p4().E() << "; p4,p=(" << eleA_.p4().Px() << ", " << eleA_.p4().Py() << ", " << eleA_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleA_.et() << "; energy=" << eleA_.energy() << "; eta=" << eleA_.eta() << "; phi=" << eleA_.phi() << std::endl;
		std::cout << "        EleB:" << std::endl;
		std::cout << "          p4,E=" << eleB_.p4().E() << "; p4,p=(" << eleB_.p4().Px() << ", " << eleB_.p4().Py() << ", " << eleB_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleB_.et() << "; energy=" << eleB_.energy() << "; eta=" << eleB_.eta() << "; phi=" << eleB_.phi() << std::endl;

	}

}// End of tsw namespace.

#endif
