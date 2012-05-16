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

			HEEPDiEle(tsw::HEEPEle ele1, tsw::HEEPEle ele2):
				cache_eleA_modEmHad1Iso_(9999.9),
				cache_eleB_modEmHad1Iso_(9999.9),
				cache_eleA_modHad1Iso_(9999.9),
				cache_eleB_modHad1Iso_(9999.9),
				eleA_modEmHad1Iso_isCached_(false),
				eleB_modEmHad1Iso_isCached_(false),
				eleA_modHad1Iso_isCached_(false),
				eleB_modHad1Iso_isCached_(false)
			{
				// Assign the higher Et electron to be eleA_, and assign the other ele to be eleB_ ...
				if(ele1.et()>ele2.et()){
					eleA_ = ele1;
					eleB_ = ele2;
				}
				else{
					eleA_ = ele2;
					eleB_ = ele1;
				}
			}

			HEEPDiEle(const std::vector<tsw::HEEPEle>& etOrderedEles, const std::vector<bool> cutsPassFailFlags);
			~HEEPDiEle(){}

			//Methods that return properties of the di-electron pair ...
			TLorentzVector totalP4()  const  {return (eleA_.p4() + eleB_.p4());}
			TLorentzVector p4()  const  { return totalP4(); }
			float pT()           const  {return totalP4().Pt();}
			float invMass()      const  {return totalP4().M();}
			float openingAngle() const  {return eleA_.p4().Angle(eleB_.p4().Vect());}
			float deltaEta()     const  {return eleB_.eta() - eleA_.eta();}
			float deltaPhi()     const  {
				double value_pi = 3.14159265;
				float dPhi = eleB_.phi() - eleA_.phi();
				while(dPhi >= value_pi)
					dPhi -= 2.0*value_pi;
				while(dPhi < -value_pi)
					dPhi += 2.0*value_pi;
				return dPhi;
			}
			float scDeltaEta() const  {return eleB_.scEta() - eleA_.scEta();}
			float scDeltaPhi() const  {return eleB_.scPhi() - eleA_.scPhi();}
			float deltaR()     const  {return sqrt( pow(deltaEta(), 2.0) + pow(deltaPhi(),2.0) ); }

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
			bool isAboveZMassRange(){
				if( invMass()>120.0 )
					return true;
				else
					return false;
			}

			//Methods for applying the modified track isolation HEEP cut - w/ modification occuring using innerIsoConeTrks information
			float eleA_modTrkIso(){ return eleA_.modTrkIso(&eleB_); }
			float eleB_modTrkIso(){ return eleB_.modTrkIso(&eleA_); }
			bool eleA_modTrkIsolCut();
			bool eleB_modTrkIsolCut();
			bool ApplyDiEleTrkIsolCut();

			// The old (simpler) methods for applying the modified track isolation HEEP cut.
			float dR_gsfEleToCtfTrk(tsw::HEEPEle*, float, float);
			float dR_gsfAToCtfB();
			float dR_gsfBToCtfA();
			float eleA_modTrkIso_v0();
			float eleB_modTrkIso_v0();
			bool ApplyDiEleTrkIsolCut_v0();

		private:
			// Member variables for the caching of the modified electron isolation values ...
			float cache_eleA_modEmHad1Iso_;
			float cache_eleB_modEmHad1Iso_;
			float cache_eleA_modHad1Iso_;
			float cache_eleB_modHad1Iso_;
			bool eleA_modEmHad1Iso_isCached_;
			bool eleB_modEmHad1Iso_isCached_;
			bool eleA_modHad1Iso_isCached_;
			bool eleB_modHad1Iso_isCached_;

		public:
			//Methods for applying the modified EmHad1 isolation HEEP cut ...
			float eleA_modEmHad1Iso(){
				if(!eleA_modEmHad1Iso_isCached_){
					cache_eleA_modEmHad1Iso_ = eleA_.modEmHad1Iso(&eleB_);
					eleA_modEmHad1Iso_isCached_ = true;
				}
				return cache_eleA_modEmHad1Iso_;
			}
			float eleB_modEmHad1Iso(){
				if(!eleB_modEmHad1Iso_isCached_){
					cache_eleB_modEmHad1Iso_ = eleB_.modEmHad1Iso(&eleA_);
					eleB_modEmHad1Iso_isCached_ = true;
				}
				return cache_eleB_modEmHad1Iso_;
			}
			float eleA_modHad1Iso(){
				if(!eleA_modHad1Iso_isCached_){
					cache_eleA_modHad1Iso_ = eleA_.modHad1Iso(&eleB_);
					eleA_modHad1Iso_isCached_ = true;
				}
				return cache_eleA_modHad1Iso_;
			}
			float eleB_modHad1Iso(){
				if(!eleB_modHad1Iso_isCached_){
					cache_eleB_modHad1Iso_ = eleB_.modHad1Iso(&eleA_);
					eleB_modHad1Iso_isCached_ = true;
				}
				return cache_eleB_modHad1Iso_;
			}
			bool eleA_modEmHad1IsoCut();
			bool eleB_modEmHad1IsoCut();
			bool ApplyDiEleEmHad1IsolCut();
//			bool ApplyDiEleHEEPIsolCuts(){
//				return ( ApplyDiEleTrkIsolCut() && ApplyDiEleEmHad1IsolCut() );
//			}

			//Methods that return the electrons ...
			tsw::HEEPEle eleA()  const  { return eleA_; }
			tsw::HEEPEle eleB()  const  { return eleB_; }

			//Method to write out general quantities of interest about the di-ele
			void PrintOutInfo();

		private:
			tsw::HEEPEle eleA_;
			tsw::HEEPEle eleB_;
	};

	HEEPDiEle::HEEPDiEle(const std::vector<tsw::HEEPEle>& etOrderedEles, const std::vector<bool> cutsPassFailFlags):
		cache_eleA_modEmHad1Iso_(9999.9),
		cache_eleB_modEmHad1Iso_(9999.9),
		cache_eleA_modHad1Iso_(9999.9),
		cache_eleB_modHad1Iso_(9999.9),
		eleA_modEmHad1Iso_isCached_(false),
		eleB_modEmHad1Iso_isCached_(false),
		eleA_modHad1Iso_isCached_(false),
		eleB_modHad1Iso_isCached_(false)
	{
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
	//Methods for applying the modified track isolation HEEP cut ...
	bool HEEPDiEle::eleA_modTrkIsolCut()
	{
		float eleA_isolPtTrks = eleA_modTrkIso();
		bool eleA_cutFlag = false;

		if( eleA().isHEEPEB() && (eleA_isolPtTrks<5.0) )
			eleA_cutFlag = true;
		else if( eleA().isHEEPEE() && (eleA_isolPtTrks<5.0) )
			eleA_cutFlag = true;
		else
			eleA_cutFlag = false;

		return eleA_cutFlag;
	}
	bool HEEPDiEle::eleB_modTrkIsolCut()
	{
		float eleB_isolPtTrks = eleB_modTrkIso();
		bool eleB_cutFlag = false;

		if( eleB().isHEEPEB() && (eleB_isolPtTrks<5.0) )
			eleB_cutFlag = true;
		else if( eleB().isHEEPEE() && (eleB_isolPtTrks<5.0) )
			eleB_cutFlag = true;
		else
			eleB_cutFlag = false;

		return eleB_cutFlag;
	}

	bool HEEPDiEle::ApplyDiEleTrkIsolCut()
	{
		const bool coutDebugTxt = false;
		if(coutDebugTxt){
			std::cout << "      **Applying the di-electron modified HEEP tracker isolation cut ...:" << std::endl;
			std::cout << "      ** isolPtTrks..." << std::endl;
			std::cout << "      **     Standard values...  eleA: " << eleA().isolPtTrks() << ", eleB: " << eleB().isolPtTrks() << std::endl;
			std::cout << "      **     Ele A calc'n ..." << std::endl;
		}

		float eleA_isolPtTrks = eleA_modTrkIso();
		if(coutDebugTxt){std::cout << "      **     Ele B calc'n ..." << std::endl;}
		float eleB_isolPtTrks = eleB_modTrkIso();
		if(coutDebugTxt){std::cout << "      **     Modified values... eleA: " << eleA_isolPtTrks << ", eleB: " << eleB_isolPtTrks << std::endl;}

		// Now, apply the HEEP track pT isolation cuts to the modified isolPtTrks values ...
		bool eleA_cutFlag = false;
		bool eleB_cutFlag = false;

		if( eleA().isHEEPEB() && (eleA_isolPtTrks<5.0) )
			eleA_cutFlag = true;
		else if( eleA().isHEEPEE() && (eleA_isolPtTrks<5.0) )
			eleA_cutFlag = true;
		else
			eleA_cutFlag = false;
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron A => " << eleA_cutFlag << std::endl;}

		if( eleB().isHEEPEB() && (eleB_isolPtTrks<5.0) )
			eleB_cutFlag = true;
		else if( eleB().isHEEPEE() && (eleB_isolPtTrks<5.0) )
			eleB_cutFlag = true;
		else
			eleB_cutFlag = false;
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron B => " << eleB_cutFlag << std::endl;}

		if(coutDebugTxt){std::cout << "      ** Value returned = " << (eleA_cutFlag && eleB_cutFlag) << std::endl;}
		return (eleA_cutFlag && eleB_cutFlag);
	}

	//-------------------------------------------------------------------***
	// OLD Methods for applying the modified track isolation HEEP cut ...
	float HEEPDiEle::dR_gsfEleToCtfTrk(tsw::HEEPEle* theGsfEle, float ctfTrk_eta, float ctfTrk_phi){
		// Calculate deltaEta ...
		float dEta = theGsfEle->eta() - ctfTrk_eta;
		// Calculate deltaPhi ...
		double value_pi = 3.14159265;
		float dPhi = theGsfEle->phi() - ctfTrk_phi;
		while(dPhi >= value_pi)
			dPhi -= 2.0*value_pi;
		while(dPhi < -value_pi)
			dPhi += 2.0*value_pi;

		// Now calculate deltaR ...
		float dR = sqrt( pow(dEta, 2.0) + pow(dPhi,2.0) );
		return dR;
	}
	float HEEPDiEle::dR_gsfAToCtfB(){
		// Calculating dR using the gsfTrk of eleA and the CTF track of eleB ...
		return dR_gsfEleToCtfTrk(&eleA_, eleB_.closestCtfTrk_eta(), eleB_.closestCtfTrk_phi());
	}
	float HEEPDiEle::dR_gsfBToCtfA(){
		// Calculating dR using the gsfTrk of eleB and the CTF track of eleA ...
		return dR_gsfEleToCtfTrk(&eleB_, eleA_.closestCtfTrk_eta(), eleA_.closestCtfTrk_phi());
	}



	float HEEPDiEle::eleA_modTrkIso_v0()
	{
		bool coutDebugTxt = false;
		double eleA_isolPtTrks = eleA().isolPtTrks();
		// If deltaR for the di-electron is < 0.3, then the two electrons lie in each other's isolation cone, ...
		// and so remove the other electron's (i.e. eleB's) track pT from that electron ...
		if(eleB().closestCtfTrk_exists() && dR_gsfAToCtfB()<0.3){
			eleA_isolPtTrks -= eleB().closestCtfTrk_pt(); //eleA_isolPtTrks -= eleB().gsfP4().Pt();
			if(coutDebugTxt){std::cout << "                       << isolPtTrks values have been modified!! >>" << std::endl;}
		}
		return eleA_isolPtTrks;
	}
	float HEEPDiEle::eleB_modTrkIso_v0()
	{
		bool coutDebugTxt = false;
		double eleB_isolPtTrks = eleB().isolPtTrks();
		// If deltaR for the di-electron is < 0.3, then the two electrons lie in each other's isolation cone, ...
		// and so remove the other electron's (i.e. eleA's) track pT from that electron ...
		if(eleA().closestCtfTrk_exists() && dR_gsfBToCtfA()<0.3){
			eleB_isolPtTrks -= eleA().closestCtfTrk_pt(); //eleA_isolPtTrks -= eleB().gsfP4().Pt();
			if(coutDebugTxt){std::cout << "                       << isolPtTrks values have been modified!! >>" << std::endl;}
		}
		return eleB_isolPtTrks;

	}

	bool HEEPDiEle::ApplyDiEleTrkIsolCut_v0()
	{
		bool coutDebugTxt = false;
		if(coutDebugTxt){std::cout << "      **Applying the di-electron modified HEEP tracker isolation cut ...:" << std::endl;}
		double eleA_isolPtTrks = eleA().isolPtTrks();
		bool eleA_cutFlag = false;
		double eleB_isolPtTrks = eleB().isolPtTrks();
		bool eleB_cutFlag = false;

		if(coutDebugTxt){
			std::cout << "      ** isolPtTrks..." << std::endl;
			std::cout << "      **     Initial values...  eleA: " << eleA_isolPtTrks << ", eleB: " << eleB_isolPtTrks << std::endl;}
		// If deltaR for the di-electron is < 0.3, then the two electrons lie in each other's isolation cone, ...
		// and so remove the other electron's track pT from that electron ...
		eleA_isolPtTrks = eleA_modTrkIso_v0(); //eleA_isolPtTrks -= eleB().gsfP4().Pt();
		eleB_isolPtTrks = eleB_modTrkIso_v0(); //eleB_isolPtTrks -= eleA().gsfP4().Pt();
		if(coutDebugTxt){std::cout << "      **     Modified values... eleA: " << eleA_isolPtTrks << ", eleB: " << eleB_isolPtTrks << std::endl;}

		// Now, apply the HEEP track pT isolation cuts to the modified isolPtTrks values ...
		if( eleA().isHEEPEB() && (eleA_isolPtTrks<7.5) )
			eleA_cutFlag = true;
		else if( eleA().isHEEPEE() && (eleA_isolPtTrks<15.0) )
			eleA_cutFlag = true;
		else
			eleA_cutFlag = false;
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron A => " << eleA_cutFlag << std::endl;}

		if( eleB().isHEEPEB() && (eleB_isolPtTrks<7.5) )
			eleB_cutFlag = true;
		else if( eleB().isHEEPEE() && (eleB_isolPtTrks<15.0) )
			eleB_cutFlag = true;
		else
			eleB_cutFlag = false;
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron B => " << eleB_cutFlag << std::endl;}

		if(coutDebugTxt){std::cout << "      ** Value returned = " << (eleA_cutFlag && eleB_cutFlag) << std::endl;}
		return (eleA_cutFlag && eleB_cutFlag);
	}

	//-------------------------------------------------------------------***
	//Method for applying the modified EmHad1 isolation HEEP cut ...
	bool HEEPDiEle::eleA_modEmHad1IsoCut()
	{
		double eleA_isolEmHad1 = eleA_modEmHad1Iso();
		bool eleA_cutFlag = eleA_.ApplyHEEPIsoCut_EmHad1( eleA_isolEmHad1 );

		return eleA_cutFlag;
	}
	bool HEEPDiEle::eleB_modEmHad1IsoCut()
	{
		double eleB_isolEmHad1 = eleB_modEmHad1Iso();
		bool eleB_cutFlag = eleB_.ApplyHEEPIsoCut_EmHad1( eleB_isolEmHad1 );

		return eleB_cutFlag;
	}

	bool HEEPDiEle::ApplyDiEleEmHad1IsolCut(){
		bool coutDebugTxt = false;
		if(coutDebugTxt){std::cout << "      **Applying the di-electron modified HEEP EmHad1 isolation cut ...:" << std::endl;}
		double eleA_isolEmHad1 = eleA().isolEmHadDepth1();
		bool eleA_cutFlag = false;
		double eleB_isolEmHad1 = eleB().isolEmHadDepth1();
		bool eleB_cutFlag = false;

		if(coutDebugTxt){
			std::cout << "      ** isolEmHad1:" << std::endl;
			std::cout << "      **     Initial values...  eleA: " << eleA_isolEmHad1 << ", eleB: " << eleB_isolEmHad1 << std::endl;}

		// Getting the modified isolation values ...
		if(coutDebugTxt){std::cout << "             *** eleA calcn ... ***" << std::endl;}
		eleA_isolEmHad1 = eleA_modEmHad1Iso();
		if(coutDebugTxt){std::cout << "             *** eleB calcn ... ***" << std::endl;}
		eleB_isolEmHad1 = eleB_modEmHad1Iso();
		if(coutDebugTxt){std::cout << "      **     Modified values... eleA: " << eleA_isolEmHad1 << ", eleB: " << eleB_isolEmHad1 << std::endl;}

		// Now, apply the HEEP EmHad1 isolation cuts to the modified isolEmHad1 values ...
		eleA_cutFlag = eleA_.ApplyHEEPIsoCut_EmHad1( eleA_isolEmHad1 );
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron A => " << eleA_cutFlag << std::endl;}

		eleB_cutFlag = eleB_.ApplyHEEPIsoCut_EmHad1( eleB_isolEmHad1 );
		if(coutDebugTxt){std::cout << "      ** Applying HEEP cut to electron B => " << eleB_cutFlag << std::endl;}

		if(coutDebugTxt){std::cout << "      ** Value returned = " << (eleA_cutFlag && eleB_cutFlag) << std::endl;}
		return (eleA_cutFlag && eleB_cutFlag);

	}

	void HEEPDiEle::PrintOutInfo(){
		std::cout << "  *** Di-electron information: " << std::endl;
		std::cout << "          p4,E=" << totalP4().E() << "; p4,p=(" << totalP4().Px() << ", " << totalP4().Py() << ", " <<totalP4().Pz() << ")" << std::endl;
		std::cout << "          pT=" << pT() << "; mass=" << invMass() << "; opening angle = " << openingAngle() << std::endl;
		std::cout << "          deltaEta=" << deltaEta() << "; deltaPhi=" << deltaPhi() << "; deltaR=" << deltaR() << std::endl;
		std::cout << "          dR_gsfAToCtfB=" << dR_gsfAToCtfB() << "; dR_gsfBToCtfA=" << dR_gsfBToCtfA() << "; dR_SCs=" << eleA_.dR_SCs(&eleB_) << "=" << eleB_.dR_SCs(&eleA_) << std::endl;
		std::cout << "        EleA:" << std::endl;
		std::cout << "          p4,E=" << eleA_.p4().E() << "; p4,p=(" << eleA_.p4().Px() << ", " << eleA_.p4().Py() << ", " << eleA_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleA_.et() << "; scEt=" << eleA_.scEt() << "; energy=" << eleA_.energy() << "; eta=" << eleA_.eta() << "; phi=" << eleA_.phi() << "; sc(eta,phi) = (" << eleA_.scEta() << ", " << eleA_.scPhi() << ")" << std::endl;
		std::cout << "          ptVtx=" << eleA_.ptVtx() << "; ptCalo=" << eleA_.ptCalo() << "; closestCtfTrk... (pT,eta,phi)=(" << eleA_.closestCtfTrk_pt() << ", " << eleA_.closestCtfTrk_eta() << ", " << eleA_.closestCtfTrk_phi() << ")" << std::endl;
		std::cout << "          gsfTrk... (eta, phi) = (" << eleA_.gsfTrk_eta() << ", " << eleA_.gsfTrk_phi() << "), vz=" << eleA_.gsfTrk_vz() << std::endl;
		std::cout << "        EleB:" << std::endl;
		std::cout << "          p4,E=" << eleB_.p4().E() << "; p4,p=(" << eleB_.p4().Px() << ", " << eleB_.p4().Py() << ", " << eleB_.p4().Pz() << std::endl;
		std::cout << "          et=" << eleB_.et() << "; scEt=" << eleB_.scEt() << "; energy=" << eleB_.energy() << "; eta=" << eleB_.eta() << "; phi=" << eleB_.phi() << "; sc(eta,phi) = (" << eleB_.scEta() << ", " << eleB_.scPhi() << ")" << std::endl;
		std::cout << "          ptVtx=" << eleB_.ptVtx() << "; ptCalo=" << eleB_.ptCalo() << "; closestCtfTrk... (pT,eta,phi)=(" << eleB_.closestCtfTrk_pt() << ", " << eleB_.closestCtfTrk_eta() << ", " << eleB_.closestCtfTrk_phi() << ")" << std::endl;
		std::cout << "          gsfTrk... (eta, phi) = (" << eleB_.gsfTrk_eta() << ", " << eleB_.gsfTrk_phi() << "), vz=" << eleB_.gsfTrk_vz() << std::endl;

	}

}// End of tsw namespace.

#endif
