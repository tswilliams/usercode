#ifndef tswDiEleDistns_h
#define tswDiEleDistns_h

namespace tsw{
	class DiEleDistns{
		public:
			//Constructors and destructors ...
			DiEleDistns(TString hNamePrefix, const TString titleStr_diEleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt);
			~DiEleDistns();

			//Method to fill the histograms
			void FillHistos(tsw::HEEPDiEle theDiEle, const Double_t dieleWeight);

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);

			//Method to write the histograms to file ...
			void WriteHistos(TFile* ptr_outFile);

			//Method to get vector of pointers to all histograms ...
			std::vector<TH1D*> GetPtrsToHistos();

		private:
			TString hNamePrefix_;
			TH1D* diEleHist_invMass_;
			TH1D* diEleHist_sumEt_;
			TH1D* diEleHist_sumEtLog_;
			TH1D* diEleHist_openingAngle_;
			TH1D* diEleHist_deltaEta_;
			TH1D* diEleHist_deltaPhi_;
			TH1D* diEleHist_deltaR_;
			TH1D* eleAHist_et_;
			TH1D* eleAHist_energy_;
			TH1D* eleAHist_eta_;
			TH1D* eleAHist_phi_;
			TH1D* eleBHist_et_;
			TH1D* eleBHist_energy_;
			TH1D* eleBHist_eta_;
			TH1D* eleBHist_phi_;
	};

	/////////////////////
	// Constructor
	DiEleDistns::DiEleDistns(TString hNamePrefix, const TString titleStr_diEleAndEleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt):
		hNamePrefix_(hNamePrefix)
	{
		//Form standard phrases ...
		TString str_eleType_GSFeles_cutsPhrase = titleStr_diEleAndEleType + " electrons " + titleStr_cutsPhrase;
		TString str_NoOfDiEles = "Number of electron pairs";
		TString str_NoOfEles = "Number of electrons";

		/*if(hNamePrefix.EndsWith("5e32_"))
			hNamePrefix.Resize(hNamePrefix.Last('5'));
		else if(hNamePrefix.EndsWith("e33_"))
			hNamePrefix.Resize(hNamePrefix.Last('e')-1);
		else if(hNamePrefix.EndsWith("TrgEt80_") || hNamePrefix.EndsWith("TrgEt100_") || hNamePrefix.EndsWith("TrgEt120_"))
			hNamePrefix.Resize(hNamePrefix.Last('T'));*/

		//Initialise the di-ele histograms ...
		diEleHist_invMass_     = new TH1D(hNamePrefix + "invMass",   "Di-electron invariant mass distribution  (" + str_eleType_GSFeles_cutsPhrase + "); M_{ee} /GeVc^{-2};" + str_NoOfDiEles + " per 2.5GeV", hNBins_mass, 45.0, 135.0);
		Double_t hMax_DiEleEt = 1200.0;
		diEleHist_sumEt_       = new TH1D(hNamePrefix + "sumEt",     "Di-electron p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, ee} /GeVc^{-1};" + str_NoOfDiEles, hNBins_pt, 0.0, hMax_pt);

		Double_t hBinLims_ptLog[121]; hBinLims_ptLog[0] = 0.1;
		for(unsigned int i=1; i<121; i++){
			hBinLims_ptLog[i] = pow(hMax_pt/hBinLims_ptLog[0], 1.0/120.0)*hBinLims_ptLog[i-1];
		}
		diEleHist_sumEtLog_    = new TH1D(hNamePrefix + "sumEtLog",     "Di-electron p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, ee} /GeVc^{-1};" + str_NoOfDiEles, hNBins_pt, hBinLims_ptLog);
		diEleHist_openingAngle_= new TH1D(hNamePrefix + "openAngle", "Di-electron opening angle distribution  (" + str_eleType_GSFeles_cutsPhrase + "); Di-electron opening angle, #theta_{ee} /rad;" + str_NoOfDiEles, 30, 0.0, 3.1416);
		diEleHist_deltaEta_    = new TH1D(hNamePrefix + "deltaEta",  "#Delta#eta_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#eta_{ee};" + str_NoOfDiEles, 50, -5.0, 5.0);
		diEleHist_deltaPhi_    = new TH1D(hNamePrefix + "deltaPhi",  "#Delta#phi_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#phi_{ee};" + str_NoOfDiEles, 30, -3.15, 3.15);
		diEleHist_deltaR_      = new TH1D(hNamePrefix + "deltaR",  "#Delta{}R_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta{}R_{ee};" + str_NoOfDiEles, 50, 0.0, 10.0);

		//Initialise the single ele histos ...
		Int_t hNBins_eleEt = 25; Double_t hMin_eleEt = 0.0; Double_t hMax_eleEt = 100.0;
		eleAHist_et_ = new TH1D(hNamePrefix + "eleA_et",  "eleA E_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A E_{T} /GeV", hNBins_eleEt, hMin_eleEt, hMax_eleEt);
		eleBHist_et_ = new TH1D(hNamePrefix + "eleB_et",  "eleB E_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B E_{T} /GeV", hNBins_eleEt, hMin_eleEt, hMax_eleEt);

		Int_t hNBins_eleEnergy = 25; Double_t hMin_eleEnergy = 0.0; Double_t hMax_eleEnergy = 100.0;
		eleAHist_energy_ = new TH1D(hNamePrefix + "eleA_energy",  "eleA energy distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A E /GeV", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);
		eleBHist_energy_ = new TH1D(hNamePrefix + "eleB_energy",  "eleB energy distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B E /GeV", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);

		Int_t hNBins_eleEta = 30; Double_t hMin_eleEta = -3.0; Double_t hMax_eleEta = +3.0;
		eleAHist_eta_ = new TH1D(hNamePrefix + "eleA_eta",  "eleA #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);
		eleBHist_eta_ = new TH1D(hNamePrefix + "eleB_eta",  "eleB #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);

		Int_t hNBins_elePhi = 30; Double_t hMin_elePhi = -3.1416; Double_t hMax_elePhi = +3.1416;
		eleAHist_phi_ = new TH1D(hNamePrefix + "eleA_phi",  "eleA #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);
		eleBHist_phi_ = new TH1D(hNamePrefix + "eleB_phi",  "eleB #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);

		SetHistAttributes(lineColorIdx, lineStyleIdx);
	}

	/////////////////////
	// Destructor
	DiEleDistns::~DiEleDistns(){
		// Delete all of the histograms from the heap ...
		delete diEleHist_invMass_;
		delete diEleHist_sumEt_;
		delete diEleHist_sumEtLog_;
		delete diEleHist_openingAngle_;
		delete diEleHist_deltaEta_;
		delete diEleHist_deltaPhi_;
		delete diEleHist_deltaR_;

		delete eleAHist_et_;
		delete eleAHist_energy_;
		delete eleAHist_eta_;
		delete eleAHist_phi_;

		delete eleBHist_et_;
		delete eleBHist_energy_;
		delete eleBHist_eta_;
		delete eleBHist_phi_;
	}

	/////////////////////////////////
	// Method to fill histograms ...
	void DiEleDistns::FillHistos(tsw::HEEPDiEle theDiEle, const Double_t dieleWeight){
		//Fill the di-ele histos ...
		diEleHist_invMass_->Fill(      theDiEle.invMass()      , dieleWeight);
		diEleHist_sumEt_->Fill(        theDiEle.et()           , dieleWeight);
		diEleHist_sumEtLog_->Fill(     theDiEle.et()           , dieleWeight);
		diEleHist_openingAngle_->Fill( theDiEle.openingAngle() , dieleWeight);
		diEleHist_deltaEta_->Fill(     theDiEle.deltaEta()     , dieleWeight);
		diEleHist_deltaPhi_->Fill(     theDiEle.deltaPhi()     , dieleWeight);
		diEleHist_deltaR_->Fill(       theDiEle.deltaR()       , dieleWeight);

		//Fill the component ele histos ...
		eleAHist_et_->Fill(     theDiEle.eleA().et()     , dieleWeight);
		eleAHist_energy_->Fill( theDiEle.eleA().energy() , dieleWeight);
		eleAHist_eta_->Fill(    theDiEle.eleA().eta()    , dieleWeight);
		eleAHist_phi_->Fill(    theDiEle.eleA().phi()    , dieleWeight);
		eleBHist_et_->Fill(     theDiEle.eleB().et()     , dieleWeight);
		eleBHist_energy_->Fill( theDiEle.eleB().energy() , dieleWeight);
		eleBHist_eta_->Fill(    theDiEle.eleB().eta()    , dieleWeight);
		eleBHist_phi_->Fill(    theDiEle.eleB().phi()    , dieleWeight);
	}

	void DiEleDistns::SetHistAttributes(Int_t lineCol, Int_t lineStyle){
		std::vector<TH1D*> ptrsToHists;
		ptrsToHists = GetPtrsToHistos();

		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->SetLineColor(lineCol);
			ptrsToHists.at(iHist)->SetLineStyle(lineStyle);
		}
	}

	/////////////////////////////////////////////////
	// Method to write out the histograms to file ...
	void DiEleDistns::WriteHistos(TFile* ptr_outFile){
		//Create a directory within the output file to store these histos in ...
		TDirectory* tmpTDirPtr;
		tmpTDirPtr = ptr_outFile->mkdir(hNamePrefix_);
		ptr_outFile->cd(hNamePrefix_);

		//Write out the di-ele histos ...
		diEleHist_invMass_->Write();
		diEleHist_sumEt_->Write();
		diEleHist_sumEtLog_->Write();
		diEleHist_openingAngle_->Write();
		diEleHist_deltaEta_->Write();
		diEleHist_deltaPhi_->Write();
		diEleHist_deltaR_->Write();

		//Write out the eleA_ histos ...
		eleAHist_et_->Write();
		eleAHist_energy_->Write();
		eleAHist_eta_->Write();
		eleAHist_phi_->Write();

		//Write out the eleB_ histos ...
		eleBHist_et_->Write();
		eleBHist_energy_->Write();
		eleBHist_eta_->Write();
		eleBHist_phi_->Write();

		//Move back to the original directory within the output file ...
		ptr_outFile->cd("..");
	}

	std::vector<TH1D*> DiEleDistns::GetPtrsToHistos(){
		std::vector<TH1D*> tmpPtrsVector;
		tmpPtrsVector.clear();

		tmpPtrsVector.push_back(diEleHist_invMass_);
		tmpPtrsVector.push_back(diEleHist_sumEt_);
		tmpPtrsVector.push_back(diEleHist_sumEtLog_);
		tmpPtrsVector.push_back(diEleHist_openingAngle_);
		tmpPtrsVector.push_back(diEleHist_deltaEta_);
		tmpPtrsVector.push_back(diEleHist_deltaPhi_);
		tmpPtrsVector.push_back(diEleHist_deltaR_);
		tmpPtrsVector.push_back(eleAHist_et_);
		tmpPtrsVector.push_back(eleAHist_energy_);
		tmpPtrsVector.push_back(eleAHist_eta_);
		tmpPtrsVector.push_back(eleAHist_phi_);
		tmpPtrsVector.push_back(eleBHist_et_);
		tmpPtrsVector.push_back(eleBHist_energy_);
		tmpPtrsVector.push_back(eleBHist_eta_);
		tmpPtrsVector.push_back(eleBHist_phi_);

		return tmpPtrsVector;
	}
}

#endif
