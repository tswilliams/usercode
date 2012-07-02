#ifndef tswEleMuDistns_h
#define tswEleMuDistns_h

namespace tsw{
	class EleMuDistns{
		public:
			//Constructors and destructors ...
			EleMuDistns(TString hNamePrefix, const TString titleStr_eleMuType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt);
			~EleMuDistns();

			//Method to fill the histograms
			void FillHistos(tsw::EleMuObject theEleMuObj, const Double_t eleMuWeight);

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);

			//Method to write the histograms to file ...
			void WriteHistos(TFile* ptr_outFile);

			//Method to get vector of pointers to all histograms ...
			std::vector<TH1D*> GetPtrsToHistos();

		private:
			TString hNamePrefix_;
			TH1D* eleMuHist_invMass_;
			TH1D* eleMuHist_sumEt_;
			TH1D* eleMuHist_sumEtLog_;
			TH1D* eleMuHist_openingAngle_;
			TH1D* eleMuHist_deltaEta_;
			TH1D* eleMuHist_deltaPhi_;
			TH1D* eleMuHist_deltaR_;
			TH1D* eleHist_et_;
			TH1D* eleHist_energy_;
			TH1D* eleHist_eta_;
			TH1D* eleHist_phi_;
			TH1D* muHist_pt_;
			TH1D* muHist_p_;
			TH1D* muHist_eta_;
			TH1D* muHist_phi_;
	};

	/////////////////////
	// Constructor
	EleMuDistns::EleMuDistns(TString hNamePrefix, const TString titleStr_eleMuObjEleAndMuType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt):
		hNamePrefix_(hNamePrefix)
	{
		//Form standard phrases ...
		TString str_eleType_GSFeles_cutsPhrase = titleStr_eleMuObjEleAndMuType + ", " + titleStr_cutsPhrase;
		TString str_NoOfEleMuObjs = "Number of e#mu pairs";
		TString str_NoOfEles = "Number of electrons";
		TString str_NoOfMuons = "Number of muons";

		//Initialise the di-ele histograms ...
		eleMuHist_invMass_     = new TH1D(hNamePrefix + "invMass",   "e#mu mass distribution (" + str_eleType_GSFeles_cutsPhrase + "); M_{e#mu} /GeVc^{-2};" + str_NoOfEleMuObjs + " per 2.5GeV", hNBins_mass, 45.0, 135.0);
		eleMuHist_sumEt_       = new TH1D(hNamePrefix + "sumEt",     "e#mu p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, e#mu} /GeVc^{-1};" + str_NoOfEleMuObjs, hNBins_pt, 0.0, hMax_pt);

		Double_t hBinLims_ptLog[121]; hBinLims_ptLog[0] = 0.1;
		for(unsigned int i=1; i<121; i++){
			hBinLims_ptLog[i] = pow(hMax_pt/hBinLims_ptLog[0], 1.0/120.0)*hBinLims_ptLog[i-1];
		}
		eleMuHist_sumEtLog_    = new TH1D(hNamePrefix + "sumEtLog",     "e#mu p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, e#mu} /GeVc^{-1};" + str_NoOfEleMuObjs, hNBins_pt, hBinLims_ptLog);
		eleMuHist_openingAngle_= new TH1D(hNamePrefix + "openAngle", "e#mu opening angle distribution  (" + str_eleType_GSFeles_cutsPhrase + "); e#mu opening angle, #theta_{e#mu} /rad;" + str_NoOfEleMuObjs, 30, 0.0, 3.1416);
		eleMuHist_deltaEta_    = new TH1D(hNamePrefix + "deltaEta",  "#Delta#eta_{e#mu} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#eta_{e#mu};" + str_NoOfEleMuObjs, 50, -5.0, 5.0);
		eleMuHist_deltaPhi_    = new TH1D(hNamePrefix + "deltaPhi",  "#Delta#phi_{e#mu} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#phi_{e#mu};" + str_NoOfEleMuObjs, 30, -3.15, 3.15);
		eleMuHist_deltaR_      = new TH1D(hNamePrefix + "deltaR",  "#Delta R_{e#mu} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta R_{e#mu};" + str_NoOfEleMuObjs, 50, 0.0, 2.5);

		//Initialise the single ele histos ...
		Int_t hNBins_eleEt = hNBins_pt; Double_t hMin_eleEt = 0.0; Double_t hMax_eleEt = hMax_pt/2.0;
		eleHist_et_ = new TH1D(hNamePrefix + "ele_et",  "Electron E_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; Electron E_{T} /GeV", hNBins_eleEt, hMin_eleEt, hMax_eleEt);
		muHist_pt_ = new TH1D(hNamePrefix + "mu_pt",  "Muon p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfMuons + "; Muon p_{T} /GeVc^{-1}", hNBins_eleEt, hMin_eleEt, hMax_eleEt);

		Int_t hNBins_eleEnergy = hNBins_pt; Double_t hMin_eleEnergy = 0.0; Double_t hMax_eleEnergy = hMax_pt/2.0 ;
		eleHist_energy_ = new TH1D(hNamePrefix + "ele_energy",  "Electron energy distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; Electron energy, E /GeV", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);
		muHist_p_ = new TH1D(hNamePrefix + "mu_p",  "Muon momentum distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfMuons + "; Muon momentum, p /GeVc^{-1}", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);

		Int_t hNBins_eleEta = 30; Double_t hMin_eleEta = -3.0; Double_t hMax_eleEta = +3.0;
		eleHist_eta_ = new TH1D(hNamePrefix + "ele_eta",  "Electron #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; Electron #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);
		muHist_eta_ = new TH1D(hNamePrefix + "mu_eta",  "Muon #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfMuons + "; Muon #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);

		Int_t hNBins_elePhi = 30; Double_t hMin_elePhi = -3.1416; Double_t hMax_elePhi = +3.1416;
		eleHist_phi_ = new TH1D(hNamePrefix + "ele_phi",  "Electron #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; Electron #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);
		muHist_phi_ = new TH1D(hNamePrefix + "mu_phi",  "Muon #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfMuons + "; Muon #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);

		SetHistAttributes(lineColorIdx, lineStyleIdx);

		// Set-up histograms so that errors are automatically calculated as the histos are filled
		std::vector<TH1D*> ptrsToHists = GetPtrsToHistos();
		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->Sumw2();
			ptrsToHists.at(iHist)->SetDirectory(0);
		}
	}

	/////////////////////
	// Destructor
	EleMuDistns::~EleMuDistns(){
		// Delete all of the histograms from the heap ...eleMuHist_invMass_;
		delete eleMuHist_invMass_;
		delete eleMuHist_sumEt_;
		delete eleMuHist_sumEtLog_;
		delete eleMuHist_openingAngle_;
		delete eleMuHist_deltaEta_;
		delete eleMuHist_deltaPhi_;
		delete eleMuHist_deltaR_;

		delete eleHist_et_;
		delete eleHist_energy_;
		delete eleHist_eta_;
		delete eleHist_phi_;

		delete muHist_pt_;
		delete muHist_p_;
		delete muHist_eta_;
		delete muHist_phi_;
	}

	/////////////////////////////////
	// Method to fill histograms ...
	void EleMuDistns::FillHistos(tsw::EleMuObject theEleMuObj, const Double_t eleMuWeight){
		//Fill the di-ele histos ...
		eleMuHist_invMass_->Fill(      theEleMuObj.mass()        , eleMuWeight);
		eleMuHist_sumEt_->Fill(        theEleMuObj.pT()          , eleMuWeight);
		eleMuHist_sumEtLog_->Fill(     theEleMuObj.pT()          , eleMuWeight);
		eleMuHist_openingAngle_->Fill( theEleMuObj.openingAngle(), eleMuWeight);
		eleMuHist_deltaEta_->Fill(     theEleMuObj.deltaEta()    , eleMuWeight);
		eleMuHist_deltaPhi_->Fill(     theEleMuObj.deltaPhi()    , eleMuWeight);
		eleMuHist_deltaR_->Fill(       theEleMuObj.deltaR()      , eleMuWeight);

		//Fill the component ele histos ...
		eleHist_et_->Fill(     theEleMuObj.GetElectron()->et()     , eleMuWeight);
		eleHist_energy_->Fill( theEleMuObj.GetElectron()->energy() , eleMuWeight);
		eleHist_eta_->Fill(    theEleMuObj.GetElectron()->eta()    , eleMuWeight);
		eleHist_phi_->Fill(    theEleMuObj.GetElectron()->phi()    , eleMuWeight);
		muHist_pt_->Fill( theEleMuObj.GetMuon()->pT() , eleMuWeight);
		muHist_p_->Fill(  theEleMuObj.GetMuon()->p()  , eleMuWeight);
		muHist_eta_->Fill(theEleMuObj.GetMuon()->eta(), eleMuWeight);
		muHist_phi_->Fill(theEleMuObj.GetMuon()->phi(), eleMuWeight);
	}

	void EleMuDistns::SetHistAttributes(Int_t lineCol, Int_t lineStyle){
		std::vector<TH1D*> ptrsToHists;
		ptrsToHists = GetPtrsToHistos();

		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->SetLineColor(lineCol);
			ptrsToHists.at(iHist)->SetLineStyle(lineStyle);
		}
	}

	/////////////////////////////////////////////////
	// Method to write out the histograms to file ...
	void EleMuDistns::WriteHistos(TFile* ptr_outFile){
		//Create a directory within the output file to store these histos in ...
		TDirectory* tmpTDirPtr;
		tmpTDirPtr = ptr_outFile->mkdir(hNamePrefix_);
		ptr_outFile->cd(hNamePrefix_);

		//Write out the di-ele histos ...
		eleMuHist_invMass_->Write();
		eleMuHist_sumEt_->Write();
		eleMuHist_sumEtLog_->Write();
		eleMuHist_openingAngle_->Write();
		eleMuHist_deltaEta_->Write();
		eleMuHist_deltaPhi_->Write();
		eleMuHist_deltaR_->Write();

		//Write out the eleA_ histos ...
		eleHist_et_->Write();
		eleHist_energy_->Write();
		eleHist_eta_->Write();
		eleHist_phi_->Write();

		//Write out the eleB_ histos ...
		muHist_pt_->Write();
		muHist_p_->Write();
		muHist_eta_->Write();
		muHist_phi_->Write();

		//Move back to the original directory within the output file ...
		ptr_outFile->cd("..");
	}

	std::vector<TH1D*> EleMuDistns::GetPtrsToHistos(){
		std::vector<TH1D*> tmpPtrsVector;
		tmpPtrsVector.clear();

		tmpPtrsVector.push_back(eleMuHist_invMass_);
		tmpPtrsVector.push_back(eleMuHist_sumEt_);
		tmpPtrsVector.push_back(eleMuHist_sumEtLog_);
		tmpPtrsVector.push_back(eleMuHist_openingAngle_);
		tmpPtrsVector.push_back(eleMuHist_deltaEta_);
		tmpPtrsVector.push_back(eleMuHist_deltaPhi_);
		tmpPtrsVector.push_back(eleMuHist_deltaR_);
		tmpPtrsVector.push_back(eleHist_et_);
		tmpPtrsVector.push_back(eleHist_energy_);
		tmpPtrsVector.push_back(eleHist_eta_);
		tmpPtrsVector.push_back(eleHist_phi_);
		tmpPtrsVector.push_back(muHist_pt_);
		tmpPtrsVector.push_back(muHist_p_);
		tmpPtrsVector.push_back(muHist_eta_);
		tmpPtrsVector.push_back(muHist_phi_);

		return tmpPtrsVector;
	}
}

#endif
