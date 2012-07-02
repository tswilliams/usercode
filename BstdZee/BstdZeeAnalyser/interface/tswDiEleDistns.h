#ifndef tswDiEleDistns_h
#define tswDiEleDistns_h

namespace tsw{
	class DiEleDistns{
		public:
			//Constructors and destructors ...
			DiEleDistns(TString hNamePrefix, const TString titleStr_diEleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt);
			~DiEleDistns();

			//Method to fill the histograms
			void FillHistos(tsw::HEEPDiEle theDiEle, const Double_t dieleWeight, const bool );

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);

			//Method to write the histograms to file ...
			void WriteHistos(TFile* ptr_outFile);

			//Method to get vector of pointers to all histograms ...
			std::vector<TH1D*> GetPtrsToHistos();

		private:
			TString hNamePrefix_;
			TH1D* diEleHist_regions_;
			TH1D* diEleHist_invMass_;
			TH1D* diEleHist_sumEt_;
			TH1D* diEleHist_coarsePt_;
			TH1D* diEleHist_sumEtLog_;
			TH1D* diEleHist_openingAngle_;
			TH1D* diEleHist_deltaEta_;
			TH1D* diEleHist_deltaPhi_;
			TH1D* diEleHist_deltaR_;
			TH1D* diEleHist_EtAminEtB_;
			TH1D* diEleHist_EtAOverEtB_;
			TH1D* diEleHist_modTrkIso_;
			TH1D* diEleHist_dRleq03modTrkIso_;
			TH1D* diEleHist_dRleq03modTrkIsoNoCtfTrk_;

			TH1D* eleAHist_et_;
			TH1D* eleAHist_energy_;
			TH1D* eleAHist_eta_;
			TH1D* eleAHist_phi_;
			TH1D* eleAHist_modTrkIso_;
			TH1D* eleAHist_dRleq03modTrkIso_;
			TH1D* eleAHist_dRleq03modTrkIsoNoCtfTrk_;
			TH1D* eleAHist_modEmHad1Iso_;
			TH1D* eleAHist_modEmHad1IsoOverCutThr_;
			TH1D* eleAHist_modHad1IsoOverModEmHad1_;
			TH1D* eleAHist_Had1Iso_;
			TH1D* eleAHist_modHad1Iso_;
			TH1D* eleAHist_modHad1IsoOverHad1_;

			TH1D* eleBHist_et_;
			TH1D* eleBHist_energy_;
			TH1D* eleBHist_eta_;
			TH1D* eleBHist_phi_;
			TH1D* eleBHist_modTrkIso_;
			TH1D* eleBHist_dRleq03modTrkIso_;
			TH1D* eleBHist_dRleq03modTrkIsoNoCtfTrk_;
			TH1D* eleBHist_modEmHad1Iso_;
			TH1D* eleBHist_modEmHad1IsoOverCutThr_;
			TH1D* eleBHist_modHad1IsoOverModEmHad1_;
			TH1D* eleBHist_Had1Iso_;
			TH1D* eleBHist_modHad1Iso_;
			TH1D* eleBHist_modHad1IsoOverHad1_;
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
		diEleHist_regions_     = new TH1D(hNamePrefix + "regions",   "Di-electron type by region (" + str_eleType_GSFeles_cutsPhrase + "); Region;" + str_NoOfDiEles, 4, -0.5, 3.5);
		diEleHist_regions_->GetXaxis()->SetBinLabel(1, "EB-EB");
		diEleHist_regions_->GetXaxis()->SetBinLabel(2, "EB-EE");
		diEleHist_regions_->GetXaxis()->SetBinLabel(3, "EE-EE");
		diEleHist_regions_->GetXaxis()->SetBinLabel(4, "Other");
		diEleHist_invMass_     = new TH1D(hNamePrefix + "invMass",   "Di-electron invariant mass distribution  (" + str_eleType_GSFeles_cutsPhrase + "); M_{ee} /GeVc^{-2};" + str_NoOfDiEles + " per GeV", hNBins_mass, 50.0, 130.0);
		diEleHist_sumEt_       = new TH1D(hNamePrefix + "sumEt",     "Di-electron p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, ee} /GeVc^{-1};" + str_NoOfDiEles, hNBins_pt, 0.0, hMax_pt);

		Int_t    hNBins_ptCoarse = 8;
		Double_t hBinLims_ptCoarse[] = {0.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 500.0, 1000.0};
		diEleHist_coarsePt_    = new TH1D(hNamePrefix + "coarsePt",     "Di-electron p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, ee} /GeVc^{-1};" + str_NoOfDiEles, hNBins_ptCoarse, hBinLims_ptCoarse);

		Double_t hBinLims_ptLog[121]; hBinLims_ptLog[0] = 0.1;
		for(unsigned int i=1; i<121; i++){
			hBinLims_ptLog[i] = pow(hMax_pt/hBinLims_ptLog[0], 1.0/120.0)*hBinLims_ptLog[i-1];
		}
		diEleHist_sumEtLog_    = new TH1D(hNamePrefix + "sumEtLog",     "Di-electron p_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); p_{T, ee} /GeVc^{-1};" + str_NoOfDiEles, hNBins_pt, hBinLims_ptLog);
		diEleHist_openingAngle_= new TH1D(hNamePrefix + "openAngle", "Di-electron opening angle distribution  (" + str_eleType_GSFeles_cutsPhrase + "); Di-electron opening angle, #theta_{ee} /rad;" + str_NoOfDiEles, 30, 0.0, 3.1416);
		diEleHist_deltaEta_    = new TH1D(hNamePrefix + "deltaEta",  "#Delta#eta_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#eta_{ee};" + str_NoOfDiEles, 200, -5.0, 5.0);
		diEleHist_deltaPhi_    = new TH1D(hNamePrefix + "deltaPhi",  "#Delta#phi_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta#phi_{ee};" + str_NoOfDiEles, 126, -3.15, 3.15);
		diEleHist_deltaR_      = new TH1D(hNamePrefix + "deltaR",  "#Delta{}R_{ee} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); #Delta{}R_{ee};" + str_NoOfDiEles, 100, 0.0, 5.0);
		diEleHist_EtAminEtB_   = new TH1D(hNamePrefix + "EtAminEtB",  "E_{T}^{A}-E_{T}^{B} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); E_{T}^{A}-E_{T}^{B} /GeV;" + str_NoOfDiEles, 120, 0.0, 1200.0);
		diEleHist_EtAOverEtB_  = new TH1D(hNamePrefix + "EtAOverEtB", "E_{T}^{A}/E_{T}^{B} distribution  (" + str_eleType_GSFeles_cutsPhrase + "); E_{T}^{A}/E_{T}^{B};" + str_NoOfDiEles, 40, 0.0, 20.0);

		//Initialise the single ele histos ...
		Int_t hNBins_eleEt = hNBins_pt; Double_t hMin_eleEt = 0.0; Double_t hMax_eleEt = hMax_pt/2.0;
		eleAHist_et_ = new TH1D(hNamePrefix + "eleA_et",  "eleA E_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A E_{T} /GeV", hNBins_eleEt, hMin_eleEt, hMax_eleEt);
		eleBHist_et_ = new TH1D(hNamePrefix + "eleB_et",  "eleB E_{T} distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B E_{T} /GeV", hNBins_eleEt, hMin_eleEt, hMax_eleEt);

		Int_t hNBins_eleEnergy = hNBins_pt; Double_t hMin_eleEnergy = 0.0; Double_t hMax_eleEnergy = hMax_pt/2.0 ;
		eleAHist_energy_ = new TH1D(hNamePrefix + "eleA_energy",  "eleA energy distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A E /GeV", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);
		eleBHist_energy_ = new TH1D(hNamePrefix + "eleB_energy",  "eleB energy distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B E /GeV", hNBins_eleEnergy, hMin_eleEnergy, hMax_eleEnergy);

		Int_t hNBins_eleEta = 30; Double_t hMin_eleEta = -3.0; Double_t hMax_eleEta = +3.0;
		eleAHist_eta_ = new TH1D(hNamePrefix + "eleA_eta",  "eleA #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);
		eleBHist_eta_ = new TH1D(hNamePrefix + "eleB_eta",  "eleB #eta distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B #eta", hNBins_eleEta, hMin_eleEta, hMax_eleEta);

		Int_t hNBins_elePhi = 30; Double_t hMin_elePhi = -3.1416; Double_t hMax_elePhi = +3.1416;
		eleAHist_phi_ = new TH1D(hNamePrefix + "eleA_phi",  "eleA #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele A #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);
		eleBHist_phi_ = new TH1D(hNamePrefix + "eleB_phi",  "eleB #phi distribution  (" + str_eleType_GSFeles_cutsPhrase + ");" + str_NoOfEles + "; ele B #phi", hNBins_elePhi, hMin_elePhi, hMax_elePhi);

		//Initialise the tracker isolation histograms ...
		diEleHist_modTrkIso_ = new TH1D(hNamePrefix + "both_modTrkIso",  "Sum modTrkIso for both eles (" + str_eleType_GSFeles_cutsPhrase + "); eleA+eleB #Sigma p_{T}^{trk};" + str_NoOfDiEles, 80, -600.0, +600.0);
		diEleHist_dRleq03modTrkIso_ = new TH1D(hNamePrefix + "both_dRleq03modTrkIso",  "Sum modTrkIso for both eles if modified (" + str_eleType_GSFeles_cutsPhrase + "); eleA+eleB #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3;" + str_NoOfDiEles, 80, -600.0, +600.0);
		diEleHist_dRleq03modTrkIsoNoCtfTrk_ = new TH1D(hNamePrefix + "both_dRleq03modTrkIsoNoCtfTrk",  "Sum modTrkIso for both eles if modified AND no CTF track for 1 ele (" + str_eleType_GSFeles_cutsPhrase + "); eleA+eleB #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3 and 1 ele hasn't got closest CTF track;" + str_NoOfDiEles, 80, -600.0, +600.0);
		eleAHist_modTrkIso_ = new TH1D(hNamePrefix + "eleA_modTrkIso",  "eleA modTrkIso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma p_{T}^{trk};" + str_NoOfEles, 80, -300.0, +300.0);
		eleAHist_dRleq03modTrkIso_ = new TH1D(hNamePrefix + "eleA_dRleq03modTrkIso",  "eleA modTrkIso if modified distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3;" + str_NoOfEles, 80, -300.0, +300.0);
		eleAHist_dRleq03modTrkIsoNoCtfTrk_ = new TH1D(hNamePrefix + "eleA_dRleq03modTrkIsoNoCtfTrk",  "eleA modTrkIso if (modified AND no eleB CTF) distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3 and eleB hasn't got closest CTF track;" + str_NoOfEles, 80, -300.0, +300.0);
		eleBHist_modTrkIso_ = new TH1D(hNamePrefix + "eleB_modTrkIso",  "eleB modTrkIso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma p_{T}^{trk};" + str_NoOfEles, 80, -300.0, +300.0);
		eleBHist_dRleq03modTrkIso_ = new TH1D(hNamePrefix + "eleB_dRleq03modTrkIso",  "eleB modTrkIso if modified distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3;" + str_NoOfEles, 80, -300.0, +300.0);
		eleBHist_dRleq03modTrkIsoNoCtfTrk_ = new TH1D(hNamePrefix + "eleB_dRleq03modTrkIsoNoCtfTrk",  "eleB modTrkIso if (modified AND no eleA CTF) distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma p_{T}^{trk} if #Delta R_{ee}<0.3 and eleA hasn't got closest CTF track;" + str_NoOfEles, 80, -300.0, +300.0);

		//Initialise the EmHad1 isolation histograms ...
		eleAHist_modEmHad1Iso_ = new TH1D(hNamePrefix + "eleA_modEmHad1Iso",  "eleA modEmHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{EmHad1};" + str_NoOfEles, 200, -100.0, +100.0);
		eleAHist_modEmHad1IsoOverCutThr_ = new TH1D(hNamePrefix + "eleA_modEmHad1IsoOverCutThr",  "eleA modEmHad1Iso/cut distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{EmHad1} / Cut thr;" + str_NoOfEles, 80, -8.0, +8.0);
		eleAHist_modHad1IsoOverModEmHad1_ = new TH1D(hNamePrefix + "eleA_modHad1IsoOverModEmHad1",  "eleA Had1Iso/modEmHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{Had1} / #Sigma E_{T}^{EmHad1};" + str_NoOfEles, 50, 0.0, +1.0);
		eleBHist_modEmHad1Iso_ = new TH1D(hNamePrefix + "eleB_modEmHad1Iso",  "eleB modEmHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{EmHad1};" + str_NoOfEles, 200, -100.0, +100.0);
		eleBHist_modEmHad1IsoOverCutThr_ = new TH1D(hNamePrefix + "eleB_modEmHad1IsoOverCutThr",  "eleB modEmHad1Iso/cut distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{EmHad1} / Cut thr;" + str_NoOfEles, 80, -8.0, +8.0);
		eleBHist_modHad1IsoOverModEmHad1_ = new TH1D(hNamePrefix + "eleB_modHad1IsoOverModEmHad1",  "eleB Had1Iso/modEmHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{Had1} / #Sigma E_{T}^{EmHad1};" + str_NoOfEles, 50, 0.0, +1.0);

		//Initialise the Had1 isolation histograms ...
		eleAHist_Had1Iso_ = new TH1D(hNamePrefix + "eleA_Had1Iso",  "eleA Had1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{Had1};" + str_NoOfEles, 200, -10.0, +10.0);
		eleBHist_Had1Iso_ = new TH1D(hNamePrefix + "eleB_Had1Iso",  "eleB Had1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{Had1};" + str_NoOfEles, 200, -10.0, +10.0);
		eleAHist_modHad1Iso_ = new TH1D(hNamePrefix + "eleA_modHad1Iso",  "eleA modHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{modHad1};" + str_NoOfEles, 200, -10.0, +10.0);
		eleBHist_modHad1Iso_ = new TH1D(hNamePrefix + "eleB_modHad1Iso",  "eleB modHad1Iso distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{modHad1};" + str_NoOfEles, 200, -10.0, +10.0);
		eleAHist_modHad1IsoOverHad1_ = new TH1D(hNamePrefix + "eleA_modHad1IsoOverHad1",  "eleA modHad1IsoOverHad1 distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleA #Sigma E_{T}^{modHad1} / #Sigma E_{T}^{Had1};" + str_NoOfEles, 200, -4.5, +1.5);
		eleBHist_modHad1IsoOverHad1_ = new TH1D(hNamePrefix + "eleB_modHad1IsoOverHad1",  "eleB modHad1IsoOverHad1 distribution (" + str_eleType_GSFeles_cutsPhrase + "); eleB #Sigma E_{T}^{modHad1} / #Sigma E_{T}^{Had1};" + str_NoOfEles, 200, -4.5, +1.5);

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
	DiEleDistns::~DiEleDistns(){
		// Delete all of the histograms from the heap ...
		delete diEleHist_regions_;
		delete diEleHist_invMass_;
		delete diEleHist_sumEt_;
		delete diEleHist_coarsePt_;
		delete diEleHist_sumEtLog_;
		delete diEleHist_openingAngle_;
		delete diEleHist_deltaEta_;
		delete diEleHist_deltaPhi_;
		delete diEleHist_deltaR_;
		delete diEleHist_EtAminEtB_;
		delete diEleHist_EtAOverEtB_;
		delete diEleHist_modTrkIso_;
		delete diEleHist_dRleq03modTrkIso_;
		delete diEleHist_dRleq03modTrkIsoNoCtfTrk_;

		delete eleAHist_et_;
		delete eleAHist_energy_;
		delete eleAHist_eta_;
		delete eleAHist_phi_;
		delete eleAHist_modTrkIso_;
		delete eleAHist_dRleq03modTrkIso_;
		delete eleAHist_dRleq03modTrkIsoNoCtfTrk_;
		delete eleAHist_modEmHad1Iso_;
		delete eleAHist_modEmHad1IsoOverCutThr_;
		delete eleAHist_modHad1IsoOverModEmHad1_;
		delete eleAHist_Had1Iso_;
		delete eleAHist_modHad1Iso_;
		delete eleAHist_modHad1IsoOverHad1_;

		delete eleBHist_et_;
		delete eleBHist_energy_;
		delete eleBHist_eta_;
		delete eleBHist_phi_;
		delete eleBHist_modTrkIso_;
		delete eleBHist_dRleq03modTrkIso_;
		delete eleBHist_dRleq03modTrkIsoNoCtfTrk_;
		delete eleBHist_modEmHad1Iso_;
		delete eleBHist_modEmHad1IsoOverCutThr_;
		delete eleBHist_modHad1IsoOverModEmHad1_;
		delete eleBHist_Had1Iso_;
		delete eleBHist_modHad1Iso_;
		delete eleBHist_modHad1IsoOverHad1_;
	}

	/////////////////////////////////
	// Method to fill histograms ...
	void DiEleDistns::FillHistos(tsw::HEEPDiEle theDiEle, const Double_t dieleWeight, const bool fillIsolVarHistos = false){
		//Fill the di-ele histos ...
		if( theDiEle.isEBEB() )
			diEleHist_regions_->Fill(0.0, dieleWeight);
		else if( theDiEle.isEBEE() )
			diEleHist_regions_->Fill(1.0, dieleWeight);
		else if( theDiEle.isEEEE() )
			diEleHist_regions_->Fill(2.0, dieleWeight);
		else
			diEleHist_regions_->Fill(3.0, dieleWeight);
		diEleHist_invMass_->Fill(      theDiEle.invMass()      , dieleWeight);
		diEleHist_sumEt_->Fill(        theDiEle.pT()           , dieleWeight);
		diEleHist_coarsePt_->Fill(     theDiEle.pT()           , dieleWeight);
		diEleHist_sumEtLog_->Fill(     theDiEle.pT()           , dieleWeight);
		diEleHist_openingAngle_->Fill( theDiEle.openingAngle() , dieleWeight);
		diEleHist_deltaEta_->Fill(     theDiEle.deltaEta()     , dieleWeight);
		diEleHist_deltaPhi_->Fill(     theDiEle.deltaPhi()     , dieleWeight);
		diEleHist_deltaR_->Fill(       theDiEle.deltaR()       , dieleWeight);

		diEleHist_EtAminEtB_->Fill(  (theDiEle.eleA().et()-theDiEle.eleB().et()), dieleWeight );
		diEleHist_EtAOverEtB_->Fill( (theDiEle.eleA().et()/theDiEle.eleB().et()), dieleWeight );

		//Fill the component ele histos ...
		eleAHist_et_->Fill(     theDiEle.eleA().et()     , dieleWeight);
		eleAHist_energy_->Fill( theDiEle.eleA().energy() , dieleWeight);
		eleAHist_eta_->Fill(    theDiEle.eleA().eta()    , dieleWeight);
		eleAHist_phi_->Fill(    theDiEle.eleA().phi()    , dieleWeight);
		eleBHist_et_->Fill(     theDiEle.eleB().et()     , dieleWeight);
		eleBHist_energy_->Fill( theDiEle.eleB().energy() , dieleWeight);
		eleBHist_eta_->Fill(    theDiEle.eleB().eta()    , dieleWeight);
		eleBHist_phi_->Fill(    theDiEle.eleB().phi()    , dieleWeight);

		if(fillIsolVarHistos){
			//Filling the tracker isolation histos ...
			diEleHist_modTrkIso_->Fill( (theDiEle.eleA_modTrkIso()+theDiEle.eleB_modTrkIso()), dieleWeight);
			eleAHist_modTrkIso_->Fill( theDiEle.eleA_modTrkIso(), dieleWeight);
			eleBHist_modTrkIso_->Fill( theDiEle.eleB_modTrkIso(), dieleWeight);
			if( theDiEle.deltaR()<0.3 ){
				diEleHist_dRleq03modTrkIso_->Fill( (theDiEle.eleA_modTrkIso()+theDiEle.eleB_modTrkIso()), dieleWeight);
				eleAHist_dRleq03modTrkIso_->Fill( theDiEle.eleA_modTrkIso(), dieleWeight);
				eleBHist_dRleq03modTrkIso_->Fill( theDiEle.eleB_modTrkIso(), dieleWeight);
				if( (!(theDiEle.eleA().closestCtfTrk_exists())) || (!(theDiEle.eleB().closestCtfTrk_exists())) )
					diEleHist_dRleq03modTrkIsoNoCtfTrk_->Fill( (theDiEle.eleA_modTrkIso()+theDiEle.eleB_modTrkIso()), dieleWeight);
				if( !(theDiEle.eleB().closestCtfTrk_exists()) )
					eleAHist_dRleq03modTrkIsoNoCtfTrk_->Fill( theDiEle.eleA_modTrkIso(), dieleWeight);
				if( !(theDiEle.eleA().closestCtfTrk_exists()) )
					eleBHist_dRleq03modTrkIsoNoCtfTrk_->Fill( theDiEle.eleB_modTrkIso(), dieleWeight);
			}

			// Filling EmHad1 isolation histos ...
			eleAHist_modEmHad1Iso_->Fill(theDiEle.eleA_modEmHad1Iso(), dieleWeight);
			eleAHist_modEmHad1IsoOverCutThr_->Fill( (theDiEle.eleA_modEmHad1Iso()/(2.0+0.03*theDiEle.eleA().et())), dieleWeight);
			//eleAHist_modHad1IsoOverModEmHad1_->Fill(  (theDiEle.eleA().isolHadDepth1()/theDiEle.eleA_modEmHad1Iso()), dieleWeight);
			eleAHist_modHad1IsoOverModEmHad1_->Fill(  (theDiEle.eleA_modHad1Iso()/theDiEle.eleA_modEmHad1Iso()), dieleWeight);
			eleBHist_modEmHad1Iso_->Fill(theDiEle.eleB_modEmHad1Iso(), dieleWeight);
			eleBHist_modEmHad1IsoOverCutThr_->Fill( (theDiEle.eleB_modEmHad1Iso()/(2.0+0.03*theDiEle.eleB().et())), dieleWeight);
			//eleBHist_modHad1IsoOverModEmHad1_->Fill(  (theDiEle.eleB().isolHadDepth1()/theDiEle.eleB_modEmHad1Iso()), dieleWeight);
			eleBHist_modHad1IsoOverModEmHad1_->Fill(  (theDiEle.eleB_modHad1Iso()/theDiEle.eleB_modEmHad1Iso()), dieleWeight);

			// Filling Had1 isolation histos ...
			eleAHist_Had1Iso_->Fill(theDiEle.eleA().isolHadDepth1(), dieleWeight);
			eleBHist_Had1Iso_->Fill(theDiEle.eleB().isolHadDepth1(), dieleWeight);
			eleAHist_modHad1Iso_->Fill(theDiEle.eleA_modHad1Iso(), dieleWeight);
			eleBHist_modHad1Iso_->Fill(theDiEle.eleB_modHad1Iso(), dieleWeight);
			eleAHist_modHad1IsoOverHad1_->Fill( (theDiEle.eleA_modHad1Iso()/theDiEle.eleA().isolHadDepth1()), dieleWeight);
			eleBHist_modHad1IsoOverHad1_->Fill( (theDiEle.eleB_modHad1Iso()/theDiEle.eleB().isolHadDepth1()), dieleWeight);
		}
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
		diEleHist_regions_->Write();
		diEleHist_invMass_->Write();
		diEleHist_sumEt_->Write();
		diEleHist_coarsePt_->Write();
		diEleHist_sumEtLog_->Write();
		diEleHist_openingAngle_->Write();
		diEleHist_deltaEta_->Write();
		diEleHist_deltaPhi_->Write();
		diEleHist_deltaR_->Write();
		diEleHist_EtAminEtB_->Write();
		diEleHist_EtAOverEtB_->Write();
		diEleHist_modTrkIso_->Write();
		diEleHist_dRleq03modTrkIso_->Write();
		diEleHist_dRleq03modTrkIsoNoCtfTrk_->Write();

		//Write out the eleA_ histos ...
		eleAHist_et_->Write();
		eleAHist_energy_->Write();
		eleAHist_eta_->Write();
		eleAHist_phi_->Write();
		eleAHist_modTrkIso_->Write();
		eleAHist_dRleq03modTrkIso_->Write();
		eleAHist_dRleq03modTrkIsoNoCtfTrk_->Write();
		eleAHist_modEmHad1Iso_->Write();
		eleAHist_modEmHad1IsoOverCutThr_->Write();
		eleAHist_modHad1IsoOverModEmHad1_->Write();
		eleAHist_Had1Iso_->Write();
		eleAHist_modHad1Iso_->Write();
		eleAHist_modHad1IsoOverHad1_->Write();

		//Write out the eleB_ histos ...
		eleBHist_et_->Write();
		eleBHist_energy_->Write();
		eleBHist_eta_->Write();
		eleBHist_phi_->Write();
		eleBHist_modTrkIso_->Write();
		eleBHist_dRleq03modTrkIso_->Write();
		eleBHist_dRleq03modTrkIsoNoCtfTrk_->Write();
		eleBHist_modEmHad1Iso_->Write();
		eleBHist_modEmHad1IsoOverCutThr_->Write();
		eleBHist_modHad1IsoOverModEmHad1_->Write();
		eleBHist_Had1Iso_->Write();
		eleBHist_modHad1Iso_->Write();
		eleBHist_modHad1IsoOverHad1_->Write();

		//Move back to the original directory within the output file ...
		ptr_outFile->cd("..");
	}

	std::vector<TH1D*> DiEleDistns::GetPtrsToHistos(){
		std::vector<TH1D*> tmpPtrsVector;
		tmpPtrsVector.clear();

		tmpPtrsVector.push_back(diEleHist_regions_);
		tmpPtrsVector.push_back(diEleHist_invMass_);
		tmpPtrsVector.push_back(diEleHist_sumEt_);
		tmpPtrsVector.push_back(diEleHist_coarsePt_);
		tmpPtrsVector.push_back(diEleHist_sumEtLog_);
		tmpPtrsVector.push_back(diEleHist_openingAngle_);
		tmpPtrsVector.push_back(diEleHist_deltaEta_);
		tmpPtrsVector.push_back(diEleHist_deltaPhi_);
		tmpPtrsVector.push_back(diEleHist_deltaR_);
		tmpPtrsVector.push_back(diEleHist_EtAminEtB_);
		tmpPtrsVector.push_back(diEleHist_EtAOverEtB_);
		tmpPtrsVector.push_back(diEleHist_modTrkIso_);
		tmpPtrsVector.push_back(diEleHist_dRleq03modTrkIso_);
		tmpPtrsVector.push_back(diEleHist_dRleq03modTrkIsoNoCtfTrk_);

		tmpPtrsVector.push_back(eleAHist_et_);
		tmpPtrsVector.push_back(eleAHist_energy_);
		tmpPtrsVector.push_back(eleAHist_eta_);
		tmpPtrsVector.push_back(eleAHist_phi_);
		tmpPtrsVector.push_back(eleAHist_modTrkIso_);
		tmpPtrsVector.push_back(eleAHist_dRleq03modTrkIso_);
		tmpPtrsVector.push_back(eleAHist_dRleq03modTrkIsoNoCtfTrk_);
		tmpPtrsVector.push_back(eleAHist_modEmHad1Iso_);
		tmpPtrsVector.push_back(eleAHist_modEmHad1IsoOverCutThr_);
		tmpPtrsVector.push_back(eleAHist_modHad1IsoOverModEmHad1_);
		tmpPtrsVector.push_back(eleAHist_Had1Iso_);
		tmpPtrsVector.push_back(eleAHist_modHad1Iso_);
		tmpPtrsVector.push_back(eleAHist_modHad1IsoOverHad1_);

		tmpPtrsVector.push_back(eleBHist_et_);
		tmpPtrsVector.push_back(eleBHist_energy_);
		tmpPtrsVector.push_back(eleBHist_eta_);
		tmpPtrsVector.push_back(eleBHist_phi_);
		tmpPtrsVector.push_back(eleBHist_modTrkIso_);
		tmpPtrsVector.push_back(eleBHist_dRleq03modTrkIso_);
		tmpPtrsVector.push_back(eleBHist_dRleq03modTrkIsoNoCtfTrk_);
		tmpPtrsVector.push_back(eleBHist_modEmHad1Iso_);
		tmpPtrsVector.push_back(eleBHist_modEmHad1IsoOverCutThr_);
		tmpPtrsVector.push_back(eleBHist_modHad1IsoOverModEmHad1_);
		tmpPtrsVector.push_back(eleBHist_Had1Iso_);
		tmpPtrsVector.push_back(eleBHist_modHad1Iso_);
		tmpPtrsVector.push_back(eleBHist_modHad1IsoOverHad1_);

		return tmpPtrsVector;
	}
}

#endif
