#ifndef tswHEEPVarDistns_h
#define tswHEEPVarDistns_h

#include "TH2.h"

namespace tsw{
	class HEEPVarDistns{
		public:
			HEEPVarDistns(TString histNamePrefix, const TString titleString_eleType, const TString titleString_cutsPhrase);
			~HEEPVarDistns();

			//Fill methods for the histograms ...
			void FillHistos(tsw::HEEPEle& , Double_t eleWeight);
			void FillNumberElesHist(unsigned int ele_number, const Double_t& evtWeight);
			void FillChargeHist(Int_t ele_charge, const Double_t& evtWeight);
			void FillEtHists(Double_t ele_Et, Double_t ele_HEEP_Et, const Double_t& evtWeight);
			void FillEtaHists(Double_t ele_Eta, Double_t ele_scEta, const Double_t& evtWeight);
			void FillEcalDrivenHists(bool ele_ecalDriven, bool ele_ecalDrivenSeed, const Double_t& evtWeight);
			void FilldEtadPhiInHists(Double_t ele_HEEP_dEtaIn, Double_t ele_HEEP_dPhiIn, const Double_t& evtWeight);
			void FillHoverEHist(Double_t HoverE, const Double_t& evtWeight);
			void FillSigmaIetaIetaHists(Double_t ele_sigmaIetaIeta, Double_t scSigmaIetaIeta, const Double_t& evtWeight);
			void FillIsoHists(Double_t ele_dr03EmIsoEt, Double_t ele_dr03HadDepth1IsoEt, Double_t ele_dr03HadDepth2IsoEt, Double_t ele_dr03TkIsoPt, const Double_t& evtWeight);
			void FilleAxBHists(Double_t ele_e2x5Max, Double_t ele_e5x5, const Double_t& evtWeight);

			//"Write" method for the histos ...
			void WriteHistos(TFile* ptr_file);

			//Method returning vector of pointers to all of the histograms (for error calc'n etc) ...
			TH1D* GetPointerToNumberElesHist();
			std::vector<TH1D*> GetPointersToHistos();

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);

		private:
			//Method that deletes all of the histos ...
			void DeleteHistos();

		private:
			TString hNamePrefix_;
			TH1D* eleHist_et_;
			TH1D* eleHist_gsfEt_;
			TH1D* eleHist_scEt_;
			TH1D* eleHist_energy_;
			//	float gsfEnergy(){return eleStr_.gsfEnergy_;}
			//	float caloEnergy(){return eleStr_.caloEnergy_;}
			TH1D* eleHist_eta_;
			TH1D* eleHist_scEta_;
			//	float detEta(){return eleStr_.detEta_;}
			//	float detEtaAbs(){return eleStr_.detEtaAbs_;}
			TH1D* eleHist_phi_;
			TH1D* eleHist_scPhi_;
			//	float detPhi(){return eleStr_.detPhi_;}
			//TH1D* eleHist_zVtx_;
			//	ROOT::Math::XYZTVector p4(){return eleStr_.p4_;}
			//	ROOT::Math::XYZTVector gsfP4(){return eleStr_.gsfP4_;}
			TH1D* eleHist_caloVsRecHitsEnergyRatio_;
			TH1D* eleHist_caloEnergyFracError_;
			TH2D* ele2DHist_caloRecHitEnergyRatioVsNumRecHits_;

			// 'Classification'
			//	int classification(){return eleStr_.classification_;}
			TH1D* eleHist_isEcalDriven_;
			//	bool isTrackerDriven(){return eleStr_.isTrackerDriven_;}
			TH1D* eleHist_isEB_;
			TH1D* eleHist_isEE_;

			// Track variables ...
			TH1D* eleHist_charge_;
			//	int trkCharge(){return eleStr_.trkCharge_;}
			//	float pVtx(){return eleStr_.pVtx_;}
			//	float pCalo(){return eleStr_.pCalo_;}
			//	float ptVtx(){return eleStr_.ptVtx_;}
			//	float ptCalo(){return eleStr_.ptCalo_;}

			// Various other variables ...
			TH1D* eleHist_hOverE_;
			TH1D* eleHist_hOverEnonzero_;
			TH1D* eleHist_dEtaIn_;
			TH1D* eleHist_dPhiIn_;
			//	float dPhiOut(){return eleStr_.dPhiOut_;}
			TH1D* eleHist_epIn_;
			TH1D* eleHist_epOut_;
			//	float fbrem(){return eleStr_.fbrem_;}
			//	float bremFrac(){return eleStr_.bremFrac_;}
			TH1D* eleHist_invEOverInvP_;

			// Shower shape variables
			//	float sigmaEtaEta(){return eleStr_.sigmaEtaEta_;}
			//	float sigmaEtaEtaUnCorr(){return eleStr_.sigmaEtaEtaUnCorr_;}
			TH1D* eleHist_sigmaIEtaIEta_;
			//	float e1x5(){return eleStr_.e1x5_;}
			// float e2x5Max(){return eleStr_.e2x5Max_;}
			// float e5x5(){return eleStr_.e5x5_;}
			// float e1x5Over5x5(){return eleStr_.e1x5Over5x5_;}
			// float e2x5MaxOver5x5(){return eleStr_.e2x5MaxOver5x5_;}

			// Isolation variables ...
			TH1D* eleHist_isolEm_;
			//	float isolHad(){return eleStr_.isolHad_;}
			TH1D* eleHist_isolHadDepth1_;
			TH1D* eleHist_isolHadDepth2_;
			TH1D* eleHist_isolPtTrks_;
			//	float isolEmHadDepth1(){return eleStr_.isolEmHadDepth1_;}


			TH1D* hist_ele_number_;
			/*TH1D* hist_ele_charge_;
			TH1D* hist_ele_Et_;
			TH1D* hist_ele_heepEt_;
			TH1D* hist_ele_Eta_;
			TH1D* hist_ele_scEta_;
			TH1D* hist_ele_ecalDriven_;
			TH1D* hist_ele_ecalDrivenSeed_;
			TH1D* hist_ele_heepdEtaIn_;
			TH1D* hist_ele_heepdPhiIn_;
			TH1D* hist_ele_HoverE_;
			TH1D* hist_ele_sigmaIetaIeta_;
			TH1D* hist_ele_scSigmaIetaIeta_;
			TH1D* hist_ele_dr03EmIsoEt_;
			TH1D* hist_ele_dr03Had1IsoEt_;
			TH1D* hist_ele_dr03Had2IsoEt_;
			TH1D* hist_ele_dr03TkIsoPt_;
			TH1D* hist_ele_e2x5Max_;
			TH1D* hist_ele_e5x5_;*/
	};

	HEEPVarDistns::HEEPVarDistns(TString hNamePrefix, const TString str_eleType, const TString str_cutsPhrase):
		hNamePrefix_(hNamePrefix)
	{
		//Initialise the histograms here ...
		TString str_eleType_GSFeles_cutsPhrase = str_eleType + " GSF electrons " + str_cutsPhrase;
		TString str_NoOfEleTypeEles = "Number of " + str_eleType + " GSF electrons";

		hist_ele_number_ = new TH1D(hNamePrefix + "number", "Histogram of number of " + str_eleType_GSFeles_cutsPhrase + "in each event; Number of " + str_eleType + " GSF electrons in event; Number of events", 6, -0.5, 5.5);
		/*hist_ele_charge_ = new TH1D(hNamePrefix + "charge", "Charge distribution of the " + str_eleType_GSFeles_cutsPhrase + "; charge /e; " + str_NoOfEleTypeEles, 3, -1.5, 1.5);

		hist_ele_Et_     = new TH1D(hNamePrefix + "Et",     "E_{T} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; electron E_{T} /GeV; Number of " + str_eleType + " GSF electrons per ...",           25, 0.0, 100.0);
		hist_ele_heepEt_ = new TH1D(hNamePrefix + "heepEt", "E_{T,HEEP} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; electron E_{T,HEEP} /GeV; Number of " + str_eleType + " GSF electrons per ...", 25, 0.0, 100.0);

		hist_ele_Eta_   = new TH1D(hNamePrefix + "Eta",   "#eta distribution of the " + str_eleType_GSFeles_cutsPhrase + "; electron pseudorapidity, #eta; " + str_NoOfEleTypeEles + " per ...",               30, -3.0, +3.0);
		hist_ele_scEta_ = new TH1D(hNamePrefix + "scEta", "#eta^{SC} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; supercluster pseudorapidity, 'eta^{SC}; " + str_NoOfEleTypeEles + " per ...", 30, -3.0, +3.0);

		hist_ele_ecalDriven_     = new TH1D(hNamePrefix + "ecalDriven",     "Histogram of whether " + str_eleType_GSFeles_cutsPhrase + "are ECAL driven or not; ecalDriven; " + str_NoOfEleTypeEles, 2, -0.5, 1.5);
		hist_ele_ecalDrivenSeed_ = new TH1D(hNamePrefix + "ecalDrivenSeed", "Histogram of whether " + str_eleType_GSFeles_cutsPhrase + "are ECAL driven or not; ecalDrivenSeed; " + str_NoOfEleTypeEles, 2, -0.5, 1.5);

		hist_ele_heepdEtaIn_ = new TH1D(hNamePrefix + "heepdEtaIn", "#Delta#eta_{In,HEEP} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Delta#eta_{In,HEEP}; " + str_NoOfEleTypeEles + " per ...", 25, -0.02, +0.3);
		hist_ele_heepdPhiIn_ = new TH1D(hNamePrefix + "heepdPhiIn", "#Delta#phi_{In,HEEP} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Delta#phi_{In,HEEP}; " + str_NoOfEleTypeEles + " per ...", 30, -0.15, +0.15);

		hist_ele_HoverE_ = new TH1D(hNamePrefix + "HoverE", "H/E distribution of the " + str_eleType_GSFeles_cutsPhrase + "; H/E; " + str_NoOfEleTypeEles + " per ...", 30, 0.0, 0.15);

		hist_ele_sigmaIetaIeta_   = new TH1D(hNamePrefix + "sigmaIetaIeta",   "#sigma_{i#eta i#eta} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #sigma_{i#eta i#eta};" + str_NoOfEleTypeEles + " per ...",           25, 0.0, 0.05);
		hist_ele_scSigmaIetaIeta_ = new TH1D(hNamePrefix + "scSigmaIetaIeta", "#sigma_{i#eta i#eta}^{SC} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #sigma_{i#eta i#eta}^{SC};" + str_NoOfEleTypeEles + " per ...", 25, 0.0, 0.05);

		hist_ele_dr03EmIsoEt_   = new TH1D(hNamePrefix + "dr03EmIsoEt",   "#Sigma{}E_{T,Em}, #Delta{}R<0.3 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Sigma{}E_{T} in ECAL within #Delta{}R=0.3 of electron /GeV; " + str_NoOfEleTypeEles + " per ...",            50, 0.0, 10.0);
		hist_ele_dr03Had1IsoEt_ = new TH1D(hNamePrefix + "dr03Had1IsoPt", "#Sigma{}E_{T,Had1}, #Delta{}R<0.3 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Sigma{}E_{T} in HCAL (depth1) within #Delta{}R=0.3 of electron /GeV; " + str_NoOfEleTypeEles + " per ...", 50, 0.0, 1.0);
		hist_ele_dr03Had2IsoEt_ = new TH1D(hNamePrefix + "dr03Had2IsoPt", "#Sigma{}E_{T,Had2}, #Delta{}R<0.3 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Sigma{}E_{T} in HCAL (depth2) within #Delta{}R=0.3 of electron /GeV; " + str_NoOfEleTypeEles + " per ...", 50, 0.0, 1.0);
		hist_ele_dr03TkIsoPt_   = new TH1D(hNamePrefix + "dr03TkIsoPt",   "#Sigma{}p_{T}, #Delta{}R<0.3 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Sigma{}p_{T} of other tracks within #Delta{}R=0.3 of electron /GeVc^{-1}; " + str_NoOfEleTypeEles + " per ...", 25, 0.0, 5.0);

		hist_ele_e2x5Max_ = new TH1D(hNamePrefix + "e2x5Max", "E^{2x5Max} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E^{2x5Max} /GeV; " + str_NoOfEleTypeEles + " per ...", 50, 0.0, 200.0);
		hist_ele_e5x5_    = new TH1D(hNamePrefix + "e5x5",    "E^{5x5} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E^{5x5} /GeV; " + str_NoOfEleTypeEles + " per ...", 50, 0.0, 200.0);*/

		Double_t hMin_et = 0.0; Double_t hMax_et = 100.0; Int_t hNBins_et = 50;
		eleHist_et_ = new TH1D(hNamePrefix + "et",         "E_{T} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E_{T} /GeV; " + str_NoOfEleTypeEles + " per ...",           hNBins_et, hMin_et, hMax_et);
		eleHist_gsfEt_ = new TH1D(hNamePrefix + "gsfEt",   "E_{T, GSF} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E_{T, GSF} /GeV; " + str_NoOfEleTypeEles + " per ...", hNBins_et, hMin_et, hMax_et);
		eleHist_scEt_ = new TH1D(hNamePrefix + "scEt",     "E_{T, SC} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E_{T, SC} /GeV; " + str_NoOfEleTypeEles + " per ...",   hNBins_et, hMin_et, hMax_et);
		eleHist_energy_ = new TH1D(hNamePrefix + "energy", "Energy distribution of the " + str_eleType_GSFeles_cutsPhrase + "; Energy /GeV; " + str_NoOfEleTypeEles + " per ...",         hNBins_et, hMin_et, hMax_et);
		Double_t hMin_eta = -3.0; Double_t hMax_eta = +3.0; Int_t hNBins_eta = 60;
		eleHist_eta_   = new TH1D(hNamePrefix + "eta", "#eta distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #eta; " + str_NoOfEleTypeEles + " per ...", hNBins_eta, hMin_eta, hMax_eta);
		eleHist_scEta_ = new TH1D(hNamePrefix + "scEta", "#eta_{SC} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #eta_{SC}; " + str_NoOfEleTypeEles + " per ...", hNBins_eta, hMin_eta, hMax_eta);
		Double_t hMin_phi = -3.1416; Double_t hMax_phi = +3.1416; Int_t hNBins_phi = 60;
		eleHist_phi_   = new TH1D(hNamePrefix + "phi", "#phi distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #phi; " + str_NoOfEleTypeEles + " per ...", hNBins_phi, hMin_phi, hMax_phi);
		eleHist_scPhi_ = new TH1D(hNamePrefix + "scPhi", "#phi_{SC} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #phi_{SC}; " + str_NoOfEleTypeEles + " per ...", hNBins_phi, hMin_phi, hMax_phi);


		eleHist_caloVsRecHitsEnergyRatio_ = new TH1D(hNamePrefix + "caloVsRecHitsEnergyRatio", "E^{calo}/E^{SC}_{recHits} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E^{calo}/E^{SC}_{recHits}; " + str_NoOfEleTypeEles + " per ...", 200, 0.95, 1.15);
		eleHist_caloEnergyFracError_ = new TH1D(hNamePrefix + "caloEnergyFracError", "E_{ecalError}/E_{caloEnergy} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E_{ecalError}/E_{caloEnergy}; " + str_NoOfEleTypeEles + " per ...", 200, 0.95, 1.15);
		ele2DHist_caloRecHitEnergyRatioVsNumRecHits_ = new TH2D(hNamePrefix+"caloRecHitEnergyRatioVsNumRecHits", "E^{calo}/E^{SC}_{recHits} vs N^{SC}_{recHits} 2-D plot of the " + str_eleType_GSFeles_cutsPhrase + "; E^{calo}/E^{SC}_{recHits}; N^{SC}_{recHits};" + str_NoOfEleTypeEles + " per ...", 40, 0.95, 1.15, 60, -0.5, 119.5);
		//eleHist_zVtx_;

		// 'Classification'
		Double_t hMin_bool = -0.5; Double_t hMax_bool = +1.5; Int_t hNBins_bool = 2;
		eleHist_isEcalDriven_ = new TH1D(hNamePrefix + "isEcalDriven", "isEcalDriven distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isEcalDriven; " + str_NoOfEleTypeEles + " per ...", hNBins_bool, hMin_bool, hMax_bool);
		eleHist_isEB_         = new TH1D(hNamePrefix + "isEB",         "isEB distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isEB; " + str_NoOfEleTypeEles + " per ...", hNBins_bool, hMin_bool, hMax_bool);
		eleHist_isEE_         = new TH1D(hNamePrefix + "isEE",         "isEE distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isEE; " + str_NoOfEleTypeEles + " per ...", hNBins_bool, hMin_bool, hMax_bool);

		// Track variables ...
		Double_t hMin_charge = -1.5; Double_t hMax_charge = +1.5; Int_t hNBins_charge = 3;
		eleHist_charge_ = new TH1D(hNamePrefix + "charge", "Charge distribution of the " + str_eleType_GSFeles_cutsPhrase + "; charge /e; " + str_NoOfEleTypeEles + " per ...", hNBins_charge, hMin_charge, hMax_charge);

		// Various other variables ...
		Double_t hMin_hOverE = 0.0; Double_t hMax_hOverE = 0.3; Int_t hNBins_hOverE = 60;
		eleHist_hOverE_ = new TH1D(hNamePrefix + "hOverE", "H/0 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; H/E; " + str_NoOfEleTypeEles + " per ...", hNBins_hOverE, hMin_hOverE, hMax_hOverE);
		eleHist_hOverEnonzero_ = new TH1D(hNamePrefix + "hOverEnonzero", "H/E (!=0.0) distribution of the " + str_eleType_GSFeles_cutsPhrase + "; H/E; " + str_NoOfEleTypeEles + " per ...", hNBins_hOverE, hMin_hOverE, hMax_hOverE);
		Double_t hMin_dEtaIn = -0.15; Double_t hMax_dEtaIn = +0.15; Int_t hNBins_dEtaIn = 60;
		eleHist_dEtaIn_ = new TH1D(hNamePrefix + "dEtaIn", "#Delta#eta_{In} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Delta#eta_{In}; " + str_NoOfEleTypeEles + " per ...", hNBins_dEtaIn, hMin_dEtaIn, hMax_dEtaIn);
		Double_t hMin_dPhiIn = -0.3; Double_t hMax_dPhiIn = +0.3; Int_t hNBins_dPhiIn = 60;
		eleHist_dPhiIn_ = new TH1D(hNamePrefix + "dPhiIn", "#Delta#phi_{In} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #Delta#phi_{In}; " + str_NoOfEleTypeEles + " per ...", hNBins_dPhiIn, hMin_dPhiIn, hMax_dPhiIn);
		Double_t hMin_eOverP = 0.0; Double_t hMax_eOverP = 25.0; Int_t hNBins_eOverP = 50;
		eleHist_epIn_ = new TH1D(hNamePrefix + "epIn", "E/p_{In} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E/p_{In}; " + str_NoOfEleTypeEles + " per ...", hNBins_eOverP, hMin_eOverP, hMax_eOverP);
		eleHist_epOut_ = new TH1D(hNamePrefix + "epOut", "E/p_{Out} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E/p_{Out}; " + str_NoOfEleTypeEles + " per ...", hNBins_eOverP, hMin_eOverP, hMax_eOverP);
		Double_t hMin_invEOverP = -0.5; Double_t hMax_invEOverP = +0.5; Int_t hNBins_invEOverP = 50;
		eleHist_invEOverInvP_ = new TH1D(hNamePrefix + "invEOverInvP", "E^{-1}/p^{-1} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; E^{-1}/p^{-1}; " + str_NoOfEleTypeEles + " per ...", hNBins_invEOverP, hMin_invEOverP, hMax_invEOverP);

		// Shower shape variables
		Double_t hMin_sIEtaIEta = 0.0; Double_t hMax_sIEtaIEta = 0.07; Int_t hNBins_sIEtaIEta = 70;
		eleHist_sigmaIEtaIEta_ = new TH1D(hNamePrefix + "sigmaIEtaIEta", "#sigma_{i#eta i#eta} distribution of the " + str_eleType_GSFeles_cutsPhrase + "; #sigma_{i#eta i#eta}; " + str_NoOfEleTypeEles + " per ...", hNBins_sIEtaIEta, hMin_sIEtaIEta, hMax_sIEtaIEta);

		// Isolation variables ...
		Double_t hMin_isolEm = 0.0; Double_t hMax_isolEm = 70.0; Int_t hNBins_isolEm = 70;
		eleHist_isolEm_ = new TH1D(hNamePrefix + "isolEm", "isolEm distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isolEm; " + str_NoOfEleTypeEles + " per ...", hNBins_isolEm, hMin_isolEm, hMax_isolEm);
		Double_t hMin_isolHad1 = 0.0; Double_t hMax_isolHad1 = 30.0; Int_t hNBins_isolHad1 = 60;
		eleHist_isolHadDepth1_ = new TH1D(hNamePrefix + "isolHadDepth1", "isolHadDepth1 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isolHadDepth1; " + str_NoOfEleTypeEles + " per ...", hNBins_isolHad1, hMin_isolHad1, hMax_isolHad1);
		Double_t hMin_isolHad2 = 0.0; Double_t hMax_isolHad2 = 10.0; Int_t hNBins_isolHad2 = 50;
		eleHist_isolHadDepth2_ = new TH1D(hNamePrefix + "isolHadDepth2", "isolHadDepth2 distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isolHadDepth2; " + str_NoOfEleTypeEles + " per ...", hNBins_isolHad2, hMin_isolHad2, hMax_isolHad2);
		Double_t hMin_isolTrks = 0.0; Double_t hMax_isolTrks = 100.0; Int_t hNBins_isolTrks = 50;
		eleHist_isolPtTrks_ = new TH1D(hNamePrefix + "isolPtTrks", "isolPtTrks distribution of the " + str_eleType_GSFeles_cutsPhrase + "; isolPtTrks; " + str_NoOfEleTypeEles + " per ...", hNBins_isolTrks, hMin_isolTrks, hMax_isolTrks);

		// Set-up histograms so that errors are automatically calculated as the histos are filled
		std::vector<TH1D*> ptrsToHists = GetPointersToHistos();
		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->Sumw2();
			ptrsToHists.at(iHist)->SetDirectory(0);
		}
	}

	HEEPVarDistns::~HEEPVarDistns(){
//		DeleteHistos();
	}

	void HEEPVarDistns::FillHistos(tsw::HEEPEle& ele, Double_t eleWeight){

		eleHist_et_->Fill(    ele.et()    , eleWeight);
		eleHist_gsfEt_->Fill( ele.gsfEt() , eleWeight);
		eleHist_scEt_->Fill(  ele.scEt()  , eleWeight);
		eleHist_energy_->Fill(ele.energy(), eleWeight);
		eleHist_eta_->Fill(   ele.eta()   , eleWeight);
		eleHist_scEta_->Fill( ele.scEta() , eleWeight);
		eleHist_phi_->Fill(   ele.phi()   , eleWeight);
		eleHist_scPhi_->Fill( ele.scPhi() , eleWeight);
		//eleHist_zVtx_->Fill( ele.zVtx() , eleWeight);
		eleHist_caloVsRecHitsEnergyRatio_->Fill( ele.caloEnergy()/ele.SC_totEnergyRecHits(), eleWeight );
		ele2DHist_caloRecHitEnergyRatioVsNumRecHits_->Fill( ele.caloEnergy()/ele.SC_totEnergyRecHits(), ele.SC_totNumRecHits(), eleWeight );
		eleHist_caloEnergyFracError_->Fill( ele.caloEnergy()/ele.ecalEnergyError(), eleWeight);

		// 'Classification'
		eleHist_isEcalDriven_->Fill( ele.isEcalDriven() , eleWeight);
		eleHist_isEB_->Fill( ele.isEB(), eleWeight);
		eleHist_isEE_->Fill( ele.isEE(), eleWeight);

		// Track variables ...
		eleHist_charge_->Fill( ele.charge(), eleWeight);

		// Various other variables ...
		eleHist_hOverE_->Fill( ele.hOverE(), eleWeight);
		if( ele.hOverE()!=0.0 )
			eleHist_hOverEnonzero_->Fill( ele.hOverE(), eleWeight);
		eleHist_dEtaIn_->Fill( ele.dEtaIn(), eleWeight);
		eleHist_dPhiIn_->Fill( ele.dPhiIn(), eleWeight);
		eleHist_epIn_->Fill( ele.epIn(), eleWeight);
		eleHist_epOut_->Fill( ele.epOut(), eleWeight);
		eleHist_invEOverInvP_->Fill( ele.invEOverInvP(), eleWeight);

		// Shower shape variables
		/*eleHist_sigmaIEtaIEta_->Fill( ele.sigmaIEtaIEta(), eleWeight);

		// Isolation variables ...
		eleHist_isolEm_->Fill( ele.isolEm(), eleWeight);
		eleHist_isolHadDepth1_->Fill( ele.isolHadDepth1(), eleWeight);
		eleHist_isolHadDepth2_->Fill( ele.isolHadDepth2(), eleWeight);
		eleHist_isolPtTrks_->Fill( ele.isolPtTrks(), eleWeight);*/
	}

	void HEEPVarDistns::FillNumberElesHist(unsigned int ele_number, const Double_t& evtWeight){
		hist_ele_number_->Fill(ele_number, evtWeight);
	}
	void HEEPVarDistns::FillChargeHist(Int_t ele_charge, const Double_t& evtWeight){
		//hist_ele_charge_->Fill(ele_charge, evtWeight);
	}
	void HEEPVarDistns::FillEtHists(Double_t ele_Et, Double_t ele_HEEP_Et, const Double_t& evtWeight){
		//hist_ele_Et_->Fill(ele_Et, evtWeight);
		//hist_ele_heepEt_->Fill(ele_HEEP_Et, evtWeight);
	}
	void HEEPVarDistns::FillEtaHists(Double_t ele_Eta, Double_t ele_scEta, const Double_t& evtWeight){
		//hist_ele_Eta_->Fill(ele_Eta, evtWeight);
		//hist_ele_scEta_->Fill(ele_scEta, evtWeight);
	}
	void HEEPVarDistns::FillEcalDrivenHists(bool ele_ecalDriven, bool ele_ecalDrivenSeed, const Double_t& evtWeight){
		//hist_ele_ecalDriven_->Fill(ele_ecalDriven, evtWeight);
		//hist_ele_ecalDrivenSeed_->Fill(ele_ecalDrivenSeed, evtWeight);
	}
	void HEEPVarDistns::FilldEtadPhiInHists(Double_t ele_HEEP_dEtaIn, Double_t ele_HEEP_dPhiIn, const Double_t& evtWeight){
		//hist_ele_heepdEtaIn_->Fill(ele_HEEP_dEtaIn, evtWeight);
		//hist_ele_heepdPhiIn_->Fill(ele_HEEP_dPhiIn, evtWeight);
	}
	void HEEPVarDistns::FillHoverEHist(Double_t ele_HoverE, const Double_t& evtWeight){
		//hist_ele_HoverE_->Fill(ele_HoverE, evtWeight);
	}
	void HEEPVarDistns::FillSigmaIetaIetaHists(Double_t ele_sigmaIetaIeta, Double_t ele_scSigmaIetaIeta, const Double_t& evtWeight){
		//hist_ele_sigmaIetaIeta_->Fill(ele_sigmaIetaIeta, evtWeight);
		//hist_ele_scSigmaIetaIeta_->Fill(ele_scSigmaIetaIeta, evtWeight);
	}
	void HEEPVarDistns::FillIsoHists(Double_t ele_dr03EmIsoEt, Double_t ele_dr03HadDepth1IsoEt, Double_t ele_dr03HadDepth2IsoEt, Double_t ele_dr03TkIsoPt, const Double_t& evtWeight){
		//hist_ele_dr03EmIsoEt_->Fill(ele_dr03EmIsoEt, evtWeight);
		//hist_ele_dr03Had1IsoEt_->Fill(ele_dr03HadDepth1IsoEt, evtWeight);
		//hist_ele_dr03Had2IsoEt_->Fill(ele_dr03HadDepth2IsoEt, evtWeight);
		//hist_ele_dr03TkIsoPt_->Fill(ele_dr03TkIsoPt, evtWeight);
	}
	void HEEPVarDistns::FilleAxBHists(Double_t ele_e2x5Max, Double_t ele_e5x5, const Double_t& evtWeight){
		//hist_ele_e2x5Max_->Fill(ele_e2x5Max, evtWeight);
		//hist_ele_e5x5_->Fill(ele_e5x5, evtWeight);
	}

	void HEEPVarDistns::WriteHistos(TFile* ptr_file)
	{
		ptr_file->mkdir(hNamePrefix_);
		ptr_file->cd(hNamePrefix_);

		hist_ele_number_->Write();

		eleHist_et_->Write();
		eleHist_gsfEt_->Write();
		eleHist_scEt_->Write();
		eleHist_energy_->Write();
		eleHist_eta_->Write();
		eleHist_scEta_->Write();
		eleHist_phi_->Write();
		eleHist_scPhi_->Write();
		//eleHist_zVtx_->Write();
		eleHist_caloVsRecHitsEnergyRatio_->Write();
		eleHist_caloEnergyFracError_->Write();
		ele2DHist_caloRecHitEnergyRatioVsNumRecHits_->Write();

		// 'Classification'
		eleHist_isEcalDriven_->Write();
		eleHist_isEB_->Write();
		eleHist_isEE_->Write();

		// Track variables ...
		eleHist_charge_->Write();

		// Various other variables ...
		eleHist_hOverE_->Write();
		eleHist_hOverEnonzero_->Write();
		eleHist_dEtaIn_->Write();
		eleHist_dPhiIn_->Write();
		eleHist_epIn_->Write();
		eleHist_epOut_->Write();
		eleHist_invEOverInvP_->Write();

		/*hist_ele_charge_->Write();
		hist_ele_Et_->Write();
		hist_ele_heepEt_->Write();
		hist_ele_Eta_->Write();
		hist_ele_scEta_->Write();
		hist_ele_ecalDriven_->Write();
		hist_ele_ecalDrivenSeed_->Write();
		hist_ele_heepdEtaIn_->Write();
		hist_ele_heepdPhiIn_->Write();
		hist_ele_HoverE_->Write();
		hist_ele_sigmaIetaIeta_->Write();
		hist_ele_scSigmaIetaIeta_->Write();
		hist_ele_dr03EmIsoEt_->Write();
		hist_ele_dr03Had1IsoEt_->Write();
		hist_ele_dr03Had2IsoEt_->Write();
		hist_ele_dr03TkIsoPt_->Write();
		hist_ele_e2x5Max_->Write();
		hist_ele_e5x5_->Write();*/

		ptr_file->cd("..");
	}

	TH1D* HEEPVarDistns::GetPointerToNumberElesHist(){
		return hist_ele_number_;
	}

	std::vector<TH1D*> HEEPVarDistns::GetPointersToHistos(){
		//Currenlty a placeholder ...
		std::vector<TH1D*> tmpVector;
		tmpVector.clear();

		tmpVector.push_back(eleHist_et_);
		tmpVector.push_back(eleHist_gsfEt_);
		tmpVector.push_back(eleHist_scEt_);
		tmpVector.push_back(eleHist_energy_);
		tmpVector.push_back(eleHist_eta_);
		tmpVector.push_back(eleHist_scEta_);
		tmpVector.push_back(eleHist_phi_);
		tmpVector.push_back(eleHist_scPhi_);
		tmpVector.push_back(eleHist_caloVsRecHitsEnergyRatio_);
		tmpVector.push_back(eleHist_caloEnergyFracError_);
		//tmpVector.push_back(eleHist_zVtx_);

		// 'Classification'
		tmpVector.push_back(eleHist_isEcalDriven_);
		tmpVector.push_back(eleHist_isEB_);
		tmpVector.push_back(eleHist_isEE_);

		// Track variables ...
		tmpVector.push_back(eleHist_charge_);

		// Various other variables ...
		tmpVector.push_back(eleHist_hOverE_);
		tmpVector.push_back(eleHist_hOverEnonzero_);
		tmpVector.push_back(eleHist_dEtaIn_);
		tmpVector.push_back(eleHist_dPhiIn_);
		tmpVector.push_back(eleHist_epIn_);
		tmpVector.push_back(eleHist_epOut_);
		tmpVector.push_back(eleHist_invEOverInvP_);

		// Shower shape variables
		tmpVector.push_back(eleHist_sigmaIEtaIEta_);

		// Isolation variables ...
		tmpVector.push_back(eleHist_isolEm_);
		tmpVector.push_back(eleHist_isolHadDepth1_);
		tmpVector.push_back(eleHist_isolHadDepth2_);
		tmpVector.push_back(eleHist_isolPtTrks_);


		return tmpVector;
	}

	void HEEPVarDistns::SetHistAttributes(Int_t lineCol, Int_t lineStyle){
		std::vector<TH1D*> ptrsToHists;
		ptrsToHists = GetPointersToHistos();

		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->SetLineColor(lineCol);
			ptrsToHists.at(iHist)->SetLineStyle(lineStyle);
		}
	}

	void HEEPVarDistns::DeleteHistos(){
		//Currently a placeholder ...
		delete hist_ele_number_;
		/*delete hist_ele_charge_;
		delete hist_ele_Et_;
		delete hist_ele_heepEt_;
		delete hist_ele_Eta_;
		delete hist_ele_scEta_;
		delete hist_ele_ecalDriven_;
		delete hist_ele_ecalDrivenSeed_;
		delete hist_ele_heepdEtaIn_;
		delete hist_ele_heepdPhiIn_;
		delete hist_ele_HoverE_;
		delete hist_ele_sigmaIetaIeta_;
		delete hist_ele_scSigmaIetaIeta_;
		delete hist_ele_dr03EmIsoEt_;
		delete hist_ele_dr03Had1IsoEt_;
		delete hist_ele_dr03Had2IsoEt_;
		delete hist_ele_dr03TkIsoPt_;
		delete hist_ele_e2x5Max_;
		delete hist_ele_e5x5_;*/
		delete eleHist_et_;
		delete eleHist_gsfEt_;
		delete eleHist_scEt_;
		delete eleHist_energy_;
		delete eleHist_eta_;
		delete eleHist_scEta_;
		delete eleHist_phi_;
		delete eleHist_scPhi_;
		delete eleHist_caloVsRecHitsEnergyRatio_;
		delete eleHist_caloEnergyFracError_;
		delete ele2DHist_caloRecHitEnergyRatioVsNumRecHits_;

		//delete eleHist_zVtx_;

		// 'Classification'
		delete eleHist_isEcalDriven_;
		delete eleHist_isEB_;
		delete eleHist_isEE_;

		// Track variables ...
		delete eleHist_charge_;

		// Various other variables ...
		delete eleHist_hOverE_;
		delete eleHist_hOverEnonzero_;
		delete eleHist_dEtaIn_;
		delete eleHist_dPhiIn_;
		delete eleHist_epIn_;
		delete eleHist_epOut_;
		delete eleHist_invEOverInvP_;

		// Shower shape variables
		delete eleHist_sigmaIEtaIEta_;

		// Isolation variables ...
		delete eleHist_isolEm_;
		delete eleHist_isolHadDepth1_;
		delete eleHist_isolHadDepth2_;
		delete eleHist_isolPtTrks_;
	}
}

#endif
