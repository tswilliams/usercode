#ifndef tswMuonDistns_h
#define tswMuonDistns_h

namespace tsw{
	class MuonDistns{
		public:
			//CTORs and DTORs
			MuonDistns(TString hNamePrefix, const TString titleStr_MuonType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Double_t hMax_mass, Int_t hNBins_pt, Double_t hMax_pt);
			~MuonDistns();

			//Method for filling the histograms
			void FillHistos(tsw::Muon , const Double_t );

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);
			//Method to write the histograms to file ...
			void WriteHistos(TFile* ptr_outFile);

			//Method to get vector of pointers to all histograms ...
			std::vector<TH1D*> GetPtrsToHistos();
		private:
			TString hNamePrefix_;
			TH1D* muH_region_;
			TH1D* muH_charge_;
			TH1D* muH_p_;
			TH1D* muH_pT_;
			TH1D* muH_eta_;
			TH1D* muH_phi_;
			TH1D* muH_isGlobMuon_;
			TH1D* muH_isTrkrMuon_;
			TH1D* muH_isStandAloneMuon_;
			TH1D* muH_nMatchedMuStns_;
			TH1D* muH_isolR03_sumPt_;
			TH1D* muH_glob_nMuHits_;
			TH1D* muH_glob_normdChi2_;
			TH1D* muH_inner_nPixHits_;
			TH1D* muH_inner_nTrkrHits_;
			TH1D* muH_inner_dxyOrigin_;
	};

	//--------------------------------//
	//---- CTOR ...
	MuonDistns::MuonDistns(TString hNamePrefix, const TString titleStr_MuonType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Double_t hMax_mass, Int_t hNBins_pt, Double_t hMax_pt) :
		hNamePrefix_(hNamePrefix)
	{
		//Form standard phrases ...
		TString str_muonType_muons_cutsPhrase = titleStr_MuonType + " muons " + titleStr_cutsPhrase;
		TString str_NoOfMuons = "Number of muons";

		//Initialise the histograms ...
		muH_region_ = new TH1D(hNamePrefix + "region", "Muon 'region' distribution (" + str_muonType_muons_cutsPhrase + ") ; Region;" + str_NoOfMuons, 4, -0.5, 3.5);
		muH_charge_ = new TH1D(hNamePrefix + "charge", "Muon charge distribution (" + str_muonType_muons_cutsPhrase + "); Charge /e;" + str_NoOfMuons, 3, -1.5, 1.5);
		muH_p_      = new TH1D(hNamePrefix + "p",      "Muon momentum distribution (" + str_muonType_muons_cutsPhrase + "); Momentum, p /GeVc^{-1};" + str_NoOfMuons, hNBins_pt, 0.0, hMax_pt);
		muH_pT_     = new TH1D(hNamePrefix + "pT",     "Muon p_{T} distribution (" + str_muonType_muons_cutsPhrase + "); Transverse momentum, p_{T} /GeVc^{-1};" + str_NoOfMuons, hNBins_pt, 0.0, hMax_pt);
		muH_eta_    = new TH1D(hNamePrefix + "eta",    "Muon #eta distribution (" + str_muonType_muons_cutsPhrase + "); Pseudorapidity, #eta;" + str_NoOfMuons, 40, -4.0, 4.0);
		muH_phi_    = new TH1D(hNamePrefix + "phi",    "Muon #phi distribution (" + str_muonType_muons_cutsPhrase + "); Azimuthal angle, #phi /rad;" + str_NoOfMuons, 40, 0, 3.1415);

		muH_isGlobMuon_       = new TH1D(hNamePrefix + "isGlobMuon",       "isGlobalMuon distribution (" + str_muonType_muons_cutsPhrase                     + "); isGlobalMuon;"                   + str_NoOfMuons, 2, -0.5, +1.5);
		muH_isTrkrMuon_       = new TH1D(hNamePrefix + "isTrkrMuon",       "isTrackerMuon distribution (" + str_muonType_muons_cutsPhrase                    + "); isTrackerMuon;"                  + str_NoOfMuons, 2, -0.5, +1.5);
		muH_isStandAloneMuon_ = new TH1D(hNamePrefix + "isStandAloneMuon", "isStandAloneMuon distribution (" + str_muonType_muons_cutsPhrase                 + "); isStandAloneMuon;"               + str_NoOfMuons, 2, -0.5, +1.5);
		muH_nMatchedMuStns_   = new TH1D(hNamePrefix + "nMatchedMuStns",   "N_{matched #mu stns} distribution (" + str_muonType_muons_cutsPhrase             + "); N_{matched #mu stns};"           + str_NoOfMuons, 36, -0.5, +35.5);
		muH_isolR03_sumPt_    = new TH1D(hNamePrefix + "isolR03_sumPt",    "#Sigma p_{T, trk} in #Delta R<0.3 distribution(" + str_muonType_muons_cutsPhrase + "); #Sigma p_{T} in #Delta R < 0.3;" + str_NoOfMuons, 25, 0.0, 100.0);

		muH_glob_nMuHits_    = new TH1D(hNamePrefix + "glob_nMuHits",    "N_{#mu hits} distribution (" + str_muonType_muons_cutsPhrase +     "); Number of muon hits, N_{#mu hits};"               + str_NoOfMuons, 26, -0.5, 25.5);
		muH_glob_normdChi2_  = new TH1D(hNamePrefix + "glob_normdChi2",  "Norm'd #chi^{2} distribution (" + str_muonType_muons_cutsPhrase +  "); Normalised #chi^2;"                               + str_NoOfMuons, 25, 0.0, 25.0);
		muH_inner_nPixHits_  = new TH1D(hNamePrefix + "inner_nPixHits",  "N_{Pixel hits} distribution (" + str_muonType_muons_cutsPhrase +   "); Number of valid pixel hits, N_{Pixel hits};"      + str_NoOfMuons, 26, -0.5, 25.5);
		muH_inner_nTrkrHits_ = new TH1D(hNamePrefix + "inner_nTrkrHits", "N_{Tracker hits} distribution (" + str_muonType_muons_cutsPhrase + "); Number of valid tracker hits, N_{Tracker hits};"  + str_NoOfMuons, 26, -0.5, 25.5);
		muH_inner_dxyOrigin_ = new TH1D(hNamePrefix + "inner_dxyOrigin", "d_{xy,origin} distribution (" + str_muonType_muons_cutsPhrase +    "); Impact parameter from origin, d_{xy,origin} /??;" + str_NoOfMuons, 25, -5.0, 5.0);

		// Set-up histograms so that errors are automatically calculated as the histos are filled
		std::vector<TH1D*> ptrsToHists = GetPtrsToHistos();
		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->Sumw2();
			ptrsToHists.at(iHist)->SetDirectory(0);
		}
	}

	//--------------------------------//
	//---- DTOR ...
	MuonDistns::~MuonDistns(){
		// Delete all of the histograms from the heap...
		delete muH_region_;
		delete muH_charge_;
		delete muH_p_;
		delete muH_pT_;
		delete muH_eta_;
		delete muH_phi_;
		delete muH_isGlobMuon_;
		delete muH_isTrkrMuon_;
		delete muH_isStandAloneMuon_;
		delete muH_nMatchedMuStns_;
		delete muH_isolR03_sumPt_;
		delete muH_glob_nMuHits_;
		delete muH_glob_normdChi2_;
		delete muH_inner_nPixHits_;
		delete muH_inner_nTrkrHits_;
		delete muH_inner_dxyOrigin_;
	}

	//--------------------------------//
	//--- 'Fill histograms' method ---//
	void MuonDistns::FillHistos(tsw::Muon theMuon, const Double_t muWeight){
		if( theMuon.isInBarrel() )
			muH_region_->Fill( 0.0, muWeight);
		else if( theMuon.isInOverlap() )
			muH_region_->Fill( 1.0, muWeight);
		else if( theMuon.isInEndcap() )
			muH_region_->Fill( 2.0, muWeight);
		else
			muH_region_->Fill( 3.0, muWeight);

		muH_charge_->Fill( theMuon.charge(), muWeight);
		muH_p_     ->Fill( theMuon.p()     , muWeight);
		muH_pT_    ->Fill( theMuon.pT()    , muWeight);
		muH_eta_   ->Fill( theMuon.eta()   , muWeight);
		muH_phi_   ->Fill( theMuon.phi()   , muWeight);

		muH_isGlobMuon_      ->Fill( theMuon.isGlobalMuon()      , muWeight);
		muH_isTrkrMuon_      ->Fill( theMuon.isTrackerMuon()     , muWeight);
		muH_isStandAloneMuon_->Fill( theMuon.isStandAloneMuon()  , muWeight);
		muH_nMatchedMuStns_  ->Fill( theMuon.numMatchedMuonStns(), muWeight);
		muH_isolR03_sumPt_   ->Fill( theMuon.isolR03_sumPt()     , muWeight);

		muH_glob_nMuHits_   ->Fill( theMuon.glob_numValidMuonHits(), muWeight);
		muH_glob_normdChi2_ ->Fill( theMuon.glob_normalisedChi2()  , muWeight);
		muH_inner_nPixHits_ ->Fill( theMuon.inner_numValidPixHits(), muWeight);
		muH_inner_nTrkrHits_->Fill( theMuon.inner_numValidTrkrHits(),muWeight);
		muH_inner_dxyOrigin_->Fill( theMuon.inner_dxyVsOrigin()    , muWeight);
	}

	void MuonDistns::SetHistAttributes(Int_t lineCol, Int_t lineStyle){
		std::vector<TH1D*> ptrsToHists;
		ptrsToHists = GetPtrsToHistos();

		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->SetLineColor(lineCol);
			ptrsToHists.at(iHist)->SetLineStyle(lineStyle);
		}
	}

	//-------------------------------------//
	//--- 'Write histos to file' method ---//
	void MuonDistns::WriteHistos(TFile* ptr_outFile){
		//Create a directory within the output file to store these histos in ...
		ptr_outFile->mkdir(hNamePrefix_);
		ptr_outFile->cd(hNamePrefix_);

		muH_region_->Write();
		muH_charge_->Write();
		muH_p_->Write();
		muH_pT_->Write();
		muH_eta_->Write();
		muH_phi_->Write();
		muH_isGlobMuon_->Write();
		muH_isTrkrMuon_->Write();
		muH_isStandAloneMuon_->Write();
		muH_nMatchedMuStns_->Write();
		muH_isolR03_sumPt_->Write();
		muH_glob_nMuHits_->Write();
		muH_glob_normdChi2_->Write();
		muH_inner_nPixHits_->Write();
		muH_inner_nTrkrHits_->Write();
		muH_inner_dxyOrigin_->Write();

		//Move back to the original directory within the output file ...
		ptr_outFile->cd("..");
	}

	//---------------------------------------------------------------//
	//--- Method for getting vector of pointers to all histograms ---//
	std::vector<TH1D*> MuonDistns::GetPtrsToHistos(){
		std::vector<TH1D*> tmpPtrsVector; tmpPtrsVector.clear();

		tmpPtrsVector.push_back( muH_region_);
		tmpPtrsVector.push_back( muH_charge_);
		tmpPtrsVector.push_back( muH_p_);
		tmpPtrsVector.push_back( muH_pT_);
		tmpPtrsVector.push_back( muH_eta_);
		tmpPtrsVector.push_back( muH_phi_);
		tmpPtrsVector.push_back( muH_isGlobMuon_);
		tmpPtrsVector.push_back( muH_isTrkrMuon_);
		tmpPtrsVector.push_back( muH_isStandAloneMuon_);
		tmpPtrsVector.push_back( muH_nMatchedMuStns_);
		tmpPtrsVector.push_back( muH_isolR03_sumPt_);
		tmpPtrsVector.push_back( muH_glob_nMuHits_);
		tmpPtrsVector.push_back( muH_glob_normdChi2_);
		tmpPtrsVector.push_back( muH_inner_nPixHits_);
		tmpPtrsVector.push_back( muH_inner_nTrkrHits_);
		tmpPtrsVector.push_back( muH_inner_dxyOrigin_);

		return tmpPtrsVector;
	}

} // End of tsw namespace

#endif
