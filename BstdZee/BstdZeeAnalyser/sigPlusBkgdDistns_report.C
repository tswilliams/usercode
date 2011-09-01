{
	// ROOT script to plot stacked signal + background distributions for the boostd
	// Distributions plotted: Mass spectrum, pT spectrum ...
	// ... then: dEta, dPhi, dR and opening angle distributions
	// ... Also: highest Et ele from diele - et distn.
	// In each case for di-eles formed from modified reco'n electrons, passing all HEEP cuts apart from isol. ones.

	// Including files that are
	//#include "interface/compareDistns.h"
	#include <vector>

	#include "tdrstyle_tsw.C";
	setTDRStyle();

	std::cout << "Starting..." << std::endl;

	//Flag for whether to scale sum of background histograms to data ...
	bool scaleBkgdMCToData = true;
	Double_t mcScaleFactor = 1.0;

	// Defining which distributions to plot ...
	std::vector<TString> histosToPlot; histosToPlot.clear();
	std::vector<Double_t> histosYLims; histosYLims.clear();
	std::vector<TString> histosXTitles; histosXTitles.clear();

	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_sumEt");   histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Transverse energy of Z, E_{T,ee} /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_sumEt");   histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Transverse energy of Z, E_{T,ee} /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_sumEt");   histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Transverse energy of Z, E_{T,ee} /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_sumEt");   histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Transverse energy of Z, E_{T,ee} /GeV");
	//2nd canvas ...
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_invMass"); histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Z candidate mass, M_{ee} /GeVc^{-2}");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_invMass"); histosYLims.push_back( 7500.0 );  histosXTitles.push_back("Z candidate mass, M_{ee} /GeVc^{-2}");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_regions"); histosYLims.push_back( 300000.0 );  histosXTitles.push_back("Di-electron regions");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_regions"); histosYLims.push_back( 100000.0 );  histosXTitles.push_back("Di-electron regions");
	//3rd canvas ...
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_openAngle"); histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Di-electron opening angle, #theta_{ee}");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_deltaR");    histosYLims.push_back( 300000.0 ); histosXTitles.push_back("#DeltaR_{ee}");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_deltaEta");  histosYLims.push_back( 300000.0 ); histosXTitles.push_back("#Delta#eta_{ee}");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_deltaPhi");  histosYLims.push_back( 300000.0 ); histosXTitles.push_back("#Delta#phi_{ee}");
	//4th canvas ...
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleA_et");     histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Transverse energy, E_{T} /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleA_energy"); histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Energy, E /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleA_eta");    histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Pseudorapidity, #eta");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleA_phi");    histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Azimuthal angle, #phi");
	//5th canvas ...
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleB_et");     histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Transverse energy, E_{T} /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleB_energy"); histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Energy, E /GeV");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleB_eta");    histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Pseudorapidity, #eta");
	histosToPlot.push_back("h_normDiEle_EB_HEEPNoIso_MZ_trgA_/h_normDiEle_EB_HEEPNoIso_MZ_trgA_eleB_phi");    histosYLims.push_back( 300000.0 ); histosXTitles.push_back("Azimuthal angle, #phi");
	//6th canvas ...
	histosToPlot.push_back("h_eleMu_EB_HEEPNoIso_tight_MZ_/h_eleMu_EB_HEEPNoIso_tight_MZ_sumEt");    histosYLims.push_back( 10000.0 ); histosXTitles.push_back( "Transverse momentum, p_{T,e#mu} /GeVc^{-1}" );
	histosToPlot.push_back("h_eleMu_EB_HEEPNoIso_tight_MZ_/h_eleMu_EB_HEEPNoIso_tight_MZ_sumEt");    histosYLims.push_back( 10000.0 ); histosXTitles.push_back( "Transverse momentum, p_{T,e#mu} /GeVc^{-1}" );
	histosToPlot.push_back("h_eleMu_EB_HEEPNoIso_tight_MZ_/h_eleMu_EB_HEEPNoIso_tight_MZ_sumEt");    histosYLims.push_back( 10000.0 ); histosXTitles.push_back( "Transverse momentum, p_{T,e#mu} /GeVc^{-1}" );
	histosToPlot.push_back("h_eleMu_EB_HEEPNoIso_tight_MZ_/h_eleMu_EB_HEEPNoIso_tight_MZ_invMass");  histosYLims.push_back( 10000.0 ); histosXTitles.push_back( "Invariant mass, M_{e#mu} /GeVc^{-2}" );

	//Defining the files for the background and signal distributions ...
	std::vector<TString> bkgdFNames; bkgdFNames.clear();
	std::vector<TString> bkgdLabels; bkgdLabels.clear();
	std::vector<Int_t> bkgdFillColours; bkgdFillColours.clear();
	std::vector<Int_t> bkgdLineColours; bkgdLineColours.clear();

	bkgdFNames.push_back("2011-08-04/histos_ZZTo2L2Nu_2011-08-04.root");  bkgdLabels.push_back("Z(ll)Z(#nu#nu)");  bkgdFillColours.push_back(8); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_WZTo3LNu_2011-08-04.root");   bkgdLabels.push_back("Z(ll)W(l'#nu)");   bkgdFillColours.push_back(3); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_WWTo2L2Nu_2011-08-04.root");  bkgdLabels.push_back("W(l#nu)W(l'#nu)"); bkgdFillColours.push_back(6); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_DYJetsToLL-TauTau_2011-08-04.root"); bkgdLabels.push_back("DY(#tau#tau)+Jets");   bkgdFillColours.push_back(9); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_TTbar_2011-08-04.root");      bkgdLabels.push_back("ttbar");           bkgdFillColours.push_back(2); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_WJetsToLNu_2011-08-04.root"); bkgdLabels.push_back("W(l#nu) + Jets");  bkgdFillColours.push_back(4); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_QCD-dblEM_2011-08-04.root"); bkgdLabels.push_back("QCD");  bkgdFillColours.push_back(16); bkgdLineColours.push_back(1);
	bkgdFNames.push_back("2011-08-04/histos_DYJetsToLL-ee_2011-08-04.root"); bkgdLabels.push_back("DY(ee)+Jets");   bkgdFillColours.push_back(7); bkgdLineColours.push_back(1);

	std::vector<TString> sigFNames; sigFNames.clear();
	std::vector<TString> sigLabels; sigLabels.clear();
	std::vector<Int_t> sigLineColours; sigLineColours.clear();
	sigFNames.push_back("2011-08-04/histos_0-75TeVu_2011-08-01.root"); sigLabels.push_back("0.75TeV u*"); sigLineColours.push_back(3);
	sigFNames.push_back("2011-08-04/histos_1-00TeVu_2011-08-01.root"); sigLabels.push_back("1.00TeV u*"); sigLineColours.push_back(6);
	sigFNames.push_back("2011-08-04/histos_2-00TeVu_2011-08-01.root"); sigLabels.push_back("2.00TeV u*"); sigLineColours.push_back(2);

	TString dataFName = "2011-08-04/histos_data-Photon_0-5invfb_2011-08-04.root";
	TString dataMuEGFName = "2011-08-04/histos_data-MuEG_0-5invfb_2011-08-04.root";

	TString outputFileName   = "2011-08-06/sigPlusBkgdDistns_2011-08-06_EB_HEEPNoIso_MZ_trgA_";

	std::vector<TCanvas*> vecOfCanvas; vecOfCanvas.clear();
	std::vector<TString> canvasTag; canvasTag.clear();
	canvasTag.push_back("DiEle1");
	canvasTag.push_back("DiEle2");
	canvasTag.push_back("DiEle3");
	canvasTag.push_back("LeadingEle");
	canvasTag.push_back("SubLeadingEle");
	canvasTag.push_back("eleMu");

	//Do plot for each distribution in turn ...
	for(unsigned int plotIdx=0; plotIdx<histosToPlot.size(); plotIdx++){
		//Set flag for whether this graph is cumulative or not ...
		bool calcCumulative = (plotIdx==1) || (plotIdx==3 || (plotIdx==22) );
		//Setting up the stack of histograms ...
		std::vector<THStack*> mcStacks; mcStacks.clear();
		THStack* ptStack  = new THStack("ptStack","Signal plus background pT distribution; p_{T ee} /GeVc^{-1}; Number of events per XX GeVc^{-1}");

		TLegend* myLegend = new TLegend(0.7,0.32,0.94,0.85);
		myLegend->SetLineColor(0);
		myLegend->SetFillColor(0);
		myLegend->SetShadowColor(0);
		myLegend->SetBorderSize(0);

		Double_t xAxisMax = 0.000;
		Double_t xAxisMin = 0.000;
		//Grabbing the background distributions and adding to stack ...
		Double_t tmp_newBinContent = -999.9;
		for(unsigned int bkgdFileIdx = 0; bkgdFileIdx<bkgdFNames.size(); bkgdFileIdx++){
			TFile* inFile = TFile::Open( bkgdFNames.at(bkgdFileIdx) );
			TH1D* histo = (TH1D*)inFile->Get(histosToPlot.at(plotIdx));

			histo->SetFillColor( bkgdFillColours.at(bkgdFileIdx) );
			histo->SetLineColor( bkgdLineColours.at(bkgdFileIdx) );
			histo->SetDirectory(0);
			xAxisMin = histo->GetXaxis()->GetXmin();
			xAxisMax = histo->GetXaxis()->GetXmax();

			if(calcCumulative){
				//Compute the integral histogram here ...
				for(int binIdx=histo->GetNbinsX(); binIdx>=0; binIdx--){
					//std::cout << "  binIdx = " << binIdx << std::endl;
					histo->AddBinContent(binIdx, histo->GetBinContent(binIdx+1) );
				}
			}

			if( (histosToPlot.at(plotIdx).EndsWith("Et") || (histosToPlot.at(plotIdx).EndsWith("EtLog")) && !calcCumulative) ){
				for(int binIdx=1; binIdx<=histo->GetNbinsX(); binIdx++){
					histo->SetBinContent(binIdx, histo->GetBinContent(binIdx)/(histo->GetBinWidth(binIdx)/20.0) );  // was 10.0
				}
			}

			//Scale histogram ...
			if(scaleBkgdMCToData && ( (plotIdx>=2 && plotIdx<20) || (plotIdx>20) ))
				histo->Scale(mcScaleFactor);

			ptStack->Add( histo );
			TString legendSuffix = "";
//			if(calcCumulative)
//				legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), histo->GetBinContent(1));
//			else
//				legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), histo->Integral()+histo->GetBinContent(histo->GetNbinsX()+1));
			myLegend->AddEntry(histo, bkgdLabels.at(bkgdFileIdx) + legendSuffix,"f");

			inFile->Close();
		}

		if(plotIdx%4==0){
			TCanvas* newCanvas = new TCanvas("sigPlusBkgd_" + canvasTag.at(plotIdx/4), "Di-ele distributions " + canvasTag.at(plotIdx/4),1200,900);
			vecOfCanvas.push_back( newCanvas );
		}
		TCanvas* ptCanvas = vecOfCanvas.at(plotIdx/4);
		if(plotIdx%4==0)
			ptCanvas->Divide(2,2);
		TPad* currentPad = (TPad*)ptCanvas->cd(1+plotIdx%4);
		if(plotIdx!=7 && plotIdx!=5)
			currentPad->SetLogy();
		currentPad->DrawFrame(xAxisMin, 0.1, xAxisMax, histosYLims.at(plotIdx), "; " + histosXTitles.at(plotIdx) + "; ");

		// Grabbing the data distribution ...
		if(!histosToPlot.at(plotIdx).Contains("eleMu"))
			TFile* inFile = TFile::Open( dataFName );
		else{
			TFile* inFile = TFile::Open( dataMuEGFName );
			std::cout << "Using MuEG input file ..." << std::endl;}
		TH1D* histo = (TH1D*)inFile->Get(histosToPlot.at(plotIdx));
		histo->SetLineColor(1);
		histo->SetMarkerStyle(7);
		histo->Sumw2();

		// Computing the ratio of the number of data points to the number of MC points ...
		if(plotIdx==0 || plotIdx==20){
			TList* listOfBkgdMCHists = ptStack->GetHists();
			dataHist = histo;

			Double_t predNumBkgdMCEvts = 0.0;
			for(unsigned int bkgdHistIdx=0; bkgdHistIdx<listOfBkgdMCHists->GetSize(); bkgdHistIdx++){
				TH1D* ithBkgdHist = (TH1D*)listOfBkgdMCHists->At(bkgdHistIdx);
				predNumBkgdMCEvts += ithBkgdHist->Integral() + ithBkgdHist->GetBinContent(ithBkgdHist->GetNbinsX()+1);
				std::cout << "   " << bkgdFNames.at(bkgdHistIdx) << ":  " << ithBkgdHist->Integral() + ithBkgdHist->GetBinContent(ithBkgdHist->GetNbinsX()+1) << " events" << std::endl;
			}

			Double_t numOfDataEvts = dataHist->Integral() + dataHist->GetBinContent(dataHist->GetNbinsX()+1);
			mcScaleFactor = numOfDataEvts/predNumBkgdMCEvts;
			std::cout << "** N.B: Num data evts / predicted num bkgd MC evts = " << mcScaleFactor << " (=" << numOfDataEvts << "/" << predNumBkgdMCEvts << ")" << std::endl;
			if(scaleBkgdMCToData)
				std::cout << " * MC will now be scaled by this factor" << std::endl;
			else
				std::cout << " * But, bkgd MC will NOT be scaled." << std::endl;
		}

		if(!histosToPlot.at(plotIdx).Contains("eleMu")){
			///Grabbing the signal distribution and adding to stack ...
			for(unsigned int sigFileIdx = 0; sigFileIdx<sigFNames.size(); sigFileIdx++){
				// Clone ptStack, and place pointer to clone in mcStacks ...
				THStack* sigPlusBkgdStack = (THStack*)ptStack->Clone();
				mcStacks.push_back( sigPlusBkgdStack );

				TFile* sig_inFile = TFile::Open( sigFNames.at(sigFileIdx) );
				TH1D* sigHisto = (TH1D*)sig_inFile->Get(histosToPlot.at(plotIdx));

				sigHisto->SetLineColor( sigLineColours.at(sigFileIdx) );
				sigHisto->SetDirectory(0);

				if(calcCumulative){
					//Compute the integral sigHistogram here ...
					for(int binIdx=sigHisto->GetNbinsX(); binIdx>=0; binIdx--){
						sigHisto->AddBinContent(binIdx, sigHisto->GetBinContent(binIdx+1) );
					}
				}

				if( histosToPlot.at(plotIdx).EndsWith("Et") || (histosToPlot.at(plotIdx).EndsWith("EtLog") && !calcCumulative) ){
					for(int binIdx=1; binIdx<=sigHisto->GetNbinsX(); binIdx++){
						sigHisto->SetBinContent(binIdx, sigHisto->GetBinContent(binIdx)/(sigHisto->GetBinWidth(binIdx)/20.0) );  // was 10.0
					}
				}
				//Scale histogram ...
				if(scaleBkgdMCToData && plotIdx>=2)
					sigHisto->Scale(mcScaleFactor);

				sigPlusBkgdStack->Add( sigHisto );
				TString legendSuffix = "";
//				if(calcCumulative)
//					legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), sigHisto->GetBinContent(1));
//				else
//					legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), sigHisto->Integral()+sigHisto->GetBinContent(sigHisto->GetNbinsX()+1));
				myLegend->AddEntry(sigHisto, sigLabels.at(sigFileIdx) + legendSuffix,"l");

				sig_inFile->Close();

				// Plot the signal + background distribution ...
				sigPlusBkgdStack->Draw("SAME");
				sigPlusBkgdStack->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
				sigPlusBkgdStack->GetHistogram()->GetXaxis()->CenterTitle();
				sigPlusBkgdStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
				sigPlusBkgdStack->GetHistogram()->GetYaxis()->CenterTitle();
				currentPad->Update(); currentPad->Modified();
				gPad->RedrawAxis();
			}
		}
		else{
			ptStack->Draw("SAME");
			ptStack->GetHistogram()->GetXaxis()->SetTitleOffset(0.75);
			ptStack->GetHistogram()->GetXaxis()->CenterTitle();
			ptStack->GetHistogram()->GetYaxis()->SetTitleOffset(0.75);
			ptStack->GetHistogram()->GetYaxis()->CenterTitle();
			currentPad->Update(); currentPad->Modified();
			gPad->RedrawAxis();
		}




		if(calcCumulative){
			//Compute the integral histogram here ...
			for(int binIdx=histo->GetNbinsX(); binIdx>=0; binIdx--){
				histo->AddBinContent(binIdx, histo->GetBinContent(binIdx+1) );
			}
		}
		//Compute the statistical error for each bin ...
		for(int binIdx=0; binIdx<=histo->GetNbinsX()+1; binIdx++)
			histo->SetBinError(binIdx, sqrt(histo->GetBinContent(binIdx)) );

		//Compute scale the
		if( histosToPlot.at(plotIdx).EndsWith("Et") || (histosToPlot.at(plotIdx).EndsWith("EtLog") && !calcCumulative)) ){
			for(int binIdx=1; binIdx<=histo->GetNbinsX(); binIdx++){
				Double_t binScaleFactor = 1.0/(histo->GetBinWidth(binIdx)/20.0); // was 10.0
				Double_t binError = histo->GetBinError(binIdx)*binScaleFactor;
				histo->SetBinContent(binIdx, histo->GetBinContent(binIdx)*binScaleFactor );
				histo->SetBinError(binIdx, binError);
			}
		}

		histo->Draw("SAME");
		histo->SetDirectory(0);
		TString legendSuffix = "";
//		if(calcCumulative)
//			legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), histo->GetBinContent(1));
//		else
//			legendSuffix = TString::Format( (const char*)(legendSuffix+" [%.2e]"), histo->Integral());
		myLegend->AddEntry(histo, "Data "+legendSuffix, "le");
		inFile->Close();

		if(plotIdx%4==1 || plotIdx==2)
			myLegend->Draw("");

		if(plotIdx%4==1){
			TPaveText* textBox = new TPaveText(0.22, 0.78, 0.45, 0.84, "NDC");
			textBox->AddText("Run2011A, 0.5fb^{-1}");
			textBox->Draw("");
		}

		if( plotIdx%4==3 || plotIdx==(histosToPlot.size()-1) ){
			ptCanvas->SaveAs(outputFileName + canvasTag.at(plotIdx/4) + ".pdf");
			ptCanvas->SaveAs(outputFileName + canvasTag.at(plotIdx/4) + ".png");}
	}


	TFile* outFile = TFile::Open(outputFileName,"RECREATE");
	outFile->cd("");
	ptCanvas->Write();
	outFile->Close();

	/*delete ptStack;
	delete myLegend;
	delete ptCanvas;
	delete textBox;*/
}
