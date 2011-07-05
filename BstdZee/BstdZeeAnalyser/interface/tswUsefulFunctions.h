#ifndef tswUsefulFunctions_h
#define tswUsefulFunctions_h

TLorentzVector ConvertToTLorentzVector(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >*);
void SetBinomialErrorsOnNumerator(TH1D* numeratorHist, TH1D* denomHist);
void SetPoissonErrorsOnHist(TH1D* myHist, Double_t eventWeight);
void CalcPreWeightEffiHist(TH1D* effiHist, TH1D* numeratorHist, TH1D* denomHist, Double_t eventWeight);


TLorentzVector ConvertToTLorentzVector(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >* p4_old)
{
	TLorentzVector p4_new;
	p4_new.SetPxPyPzE(p4_old->Px(), p4_old->Py(), p4_old->Pz(), p4_old->E());
	return p4_new;
}

void SetBinomialErrorsOnNumerator(TH1D* numeratorHist, TH1D* denomHist){
	Int_t numBins = numeratorHist->GetNbinsX();
	Double_t binError = -999.9;
	if(numeratorHist->GetNbinsX()!=denomHist->GetNbinsX())
		std::cout << std::endl << "****ERROR in assigning errors in SetBinomialErrorsOnNumerator function***" << std::endl;

	for(int binIdx = 0; binIdx<=numBins+1; binIdx++){
		if(numeratorHist->GetBinContent(binIdx)==0.0 || denomHist->GetBinContent(binIdx)==0.0)
			binError=0.0;
		else
			binError = sqrt( numeratorHist->GetBinContent(binIdx)*(1.0 - numeratorHist->GetBinContent(binIdx)/denomHist->GetBinContent(binIdx)) );
		numeratorHist->SetBinError(binIdx, binError);
	}

	return;
}

void SetPoissonErrorsOnHist(TH1D* myHist, Double_t eventWeight){
	Int_t numBins = myHist->GetNbinsX();
	Double_t binError = -999.9;

	for(int binIdx = 0; binIdx<=numBins+1; binIdx++){
		binError = eventWeight*sqrt( myHist->GetBinContent(binIdx)/eventWeight );
		myHist->SetBinError(binIdx, binError);
	}

	return;
}

void CalcPreWeightEffiHist(TH1D* effiHist, TH1D* numeratorHist, TH1D* denomHist, Double_t eventWeight){
	Int_t numBins = effiHist->GetNbinsX();
	Double_t binValue = -999.9;
	Double_t binError = -999.9;

	if(numeratorHist->GetNbinsX()!=denomHist->GetNbinsX() || numeratorHist->GetNbinsX()!=effiHist->GetNbinsX())
		std::cout << std::endl << "****ERROR in assigning errors in SetBinomialErrorsOnNumerator function***" << std::endl;

	for(int binIdx = 0; binIdx<=numBins+1; binIdx++){
		if(numeratorHist->GetBinContent(binIdx)==0.0 || denomHist->GetBinContent(binIdx)==0.0)
			binError=0.0;
		else{
			binValue = numeratorHist->GetBinContent(binIdx)/denomHist->GetBinContent(binIdx);
			binError = sqrt( binValue*(1.0-binValue)/(denomHist->GetBinContent(binIdx)/eventWeight) );
			effiHist->SetBinContent(binIdx, binValue);
			effiHist->SetBinError(binIdx, binError);
		}
	}

	return;
}

#endif
