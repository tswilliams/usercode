#ifndef tswDiEleDistnsByRgn_h
#define tswDiEleDistnsByRgn_h

#include "tswDiEleDistns.h"

namespace tsw{
	class DiEleDistnsByRgn{
		public:
			DiEleDistnsByRgn(TString hNamePrefix, const TString titleStr_diEleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt);
			~DiEleDistnsByRgn(){}

			// Fill method for histograms ...
			void FillHistos(tsw::HEEPDiEle diElectron, const Double_t diEleWeight);

			//"Write" method for the histos ...
			void WriteHistos(TFile* );

		private:
			//Member variables ...
			DiEleDistns diEleHistos_AllRgns_;
			DiEleDistns diEleHistos_EBEB_;
			DiEleDistns diEleHistos_EBEE_;
			DiEleDistns diEleHistos_noEEEE_;
			DiEleDistns diEleHistos_EEEE_;
	};

	DiEleDistnsByRgn::DiEleDistnsByRgn(TString hNamePrefix, const TString titleStr_diEleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_mass, Int_t hNBins_pt, Double_t hMax_pt):
		diEleHistos_AllRgns_(hNamePrefix+"All_", titleStr_diEleType, titleStr_cutsPhrase, lineColorIdx, lineStyleIdx, hNBins_mass, hNBins_pt, hMax_pt),
		diEleHistos_EBEB_(   hNamePrefix+"EBEB_", titleStr_diEleType+" EB-EB", titleStr_cutsPhrase, lineColorIdx, lineStyleIdx, hNBins_mass, hNBins_pt, hMax_pt),
		diEleHistos_EBEE_(   hNamePrefix+"EBEE_", titleStr_diEleType+" EB-EE", titleStr_cutsPhrase, lineColorIdx, lineStyleIdx, hNBins_mass, hNBins_pt, hMax_pt),
		diEleHistos_noEEEE_( hNamePrefix+"noEEEE_", titleStr_diEleType+" excl. EE-EE", titleStr_cutsPhrase, lineColorIdx, lineStyleIdx, hNBins_mass, hNBins_pt, hMax_pt),
		diEleHistos_EEEE_(   hNamePrefix+"EEEE_", titleStr_diEleType+" EE-EE", titleStr_cutsPhrase, lineColorIdx, lineStyleIdx, hNBins_mass, hNBins_pt, hMax_pt)
	{}

	void DiEleDistnsByRgn::FillHistos(tsw::HEEPDiEle diElectron, const Double_t diEleWeight){
		// All regions histogram is filled for all di-ele's ...
		diEleHistos_AllRgns_.FillHistos(diElectron, diEleWeight);

		if( diElectron.isEBEB() )
			diEleHistos_EBEB_.FillHistos(diElectron, diEleWeight);
		else if( diElectron.isEBEE() )
			diEleHistos_EBEE_.FillHistos(diElectron, diEleWeight);
		else if( diElectron.isEEEE() )
			diEleHistos_EEEE_.FillHistos(diElectron, diEleWeight);
		if( !(diElectron.isEEEE()) )
			diEleHistos_noEEEE_.FillHistos(diElectron, diEleWeight);
	}

	void DiEleDistnsByRgn::WriteHistos(TFile* p_outFile){
		diEleHistos_AllRgns_.WriteHistos(p_outFile);
		diEleHistos_EBEB_.WriteHistos(p_outFile);
		diEleHistos_EBEE_.WriteHistos(p_outFile);
		diEleHistos_noEEEE_.WriteHistos(p_outFile);
		diEleHistos_EEEE_.WriteHistos(p_outFile);
	}
}

#endif
