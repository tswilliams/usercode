#ifndef tswReconValidationHistos_h
#define tswReconValidationHistos_h

#include "tswHEEPVarDistns.h"

namespace tsw{

std::vector<unsigned int> IndicesOfElesPassingCuts(std::vector<tsw::HEEPEle> theElectrons, const std::vector<bool>& passCutsFlags, const int ecalRegionFlag);

	class ReconValidationHistos{
		public:
			ReconValidationHistos(const TString hNamePrefix, const TString titleStr_eleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, bool flag_writeHighestEtEleHists);
			~ReconValidationHistos(){}

			//Fill methods for the histograms ...
			void FillNumberElesHist(unsigned int numEles, const Double_t eventWeight);
			void FillHistos(tsw::HEEPEle, Double_t weight);
			void FillHistos(std::vector<tsw::HEEPEle> etOrderedEles, std::vector<bool> passCutsFlags, Double_t weight);

			//"Write" method for the histos ...
			void WriteHistos(TFile* outFile){
				eleHistos_EB_.WriteHistos(outFile);
				eleHistos_EE_.WriteHistos(outFile);
				//eleHistos_EBplusEE_.WriteHistos();
				if(writeHighestEtEleHists_){
					eleHistos_highestEtEle_EB_.WriteHistos(outFile);
					eleHistos_highestEtEle_EE_.WriteHistos(outFile);
					eleHistos_2ndHighestEtEle_EB_.WriteHistos(outFile);
					eleHistos_2ndHighestEtEle_EE_.WriteHistos(outFile);
					eleHistos_EBplusEE_.GetPointerToNumberElesHist()->Write();
				}
			}

		private:
			//Member variables ...
			HEEPVarDistns eleHistos_EB_;
			HEEPVarDistns eleHistos_EE_;
			HEEPVarDistns eleHistos_EBplusEE_;
			HEEPVarDistns eleHistos_highestEtEle_EB_;
			HEEPVarDistns eleHistos_highestEtEle_EE_;
			HEEPVarDistns eleHistos_2ndHighestEtEle_EB_;
			HEEPVarDistns eleHistos_2ndHighestEtEle_EE_;

			//Flag for whether to write out histograms for highest and 2nd highest Et eles ...
			bool writeHighestEtEleHists_;

	};

	ReconValidationHistos::ReconValidationHistos(const TString hNamePrefix, const TString titleStr_eleType, const TString titleStr_cutsPhrase, Int_t lineColorIdx, bool flag_writeHighestEtEleHists):
			eleHistos_EB_(                hNamePrefix+"EB_",      titleStr_eleType+" EB",                   titleStr_cutsPhrase),
			eleHistos_EE_(                hNamePrefix+"EE_",      titleStr_eleType+" EE",                   titleStr_cutsPhrase),
			eleHistos_EBplusEE_(          hNamePrefix+"both_",    titleStr_eleType+" (EB+EE)",              titleStr_cutsPhrase),
			eleHistos_highestEtEle_EB_(   hNamePrefix+"Et1_EB_",  titleStr_eleType+" highest EB E_{T}",     titleStr_cutsPhrase),
			eleHistos_highestEtEle_EE_(   hNamePrefix+"Et1_EE_",  titleStr_eleType+" highest EE E_{T}",     titleStr_cutsPhrase),
			eleHistos_2ndHighestEtEle_EB_(hNamePrefix+"Et2_EB_",  titleStr_eleType+" 2nd highest EB E_{T}", titleStr_cutsPhrase),
			eleHistos_2ndHighestEtEle_EE_(hNamePrefix+"Et2_EE_",  titleStr_eleType+" 2nd highest EE E_{T}", titleStr_cutsPhrase),
			writeHighestEtEleHists_(flag_writeHighestEtEleHists)
	{
		// Set histogram styles ...
		Int_t barrelLineStyle = 1; //1=>Solid line
		Int_t endcapLineStyle = 7; //7=>Dashed line
		eleHistos_EB_.SetHistAttributes(lineColorIdx, barrelLineStyle);
		eleHistos_EE_.SetHistAttributes(lineColorIdx, endcapLineStyle);
		eleHistos_EBplusEE_.SetHistAttributes(lineColorIdx, 1);
		eleHistos_highestEtEle_EB_.SetHistAttributes(lineColorIdx, barrelLineStyle);
		eleHistos_highestEtEle_EE_.SetHistAttributes(lineColorIdx, endcapLineStyle);
		eleHistos_2ndHighestEtEle_EB_.SetHistAttributes(lineColorIdx, barrelLineStyle);
		eleHistos_2ndHighestEtEle_EE_.SetHistAttributes(lineColorIdx, endcapLineStyle);
	}

	void ReconValidationHistos::FillNumberElesHist(unsigned int numEles, const Double_t eventWeight){
		eleHistos_EBplusEE_.FillNumberElesHist(numEles, eventWeight);
	}

	void ReconValidationHistos::FillHistos(tsw::HEEPEle theEle, Double_t weight){
		eleHistos_EBplusEE_.FillHistos( theEle, weight);
		if(theEle.isEB()){
			eleHistos_EB_.FillHistos( theEle, weight);}
		if(theEle.isEE())
			eleHistos_EE_.FillHistos( theEle, weight);
	}

	void ReconValidationHistos::FillHistos(std::vector<tsw::HEEPEle> etOrderedEles, std::vector<bool> passCutsFlags, Double_t weight){
		//General variable declarations
		tsw::HEEPEle ithEle;
		// Declare variables that acts as a counters for the number of electrons that have passed cuts ...
		unsigned int numElesPassed_EB = 0;
		unsigned int numElesPassed_EE = 0;
		// Declare variables that store the indices
		std::vector<unsigned int> idxs_elesPassedCuts   = IndicesOfElesPassingCuts(etOrderedEles, passCutsFlags, 0);
		//std::vector<unsigned int> idxs_EBelesPassedCuts = IndicesOfElesPassingCuts(etOrderedEles, passCutsFlags, 1);
		//std::vector<unsigned int> idxs_EEelesPassedCuts = IndicesOfElesPassingCuts(etOrderedEles, passCutsFlags, 2);

		// Run over all electrons ...
		for(unsigned int iEle=0; iEle<etOrderedEles.size(); iEle++){
			ithEle = etOrderedEles.at(iEle);
			//Fill the 'all-ele' histograms if electron passed the cuts ...
			if( passCutsFlags.at(iEle) ){
				eleHistos_EBplusEE_.FillHistos( ithEle, weight);
				if(ithEle.isEB()){
					eleHistos_EB_.FillHistos( ithEle, weight);
					numElesPassed_EB++;}
				if(ithEle.isEE()){
					eleHistos_EE_.FillHistos( ithEle, weight);
					numElesPassed_EE++;}
			}
		}

		//Find the indices of the highest two Et eles ...
		int idx_highestEtEle = -1;
		tsw::HEEPEle highestEtEle;
		int idx_2ndHighestEtEle = -1;
		tsw::HEEPEle secondHighestEtEle;

		// There can only be a highest E_T electron if at least one electron has passed the cuts ...
		if(idxs_elesPassedCuts.size()>0){
			idx_highestEtEle = idxs_elesPassedCuts.at(0);
			highestEtEle = etOrderedEles.at(idx_highestEtEle);
			// Fill the highest E_T electron histograms - with EE and EB plots corresponding to the whole-detector highest E_T electrons that happened to be in that region ...
			if( highestEtEle.isEB() )
				eleHistos_highestEtEle_EB_.FillHistos( highestEtEle, weight);
			else if( highestEtEle.isEE() )
				eleHistos_highestEtEle_EE_.FillHistos( highestEtEle, weight);

			// There can only be a 2nd highest E_T electron if at least two electrons have passed the cuts ...
			if(idxs_elesPassedCuts.size()>1){
				idx_2ndHighestEtEle = idxs_elesPassedCuts.at(1);
				secondHighestEtEle = etOrderedEles.at( idx_2ndHighestEtEle );
				// Fill the 2nd highest E_T electron histograms - with EE and EB plots corresponding to the whole-detector 2nd highest E_T electrons that happened to be in that region ...
				if( secondHighestEtEle.isEB() )
					eleHistos_2ndHighestEtEle_EB_.FillHistos( secondHighestEtEle, weight);
				else if( secondHighestEtEle.isEE() )
					eleHistos_2ndHighestEtEle_EE_.FillHistos( secondHighestEtEle, weight);
			}
		}

		if( idxs_elesPassedCuts.size()!=(numElesPassed_EB + numElesPassed_EE) ){
			std::cout << "  **********************************************************************" << std::endl;
			std::cout << "  *** WARNING - Error in FillHistos method of ReconValidationHistos ..." << std::endl;
			std::cout << "  **********************************************************************" << std::endl;
		}

		eleHistos_EBplusEE_.GetPointerToNumberElesHist()->Fill(idxs_elesPassedCuts.size());

		//Fill the number of electrons histogram for each region ...
		eleHistos_EB_.FillNumberElesHist(numElesPassed_EB, 1.0);
		eleHistos_EE_.FillNumberElesHist(numElesPassed_EE, 1.0);
	}

	std::vector<unsigned int> IndicesOfElesPassingCuts(std::vector<tsw::HEEPEle> theElectrons, const std::vector<bool>& passCutsFlags, const int ecalRegionFlag){
		//ecalRegionFlag:
		//			=0 => Consider electrons in all regions
		//			=1 => Consider only EB electrons
		//			=2 => Consider only EE electrons
		std::vector<unsigned int> indicesVector;
		indicesVector.clear();

		if(ecalRegionFlag==0){
			for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
				if(passCutsFlags.at(idx))
					indicesVector.push_back(idx);
			}
		}
		else if(ecalRegionFlag==1){
			for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
				if( passCutsFlags.at(idx) && theElectrons.at(idx).isEB() )
					indicesVector.push_back(idx);
			}
		}
		else if(ecalRegionFlag==2){
			for(unsigned int idx=0; idx<passCutsFlags.size(); idx++){
				if( passCutsFlags.at(idx) && theElectrons.at(idx).isEE() )
					indicesVector.push_back(idx);
			}
		}

		return indicesVector;
	}

	unsigned int NumPassingCuts(std::vector<bool> passCutsFlags){
		unsigned int numberElesPassingCuts = 0;

		//Loop over the flags and count how many occurences of "true" there are ...
		for(unsigned int iEle=0; iEle<passCutsFlags.size(); iEle++){
			if(passCutsFlags.at(iEle))
				numberElesPassingCuts++;
		}

		return numberElesPassingCuts;
	}

}
#endif
