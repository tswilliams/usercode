#ifndef tswMCParticleDistns_h
#define tswMCParticleDistns_h

#include "BstdZeeFirst/Analyser/interface/tswMCParticle.h"

namespace tsw{
	class MCParticleDistns{
		public:
			//Constructors and destructors ...
			MCParticleDistns(TString hNamePrefix, const TString titleStr_particle, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_pt, Double_t hMax_pt);
			~MCParticleDistns();

			//Method to fill the histograms
			void FillHistos(tsw::MCParticle theMCParticle, const Double_t particleWeight);
			void FillHistos(tsw::MCParticle theMCParticle, Double_t daughter_dR, Double_t daughter_dEta, Double_t daughter_dPhi, Double_t daughter_openingAngle, const Double_t particleWeight);

			//Method for setting the line attributes for the histograms ...
			void SetHistAttributes(Int_t lineCol, Int_t lineStyle);

			//Method to write the histograms to file ...
			void WriteHistos(TFile* ptr_outFile);

			//Method to get vector of pointers to all histograms ...
			std::vector<TH1D*> GetPtrsToHistos();

		private:
			TString hNamePrefix_;
			TH1D* hist_pdgId_;
			TH1D* hist_hepMCstatus_;
			TH1D* hist_charge_;
			TH1D* hist_p_;
			TH1D* hist_pT_;
			TH1D* hist_eta_;
			TH1D* hist_phi_;
			TH1D* hist_daughter_dR_;
			TH1D* hist_daughter_dEta_;
			TH1D* hist_daughter_dPhi_;
			TH1D* hist_daughter_openingAngle_;
	};

	/////////////////////
	// Constructor
	MCParticleDistns::MCParticleDistns(TString hNamePrefix, const TString titleStr_particle, const TString titleStr_cutsPhrase, Int_t lineColorIdx, Int_t lineStyleIdx, Int_t hNBins_pt, Double_t hMax_pt):
		hNamePrefix_(hNamePrefix)
	{
		//Form standard phrases ...
		TString str_particleType_cutsPhrase = titleStr_particle + ", " + titleStr_cutsPhrase;
		TString str_NoOfParticles = "Number";

		//Initialise MC particle histograms ...
		hist_pdgId_       = new TH1D(hNamePrefix + "pdgId",   "PDG ID distribution  (" + str_particleType_cutsPhrase + "); PDG ID;" + str_NoOfParticles + "", 41, -20.5, 20.5);
		hist_hepMCstatus_ = new TH1D(hNamePrefix + "hepMCstatus",     "HEP MC status distribution  (" + str_particleType_cutsPhrase + "); HEP MC status;" + str_NoOfParticles, 11, -0.5, 10.5);
		hist_charge_      = new TH1D(hNamePrefix + "charge", "MC charge distribution ("   + str_particleType_cutsPhrase + "); MC particle charge /e;"  + str_NoOfParticles,                   5, -2.5, 2.5);
		hist_p_           = new TH1D(hNamePrefix + "p",      "MC momentum distribution (" + str_particleType_cutsPhrase + "); MC momentum /GeVc^{-1};" + str_NoOfParticles,                   hNBins_pt, 0.0, hMax_pt);
		hist_pT_          = new TH1D(hNamePrefix + "pT",     "MC p_{T} distribution ("    + str_particleType_cutsPhrase + "); MC transverse momentum, p_{T} /GeVc^{-1};" + str_NoOfParticles, hNBins_pt, 0.0, hMax_pt);
		hist_eta_         = new TH1D(hNamePrefix + "eta",    "MC #eta distribution ("     + str_particleType_cutsPhrase + "); MC pseudorapidity, #eta;"                  + str_NoOfParticles, 40, -4.0, +4.0);
		hist_phi_         = new TH1D(hNamePrefix + "phi",    "MC #phi distribution ("     + str_particleType_cutsPhrase + "); MC azimuthal angle, #phi;"                 + str_NoOfParticles, 40, -3.1416, +3.1416);

		hist_daughter_dR_           = new TH1D(hNamePrefix + "daughter_dR",        "MC daughter #Delta R ("      + str_particleType_cutsPhrase + "); MC "+titleStr_particle+" daughter #Delta R;"              + str_NoOfParticles, 50, 0.0, 5.0);
		hist_daughter_dEta_         = new TH1D(hNamePrefix + "daughter_dEta",      "MC daughter #Delta#eta ("    + str_particleType_cutsPhrase + "); MC "+titleStr_particle+" daughter #Delta#eta;"            + str_NoOfParticles, 50, -5.0, +5.0);
		hist_daughter_dPhi_         = new TH1D(hNamePrefix + "daughter_dPhi",      "MC daughter #Delta#phi ("    + str_particleType_cutsPhrase + "); MC "+titleStr_particle+" daughter #Delta#phi;"            + str_NoOfParticles, 50, -3.1416, +3.1416);
		hist_daughter_openingAngle_ = new TH1D(hNamePrefix + "daughter_openAngle", "MC daughter opening angle (" + str_particleType_cutsPhrase + "); MC "+titleStr_particle+" daughter opening angle, #theta;" + str_NoOfParticles, 50, 0.0, 3.1416);

		SetHistAttributes(lineColorIdx, lineStyleIdx);
	}

	/////////////////////
	// Destructor
	MCParticleDistns::~MCParticleDistns(){
		// Delete all of the histograms from the heap ...
		delete hist_pdgId_;
		delete hist_hepMCstatus_;
		delete hist_charge_;
		delete hist_p_;
		delete hist_pT_;
		delete hist_eta_;
		delete hist_phi_;
		delete hist_daughter_dR_;
		delete hist_daughter_dEta_;
		delete hist_daughter_dPhi_;
		delete hist_daughter_openingAngle_;
	}

	/////////////////////////////////
	// Methods to fill histograms ...
	void MCParticleDistns::FillHistos(tsw::MCParticle parentParticle, const Double_t particleWeight){
		//Fill each of the histograms ...
		hist_pdgId_->Fill(parentParticle.pdgId(), particleWeight);
		hist_hepMCstatus_->Fill(parentParticle.hepMCstatus(), particleWeight);
		hist_charge_->Fill(parentParticle.charge(), particleWeight);
		hist_p_->Fill(parentParticle.p4().P(), particleWeight);
		hist_pT_->Fill(parentParticle.pT(), particleWeight);
		hist_eta_->Fill(parentParticle.eta(), particleWeight);
		hist_phi_->Fill(parentParticle.phi(), particleWeight);
	}
	void MCParticleDistns::FillHistos(tsw::MCParticle parent, Double_t daughter_dR, Double_t daughter_dEta, Double_t daughter_dPhi, Double_t daughter_openingAngle, const Double_t particleWeight){
		//Fill each of the parent particle's histograms ...
		FillHistos(parent, particleWeight);

		//Fill each of the daughter particles' histograms ...
		hist_daughter_dR_->Fill(  daughter_dR, particleWeight);
		hist_daughter_dEta_->Fill(daughter_dEta, particleWeight);
		hist_daughter_dPhi_->Fill(daughter_dPhi, particleWeight);
		hist_daughter_openingAngle_->Fill(daughter_openingAngle, particleWeight);
	}

	void MCParticleDistns::SetHistAttributes(Int_t lineCol, Int_t lineStyle){
		std::vector<TH1D*> ptrsToHists;
		ptrsToHists = GetPtrsToHistos();

		for(unsigned int iHist=0; iHist<ptrsToHists.size() ; iHist++){
			ptrsToHists.at(iHist)->SetLineColor(lineCol);
			ptrsToHists.at(iHist)->SetLineStyle(lineStyle);
		}
	}

	/////////////////////////////////////////////////
	// Method to write out the histograms to file ...
	void MCParticleDistns::WriteHistos(TFile* ptr_outFile){
		//Create a directory within the output file to store these histos in ...
		TDirectory* tmpTDirPtr;
		tmpTDirPtr = ptr_outFile->mkdir(hNamePrefix_);
		ptr_outFile->cd(hNamePrefix_);

		//Write out the di-ele histos ...
		hist_pdgId_->Write();
		hist_hepMCstatus_->Write();
		hist_charge_->Write();
		hist_p_->Write();
		hist_pT_->Write();
		hist_eta_->Write();
		hist_phi_->Write();
		hist_daughter_dR_->Write();
		hist_daughter_dEta_->Write();
		hist_daughter_dPhi_->Write();
		hist_daughter_openingAngle_->Write();

		//Move back to the original directory within the output file ...
		ptr_outFile->cd("..");
	}

	std::vector<TH1D*> MCParticleDistns::GetPtrsToHistos(){
		std::vector<TH1D*> tmpPtrsVector;
		tmpPtrsVector.clear();

		tmpPtrsVector.push_back(hist_pdgId_);
		tmpPtrsVector.push_back(hist_hepMCstatus_);
		tmpPtrsVector.push_back(hist_charge_);
		tmpPtrsVector.push_back(hist_p_);
		tmpPtrsVector.push_back(hist_pT_);
		tmpPtrsVector.push_back(hist_eta_);
		tmpPtrsVector.push_back(hist_phi_);
		tmpPtrsVector.push_back(hist_daughter_dR_);
		tmpPtrsVector.push_back(hist_daughter_dEta_);
		tmpPtrsVector.push_back(hist_daughter_dPhi_);
		tmpPtrsVector.push_back(hist_daughter_openingAngle_);

		return tmpPtrsVector;
	}
}

#endif
