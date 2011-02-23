// -*- C++ -*-
//
// Package:    EvtKinAnr
// Class:      EvtKinAnr
// 
/**\class EvtKinAnr EvtKinAnr.cc Kinematics/EvtKinAnr/src/EvtKinAnr.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams
//         Created:  Tue Feb 22 16:28:14 GMT 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//To use the GenParticle MC-truth collection...
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//...for histograms creation
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

//
// class declaration
//

class EvtKinAnr : public edm::EDAnalyzer {
   public:
      explicit EvtKinAnr(const edm::ParameterSet&);
      ~EvtKinAnr();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
		TH1D *numZs_status2Hist;
		TH1D *numZs_status3Hist;
		TH1D *numZs_statusOtherHist;
		TH1D *status2Zs_PtHist;
		TH1D *eleHighestEt_EtHist;
		TH1D *ele2ndHighestEt_EtHist;
		TH1D *numAboveSingleEleEtThreshold_Hist;
		edm::Service<TFileService> fHistos;

		int numEvts;

		double histoEtMin;
		double histoEtMax;
		int    histoEtNBins;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EvtKinAnr::EvtKinAnr(const edm::ParameterSet& iConfig)

{
	//now do what ever initialization is needed
	numEvts = 0;

	histoEtMin = 0.0;
	histoEtMax = 80.0;
	histoEtNBins = 40;
	double histoEtBinWidth = (histoEtMax-histoEtMin)/static_cast<double>(histoEtNBins);
	
	numZs_status2Hist      = fHistos->make<TH1D>("numZs_status2Hist",     "Number of Z bosons with status=2 in each event", 8, -0.5, 7.5);
	numZs_status3Hist      = fHistos->make<TH1D>("numZs_status3Hist",     "Number of Z bosons with status=3 in each event", 8, -0.5, 7.5);
	numZs_statusOtherHist  = fHistos->make<TH1D>("numZs_statusOtherHist", "Number of Z bosons with status!=2or3 in each event", 8, -0.5, 7.5);
	status2Zs_PtHist       = fHistos->make<TH1D>("status2Zs_PtHist",      "p_{T} distribution of decaying Z bosons; Z boson p_{T} /GeVc^{-1}; Number per 2GeVc^{-1}",          histoEtNBins, histoEtMin, histoEtMax);
	eleHighestEt_EtHist    = fHistos->make<TH1D>("eleHighestEt_EtHist",   "p_{T} distribution of highest p_{T} electron; Electron p_{T} /GeVc^{-1}; Number per 2GeVc^{-1}",    histoEtNBins, histoEtMin, histoEtMax);
	ele2ndHighestEt_EtHist = fHistos->make<TH1D>("ele2ndHighestEt_EtHist","p_{T} distribution of 2nd highest p_{T} electron; Electron p_{T} /GeVc^{-1}; Number per 2GeVc^{-1}",histoEtNBins, histoEtMin, histoEtMax);
	numAboveSingleEleEtThreshold_Hist = fHistos->make<TH1D>("numAboveSingleEleEtThreshold_Hist", "Fraction of events passing a single electron Et threshold; Threshold E_{T} /GeV; Fraction", histoEtNBins, histoEtMin-histoEtBinWidth/2.0, histoEtMax-histoEtBinWidth/2.0);

   numAboveSingleEleEtThreshold_Hist->SetLineColor(4);
   numAboveSingleEleEtThreshold_Hist->SetMarkerColor(4);
   numAboveSingleEleEtThreshold_Hist->SetMarkerStyle(20);
   numAboveSingleEleEtThreshold_Hist->SetMarkerSize(0.4);
   numAboveSingleEleEtThreshold_Hist->GetXaxis()->CenterTitle();
   numAboveSingleEleEtThreshold_Hist->GetYaxis()->CenterTitle();
	
}


EvtKinAnr::~EvtKinAnr()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
EvtKinAnr::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//Variable declarations...
	int numZs_status2 = 0;
	int numZs_status3 = 0;
	int numZs_statusOther = 0;
	int particlePDGid = -999;
	int particleStatus = -999;
	reco::Candidate::LorentzVector particle4mom(0.0, 0.0, 0.0, 0.0);
	unsigned int numDaughters = -999;
	//const reco::GenParticle mcParticle;
	int idx_eleHighestEt = -999;
	double eleHighestEt = -999.9; //Must be set to less than 0.0
	int idx_ele2ndHighestEt = -999; 
	double ele2ndHighestEt = -9999.9; //Should be set to less than eleHighestEt 


   using namespace edm;
	//Handles, and grabbing various CMSSW-specific data types from event...
	Handle<reco::GenParticleCollection> genParticles;
	iEvent.getByLabel("genParticles", genParticles);

	//std::cout << "Analysing event no. " << iEvent.id() << " ..." << std::endl;

	numEvts++;
	std::cout << "   There are " << genParticles->size() << " generator-particles in this event." << std::endl;
	
	//Running over all of the particles in the MC truth...
	for(unsigned int idx=0; idx<genParticles->size(); idx++){
		//Extracting the properties of each of the MC truth particles....
		const reco::GenParticle & mcParticle = (*genParticles)[idx];
		particlePDGid = mcParticle.pdgId();
		particleStatus = mcParticle.status();
		particle4mom = mcParticle.p4();
		numDaughters = mcParticle.numberOfDaughters();
		//Finding the Z bosons, and adding to various counts...
		if(particlePDGid==23){
			std::cout << "   Generator-particle no. " << idx << " is a Z boson. (pdgId= " << particlePDGid << "; status=" << particleStatus << ")" << std::endl;
			std::cout << "      px=" << particle4mom.Px() << "GeV?; py=" << particle4mom.Py() << "GeV?; pz=" << particle4mom.Pz() << "GeV?" << std::endl;
			std::cout << "      pt=" << particle4mom.Pt() << "GeV?; eta=" << particle4mom.Eta() << "; phi=" << particle4mom.Phi() << std::endl;
			std::cout << "      It has " << numDaughters << " daughters." << std::endl;
			//For Z bosons with status=2 (i.e. that are decaying)
			if(particleStatus==2){
				numZs_status2++;
				status2Zs_PtHist->Fill(particle4mom.Pt());
				for(size_t j=0; j<numDaughters; j++){
       			const reco::Candidate * mcDaughter =mcParticle.daughter(j);
       			std::cout << "         Daughter " << j << ": pdgId=" << mcDaughter->pdgId() << "; pT=" << mcDaughter->p4().Pt() << std::endl;
     			}
			}
			//For Z bosons with status=3 (i.e. that are in the 'hard part' of the interaction - i.e. used in the matrix element caLculation)
			else if(particleStatus==3){
				numZs_status3++;
				for(size_t j=0; j<numDaughters; j++){
       			const reco::Candidate * mcDaughter =mcParticle.daughter(j);
       			std::cout << "         Daughter " << j << ": pdgId=" << mcDaughter->pdgId() << "; pT=" << mcDaughter->p4().Pt() << std::endl;
     			}
			}
			//For Z bosons with other status numbers
			else
				numZs_statusOther++;
		}
		//Finding the electrons and positrons...
		if(abs(particlePDGid)==11){
			if(particlePDGid==11)
				std::cout << "   Generator-particle no. " << idx << " is an electron. (pdgId= " << particlePDGid << "; status=" << particleStatus << "; Pt=" << particle4mom.Pt() << ")" << std::endl;
			else
				std::cout << "   Generator-particle no. " << idx << " is a positron.  (pdgId= " << particlePDGid << "; status=" << particleStatus << "; Pt=" << particle4mom.Pt() << ")" << std::endl;
			
			if(particleStatus!=3){
				if(particle4mom.Pt()>eleHighestEt){ //If electron has greater Pt than any electron before it...
					idx_ele2ndHighestEt = idx_eleHighestEt; 	// Switch the old highest Pt electron to be the 2nd highest Pt ele...
					ele2ndHighestEt = eleHighestEt;
					idx_eleHighestEt = idx; 						// .. and set this ele as the highest Pt ele.
					eleHighestEt = particle4mom.Pt();
					std::cout << "      (Now set as highest Pt electron...)" << std::endl;
					std::cout << "      (... and #" << idx_ele2ndHighestEt << " set as 2nd highest Pt ele...)" << std::endl;
				}
				else if(particle4mom.Pt()>ele2ndHighestEt){ 	// Otherwise, if ele Pt is still larger than previous 2nd highest ele Pt...
					idx_ele2ndHighestEt = idx;						// ...set this ele as the 2nd highest PT ele.
					ele2ndHighestEt = particle4mom.Pt();
					std::cout << "      (Now set as 2nd highest Pt electron...)" << std::endl;
				}
			}
		}
	}
	std::cout << " ->There are " << numZs_status2 << " Z bosons with status=2 in this event." << std::endl;
	numZs_status2Hist->Fill(numZs_status2);
	std::cout << "   There are " << numZs_status3 << " Z bosons with status=3 in this event." << std::endl;
	numZs_status3Hist->Fill(numZs_status3);
	std::cout << "   There are " << numZs_statusOther << " Z bosons with other status numbers in this event." << std::endl;
	numZs_statusOtherHist->Fill(numZs_statusOther);

	std::cout << " ->The highest Et electron in the event was genParticle#" << idx_eleHighestEt << ", with Et=" << eleHighestEt << std::endl;
	eleHighestEt_EtHist->Fill(eleHighestEt);
	std::cout << "   The 2nd highest Et electron in the event was genParticle#" << idx_ele2ndHighestEt << ", with Et=" << ele2ndHighestEt << std::endl;
	ele2ndHighestEt_EtHist->Fill(ele2ndHighestEt);
	
	std::cout << std::endl;

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
EvtKinAnr::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
EvtKinAnr::endJob() {
	int totalNumHighestEtEles = numEvts;
	double dbl_totalNumHighestEtEles = static_cast<double>(totalNumHighestEtEles); 
	double dbl_numHighestEtElesAboveThreshold = -999.9;
	double effiOfEtThreshold = -999.9;
	double error_effiOfEtThreshold = -999.9;
	
	//Setting errors on histogram bins to be sum of weights for bin entries (i.e. Poisson errors)
	status2Zs_PtHist->Sumw2();
	eleHighestEt_EtHist->Sumw2();
	ele2ndHighestEt_EtHist->Sumw2();
	
	//Calculating the fraction of events passing a single electron Et threshold for the numAboveSingleEleEtThreshold histogram...
	std::cout << std::endl << "Program is now calculating the fraction of events passing a single electron Et threshold..." << std::endl;
	for(int binIdx=0; binIdx<=histoEtNBins+1; binIdx++){
		std::cout << "  *binIdx=" << binIdx << std::endl;
		dbl_numHighestEtElesAboveThreshold = eleHighestEt_EtHist->Integral(binIdx,histoEtNBins+1);
		effiOfEtThreshold = dbl_numHighestEtElesAboveThreshold/dbl_totalNumHighestEtEles;
		error_effiOfEtThreshold = sqrt(effiOfEtThreshold*(1.0-effiOfEtThreshold)/dbl_totalNumHighestEtEles);
		std::cout << "     dbl_numHighestEtElesAboveThreshold = " << dbl_numHighestEtElesAboveThreshold << std::endl;
		std::cout << "     dbl_totalNumHighestEtEles = " << dbl_totalNumHighestEtEles << std::endl;
		std::cout << "     effiOfEtThreshold = " << effiOfEtThreshold << std::endl;
		
		numAboveSingleEleEtThreshold_Hist->SetBinContent(binIdx,effiOfEtThreshold);
		numAboveSingleEleEtThreshold_Hist->SetBinError(binIdx,error_effiOfEtThreshold);
	}
	
	std::cout << std::endl <<  "**histoEtBinWidth=" << (histoEtMax-histoEtMin)/static_cast<double>(histoEtNBins) << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(EvtKinAnr);
