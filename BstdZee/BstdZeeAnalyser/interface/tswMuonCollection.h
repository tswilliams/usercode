#ifndef tswMuonCollection_h
#define tswMuonCollection_h

#include "tswMuon.h"

namespace tsw{
	class MuonCollection{
		private:
			std::vector<Muon> vecOfMuons_;
			bool orderedByPt_;
			std::string muonType_;

		public:
			// CTORs and DTOR ...
			MuonCollection(): vecOfMuons_(), orderedByPt_(false), muonType_()
				{}
			MuonCollection(std::vector<Muon>* theMuons, std::string typeOfMuons):
				vecOfMuons_(*theMuons), orderedByPt_(false), muonType_(typeOfMuons)
				{}
			~MuonCollection(){}

			// Method for printing out the information of the muons in this
			void Print();

			// Methods for extracting information about the muons ...
			int NumOfMuons(){ return vecOfMuons_.size(); }
			std::vector<Muon> ReturnVecOfMuons(){ return vecOfMuons_; }
			Muon at(int muonIdx){
				return vecOfMuons_.at(muonIdx);}

			Muon HighestPtMuon(){
				if(orderedByPt_)
					return vecOfMuons_.at(0);
				else{
					MuonCollection tmpOrderedMuColln(&vecOfMuons_, muonType_);
					tmpOrderedMuColln.OrderBypT();
					return tmpOrderedMuColln.at(0);
				}
			}
			Muon SecondHighestPtMuon(){
				if(orderedByPt_)
					return vecOfMuons_.at(1);
				else{
					MuonCollection tmpOrderedMuColln(&vecOfMuons_, muonType_);
					tmpOrderedMuColln.OrderBypT();
					return tmpOrderedMuColln.at(1);
				}
			}

			// Method for ordering the muons by pT ...
			void OrderBypT();

			// Methods for applying cuts to the muon collection ...
			MuonCollection GetMuonsWithinAcc(){
				std::vector<Muon> tmpMuonsVec; tmpMuonsVec.clear();
				for(unsigned int muIdx=0; muIdx<vecOfMuons_.size(); muIdx++){
					if(vecOfMuons_.at(muIdx).isWithinAcc())
						tmpMuonsVec.push_back( vecOfMuons_.at(muIdx) );
				}
				return MuonCollection(&tmpMuonsVec, muonType_+", within Acc");
			}
			MuonCollection GetTightMuons(){
				std::vector<Muon> tightMuonsVec; tightMuonsVec.clear();
				for(unsigned int muIdx=0; muIdx<vecOfMuons_.size(); muIdx++){
					if(vecOfMuons_.at(muIdx).isTightMuon())
						tightMuonsVec.push_back( vecOfMuons_.at(muIdx) );
				}
				return MuonCollection(&tightMuonsVec, muonType_+", afterTightCuts");
			}
	};

	void MuonCollection::Print(){
		std::cout << "  ->There are " << vecOfMuons_.size() << " muons of type '" << muonType_ << "' in this event" << std::endl;
		for(unsigned int muonIdx=0; muonIdx<vecOfMuons_.size(); muonIdx++){
			Muon ithMuon = vecOfMuons_.at(muonIdx);
			// General information ...
			std::cout << "    **Muon " << muonIdx << ": charge=" << ithMuon.charge() << std::endl;
			std::cout << "        p=" << ithMuon.p4().P() << "; (pT, eta, phi) = (" << ithMuon.pT() << ", " << ithMuon.eta() << ", " << ithMuon.phi() << ")" << std::endl;
			std::cout << "        is{Global,Tracker,StandAlone}Muon = {" << ithMuon.isGlobalMuon() << ", " << ithMuon.isTrackerMuon() << ", " << ithMuon.isStandAloneMuon() << "}" << std::endl;
			std::cout << "        numberOfMatchedMuonStations = " << ithMuon.numMatchedMuonStns() << std::endl;
			std::cout << "        isolR03_sumPt = " << ithMuon.isolR03_sumPt() << std::endl;
			// Kin. variables from the global track ...
			std::cout << "       *GlobTrk:" << std::endl;
			if(ithMuon.glob_exists()){
				std::cout << "             (pT, eta, phi) = (" << ithMuon.glob_pT() << ", " << ithMuon.glob_eta() << ", " << ithMuon.glob_phi() << "); charge=" << ithMuon.glob_charge() << std::endl;
				std::cout << "             numValidMuonHits=" << ithMuon.glob_numValidMuonHits() << "; normalisedChi2=" << ithMuon.glob_normalisedChi2() << std::endl;
			}
			else
				std::cout << "           *** global track does NOT 'exist'. ***" << std::endl;
			// ... and from the inner track ...
			std::cout << "       *InnerTrk:" << std::endl;
			if(ithMuon.inner_exists()){
				std::cout << "             (pT, eta, phi) = (" << ithMuon.inner_pT() << ", " << ithMuon.inner_eta() << ", " << ithMuon.inner_phi() << "); charge=" << ithMuon.inner_charge() << std::endl;
				std::cout << "             numValidPixHits=" << ithMuon.inner_numValidPixHits() << "; numValidTrkrHits=" << ithMuon.inner_numValidTrkrHits() << std::endl;
				std::cout << "             d_xy vs origin = " << ithMuon.inner_dxyVsOrigin() << std::endl;
			}
			else
				std::cout << "           *** inner track does NOT 'exist'. ***" << std::endl;
			// ... and from the outer track ...
			std::cout << "       *OuterTrk:" << std::endl;
			if(ithMuon.outer_exists()){
				std::cout << "             (pT, eta, phi) = (" << ithMuon.outer_pT() << ", " << ithMuon.outer_eta() << ", " << ithMuon.outer_phi() << ")" << std::endl;
				std::cout << "             charge=" << ithMuon.outer_charge() << std::endl;
			}
			else
				std::cout << "           *** outer track does NOT 'exist'. ***" << std::endl;
		}
	}

	void MuonCollection::OrderBypT(){
		bool beVerbose = false;
		// Declare variables used to temporarily store i'th and (i+1)'th muon ...
		tsw::Muon ithMuon;
		tsw::Muon iPlusOnethMuon;

		// Declare boolean used as flag for whether any muons were swapped in latest iteration ...
		bool swappedMuons = true;

		if(beVerbose){
			std::cout << "   ***** Muons being ordered:" << std::endl;
			std::cout << "   * Current ordering:" << std::endl << "   *   pT= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).pT();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl << "   *   eta= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).eta();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl << "   *   phi= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).phi();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
		}

		// Can only compare two adjacent muons' pT's if there is more than one muon ...
		if(vecOfMuons_.size()>1){
			while(swappedMuons==true){
				//Reset muon-swap boolean to false ...
				swappedMuons = false;

				// Go through all elements of the vector, looking at the pairs of muons that are adjacent under the current order ...
				for(unsigned int iMuon=0; iMuon < vecOfMuons_.size()-1 ; iMuon++){
					ithMuon = vecOfMuons_.at(iMuon);
					iPlusOnethMuon = vecOfMuons_.at(iMuon+1);

					// If the i'th muon has less pT than the (i+1)'th muon then swap these two muons so that they are in descending E_T order ...
					if(ithMuon.pT()<iPlusOnethMuon.pT()){
						swappedMuons = true;
						vecOfMuons_.at(iMuon)   = iPlusOnethMuon;
						vecOfMuons_.at(iMuon+1) = ithMuon;
					}
				}

				if(beVerbose){
					std::cout << "   * Reordering =>" << std::endl << "   *   pT= ";
					for(unsigned int i=0; i<vecOfMuons_.size(); i++){
						std::cout << vecOfMuons_.at(i).pT();
						if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
					}
					std::cout << std::endl;
				}

			}
		}

		if(beVerbose){
			std::cout << "   * Final ordering:" << std::endl;
			std::cout << "   *   pT= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).pT();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl << "   *   eta= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).eta();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl << "   *   phi= ";
			for(unsigned int i=0; i<vecOfMuons_.size(); i++){
				std::cout << vecOfMuons_.at(i).phi();
				if(i!=( vecOfMuons_.size()-1 )){ std::cout << ", "; }
			}
			std::cout << std::endl;
			std::cout << "   ***** Muons should now be ordered in descending Et" << std::endl;
		}
		// Set the flag to denote that muons are now ordered by pT ...
		orderedByPt_ = true;
	}


} // End of tsw namespace

#endif
