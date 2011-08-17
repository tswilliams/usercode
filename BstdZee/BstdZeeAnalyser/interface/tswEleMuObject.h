#ifndef tswEleMuObject_h
#define tswEleMuObject_h

#include "BstdZeeFirst/Analyser/interface/tswMuon.h"

namespace tsw{
	class EleMuObject{
		private:
			tsw::HEEPEle ele_;
			tsw::Muon muon_;
			TLorentzVector p4_;

		public:
			// CTORs and DTOR ...
			EleMuObject() :
				ele_(), muon_(), p4_()
				{}
			EleMuObject(tsw::HEEPEle ele, tsw::Muon muon):
				ele_(ele), muon_(muon), p4_( muon_.p4()+ele_.p4() )
				{}
			~EleMuObject(){}

			// Methods to access the electron and muon ...
			tsw::HEEPEle* GetElectron(){
				return &ele_;}
			tsw::Muon* GetMuon(){
				return &muon_;}

			// Methods to access kinematic variables of the ele-mu object ...
			int charge(){
				return ( muon_.charge()+ele_.charge() );}
			TLorentzVector p4(){
				return p4_;}
			double p(){
				return p4().P();}
			double pT(){
				return p4().Pt();}
			double eta(){
				return p4().Eta();}
			double phi(){
				return p4().Phi();}
			double mass(){
				return p4().M();}

			double deltaR(){
				return ele_.p4().DeltaR(muon_.p4());}
			double deltaEta(){
				return ( ele_.p4().Eta()-muon_.p4().Eta() );}
			double deltaPhi(){
				return ele_.p4().DeltaPhi(muon_.p4());}
			double openingAngle(){
				return ele_.p4().Angle(muon_.p4().Vect());}

			bool isInZMassRange(){
				if( mass()<120 && mass()>60)
					return true;
				else
					return false;
			}

	};
}

#endif
