#ifndef tswEleMuObject_h
#define tswEleMuObject_h

#include "TSWilliams/BstdZeeAnalyser/interface/tswMuon.h"

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
			{ }
			EleMuObject(tsw::HEEPEle ele, tsw::Muon muon):
				ele_(ele), muon_(muon), p4_( muon_.p4()+ele_.p4() )
			{ }
			~EleMuObject(){}

			// Methods to access the electron and muon ...
			const tsw::HEEPEle& ele() const  { return ele_; }
			const tsw::Muon& muon() const  { return muon_; }

			// Methods to access kinematic variables of the ele-mu object ...
			int charge() const { return ( muon_.charge()+ele_.charge() ); }

			const TLorentzVector& p4() const { return p4_; }
			double p()    const { return p4().P();   }
			double pT()   const { return p4().Pt();  }
			double eta()  const { return p4().Eta(); }
			double phi()  const { return p4().Phi(); }
			double mass() const { return p4().M();   }

			double deltaR()   const { return ele_.p4().DeltaR(muon_.p4());}
			double deltaEta() const { return ( ele_.p4().Eta()-muon_.p4().Eta() );}
			double deltaPhi() const { return ele_.p4().DeltaPhi(muon_.p4()); }
			double openingAngle() const { return ele_.p4().Angle(muon_.p4().Vect()); }

			bool isInZMassRange() const { return (mass()<120 && mass()>60); }
	};
}

#endif
