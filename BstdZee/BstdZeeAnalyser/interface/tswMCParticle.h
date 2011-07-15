#ifndef tswMCParticle_h
#define tswMCParticle_h

namespace tsw{
	class MCParticle{
		public:
			// Constructors and destructors ...
			MCParticle(Int_t particle_pdgId, Int_t particle_status, Int_t particle_charge, ROOT::Math::XYZTVector fourMomentum);
			MCParticle(){
				particles_pdgId = -999;
				particles_hepMCstatus = -999;
				particles_charge = -999;
				particles_p4.SetPxPyPzE(99999.9, 0.0, 0.0, 99999.9);
			}
			~MCParticle(){}
			// Methods for access to member variables ...
			Int_t pdgId(){return particles_pdgId;}
			Int_t hepMCstatus(){return particles_hepMCstatus;}
			Int_t charge(){return particles_charge;}
			TLorentzVector p4(){return particles_p4;}
			Double_t pT(){return p4().Pt();}
			Double_t eta(){return p4().Eta();}
			Double_t phi(){return p4().Phi();}

		private:
			Int_t particles_pdgId;
			Int_t particles_hepMCstatus;
			Int_t particles_charge;
			TLorentzVector particles_p4;
	};

	MCParticle::MCParticle(Int_t particle_pdgId, Int_t particle_status, Int_t particle_charge, ROOT::Math::XYZTVector fourMomentum):
		particles_pdgId(particle_pdgId),
		particles_hepMCstatus(particle_status),
		particles_charge(particle_charge),
		particles_p4( ConvertToTLorentzVector(&fourMomentum) )
	{
		//Nothing needed here ....
	}
} // End of tsw namespace


#endif
