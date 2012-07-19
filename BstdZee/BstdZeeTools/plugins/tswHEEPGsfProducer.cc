
#include "TSWilliams/BstdZeeTools/interface/tswHEEPGsfProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <iostream>

tsw::HEEPGsfProducer::HEEPGsfProducer(const edm::ParameterSet& iPara) :
	cutValueMapTag_( iPara.getParameter<edm::InputTag>("cutValueMap") ),
   inputGsfCollnTag_( iPara.getParameter<edm::InputTag>("inputGsfEles") )
{
	produces< std::vector<reco::GsfElectron> >();
}


void tsw::HEEPGsfProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	edm::Handle<reco::GsfElectronCollection> h_inputEles;
	iEvent.getByLabel(inputGsfCollnTag_,h_inputEles);

	edm::Handle<edm::ValueMap<int> > h_heepNoIsoBits;
	iEvent.getByLabel(cutValueMapTag_,h_heepNoIsoBits);


	std::auto_ptr< std::vector<reco::GsfElectron> > heepEles( new std::vector<reco::GsfElectron> );
	for(size_t eleNr=0;eleNr<h_inputEles->size();eleNr++){
		reco::GsfElectronRef gsfRef(h_inputEles,eleNr);
		if((*h_heepNoIsoBits)[gsfRef]==0x0) heepEles->push_back(*gsfRef);
	}

	iEvent.put(heepEles);

}

//define this as a plug-in
DEFINE_FWK_MODULE(tsw::HEEPGsfProducer);
