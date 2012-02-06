#ifndef tsw_HEEPGsfProducer_h
#define tsw_HEEPGsfProducer_h

#include "DataFormats/Provenance/interface/ParameterSetID.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <string>


namespace edm{
  class Event;
  class EventSetup;
  class ParameterSet;
}
namespace tsw{
	class HEEPGsfProducer : public edm::EDProducer {


	private:
	  edm::InputTag cutValueMapTag_;
	  edm::InputTag inputGsfCollnTag_;


	public:
	  explicit HEEPGsfProducer(const edm::ParameterSet& iPara);
	  ~HEEPGsfProducer(){}

	 private:
	  virtual void beginJob(){}
	  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
	  virtual void endJob(){}

};
}

#endif
