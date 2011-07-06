
#include "NTupler/BstdZeeNTupler/interface/tswEvent.h"
#include <iostream>

namespace tsw{

	Event::Event() :
		runNum_(0),
		lumiSec_(0),
		evtNum_(0)
	{

	}
	Event::~Event() {}

	void Event::SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber){
		runNum_  = runNumber;
		lumiSec_ = lumiSection;
		evtNum_  = eventNumber;
	}
	void Event::PrintBasicEventInformation(){
		std::cout << "Run " << runNum_ << ", LumiSec " << lumiSec_ << ", event " << evtNum_ << std::endl;
	}
}


#if !defined(__CINT__)
  ClassImp(tsw::Event)
#endif
