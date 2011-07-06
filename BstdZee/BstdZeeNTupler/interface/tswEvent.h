
#ifndef tswEvent_h
#define tswEvent_h

#include "TObject.h"

#include <vector>

namespace tsw{
	class Event : public TObject {
		public:
			Event();
			~Event();

			//Methods to set information ...
			void SetBasicEventInformation(unsigned int runNumber, unsigned int lumiSection, unsigned int eventNumber);
			void PrintBasicEventInformation();

		private:
			//General event information ...
			unsigned int runNum_;
			unsigned int lumiSec_;
			unsigned int evtNum_;

			ClassDef(tsw::Event,2);
	};
}

#endif
