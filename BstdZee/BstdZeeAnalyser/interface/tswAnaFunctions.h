#ifndef tswAnaFunctions_h
#define tswAnaFunctions_h

namespace tsw{
	/// Returns number of instances of true in vector of bools
	unsigned int NumPassingCuts(std::vector<bool> passCutsFlags);
}





// -------------------------------------------------------------------//
//            -->            IMPLEMENTATION            <--            //
// -------------------------------------------------------------------//

unsigned int tsw::NumPassingCuts(std::vector<bool> passCutsFlags){
	unsigned int numberElesPassingCuts = 0;

	//Loop over the flags and count how many occurences of "true" there are ...
	for(unsigned int iEle=0; iEle<passCutsFlags.size(); iEle++){
		if(passCutsFlags.at(iEle))
			numberElesPassingCuts++;
	}

	return numberElesPassingCuts;
}


#endif
