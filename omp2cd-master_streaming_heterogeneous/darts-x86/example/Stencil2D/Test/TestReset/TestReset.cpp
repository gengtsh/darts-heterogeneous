#define __STDC_FORMAT_MACROS
#include <cstdint>
//#include <cmath>
#include <inttypes.h>
#include "TestReset.h"
#include <sstream>
#include <iostream>

void 
TestResetComputerCD::fire(void) 
{
	LOAD_FRAME(TestResetTP);

    uint64_t Id	    = getID();	
	uint64_t nIter  = FRAME(nIter);
	if(nIter > 0 ){
		RESET(compute[Id]);
	}
	std::stringstream ss;
	ss<<"C"<<Id<<"],ts:"<<nIter<<"\n";
	std::cout<<ss.str();

	SYNC(controller);

	EXIT_TP();
}

void
TestResetControllerCD::fire(void)
{
	LOAD_FRAME(TestResetTP);

    uint64_t Id	    = getID();	
	uint64_t nIter  = FRAME(nIter);

	std::stringstream ss;
	ss<<"S"<<Id<<"],ts:"<<nIter<<"\n";
	std::cout<<ss.str();
	
	if ( FRAME(nIter)-- > 0 ) {
		RESET(controller);
        uint64_t nComputer = FRAME(nComputer);	
		for(size_t i = 0; i < nComputer; ++i)
			SYNC(compute[i]);
    } else {
		SIGNAL(signalUp);
    }



	EXIT_TP();
}
