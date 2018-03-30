#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "HelloTP.h"
#include "Hello_cuda.h"
#include <cassert>
//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
#include <iostream>


void 
InvokeCudaCD::fire(void) 
{
	LOAD_FRAME(HelloTP);

	std::cout<<"Invoke Cuda"<<std::endl;	

	Hello_cuda();

	SYNC(sync);
	EXIT_TP();
}


void
SyncCD::fire(void)
{
    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(HelloTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





