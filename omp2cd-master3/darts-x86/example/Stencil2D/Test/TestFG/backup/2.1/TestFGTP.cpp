#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "TestFGTP.h"
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include <pthread.h>
pthread_mutex_t mutex;
//#include <sstream>
//#include <iostream>

void 
ComputeCD::fire(void) 
{
	LOAD_FRAME(TestFGTP);

	uint64_t Id	= getID();	

	RESET(compute[Id]);
	
	uint64_t *ts = FRAME(tSteps);
	
	if(ts[Id]==0){
		EXIT_TP();
	}
	
	pthread_mutex_lock(&mutex);
	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<std::endl;
	pthread_mutex_unlock(&mutex);

//	std::stringstream ss;
//	ss <<"C"<<Id<<",ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();


	uint64_t nCut = FRAME(nCut);

	if(nCut==1){
		SYNC(swap[Id]);
		SYNC(swap[Id]);
	}else{
		if(Id ==0){
			SYNC(swap[Id]);
			SYNC(swap[Id+1]);
		}else if(Id==(nCut-1)){
			SYNC(swap[Id]);
			SYNC(swap[Id-1]);
		}else {
			SYNC(swap[Id+1]);
			SYNC(swap[Id-1]);
		}
	}

	SYNC(swap[Id]);


	EXIT_TP();
}

void
SwapCD::fire(void)
{
	LOAD_FRAME(TestFGTP);
	uint64_t Id = getID();	
	
	RESET(swap[Id]);
	uint64_t *ts = FRAME(tSteps);

	pthread_mutex_lock(&mutex);
	std::cout <<"S"<<Id<<",ts:"<<ts[Id]<<std::endl;
	pthread_mutex_unlock(&mutex);

//	std::stringstream ss;
//	ss <<"S"<<Id<<",ts:"<<ts[Id]<<"\n";
//	std::cout<<ss.str();
	
	--(ts[Id]);
	
	uint64_t nCut = FRAME(nCut);	
	if(nCut==1){
		SYNC(compute[Id]);
		SYNC(compute[Id]);
	}else{
		if(Id ==0){
			SYNC(compute[Id]);
			SYNC(compute[Id+1]);
		}else if(Id==(nCut-1)){
			SYNC(compute[Id]);
			SYNC(compute[Id-1]);
		}else {
			SYNC(compute[Id+1]);
			SYNC(compute[Id-1]);
		}
	}		

	SYNC(compute[Id]);
	
	if(ts[Id]==0){
		SYNC(sync);
		//std::cout <<"Z"<<Id<<std::endl;
	}

	EXIT_TP();
	
}


void
SyncCD::fire(void)
{
    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(TestFGTP);
	SIGNAL(signalUp);
    EXIT_TP();
}

