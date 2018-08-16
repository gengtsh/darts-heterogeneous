#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "Stencil2DRowDecomposition.h"
#include "Stencil2DKernel.h"
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
//#include <iostream>

void 
Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);

	uint64_t Id	= getID();	

	RESET(compute[Id]);
	
	uint64_t *ts = FRAME(tSteps);
	
	if(ts[Id]==0){
		SYNC(sync);
		EXIT_TP();
	}
	
//	pthread_mutex_lock(&mutex);
//	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"C"<<Id<<",ts:"<<ts[Id]<<"\n"<<std::endl;
//	std::cout <<ss.str();

	uint64_t nRowsCut = FRAME(nRowsCut);	
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	double *dst=FRAME(dstPtr[Id]);
	double *src=FRAME(srcPtr[Id]);
	uint64_t chunk = FRAME(chunk);
	uint64_t rpos2 = (Id==(nRowsCut-1))? (n_rows-2-chunk*(N_CORES-1)):chunk ;
	
	computeInner_stencil2d_v2(dst,src,rpos2,n_cols);
   
	
	//std::cout<<"core num: "<<sched_getcpu()<<" compute: "<<Id<< " timestep: "<<FRAME(tSteps)[Id] <<"\n"<<std::endl;		

	//std::cout<<" compute: "<<Id<< " timestep: "<<FRAME(tSteps)[Id] <<"\n"<<std::endl;		

	if(nRowsCut==1){
		SYNC(syncSwap[Id]);
		SYNC(syncSwap[Id]);
	}else{
		if(Id ==0){
			SYNC(syncSwap[Id]);
			SYNC(syncSwap[Id+1]);
		}else if(Id==(nRowsCut-1)){
			SYNC(syncSwap[Id]);
			SYNC(syncSwap[Id-1]);
		}else {
			SYNC(syncSwap[Id+1]);
			SYNC(syncSwap[Id-1]);
		}
	}

	SYNC(syncSwap[Id]);


	EXIT_TP();
}

void
Stencil2DRowSyncSwap::fire(void)
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	uint64_t Id = getID();	
	RESET(syncSwap[Id]);
	uint64_t *ts = FRAME(tSteps);

//	pthread_mutex_lock(&mutex);
//	std::cout <<"S"<<Id<<",ts:"<<ts[Id]<<std::endl;
//	pthread_mutex_unlock(&mutex);

//	std::stringstream ss;
//	ss <<"S"<<Id<<",ts:"<<ts[Id]<<"\n"<<std::endl;
//	std::cout<<ss.str();
	
	//std::cout<<"core num: "<<sched_getcpu()<<" swap: "<<Id <<" timestep: "<<FRAME(tSteps)[Id]<<"\n"<<std::endl;		
	
	
	double *dst=FRAME(dstPtr[Id]);
	double *src=FRAME(srcPtr[Id]);
	SWAP_PTR(&dst,&src);
	
	--(ts[Id]);

	
	uint64_t nRowsCut = FRAME(nRowsCut);	
	if(nRowsCut==1){
		SYNC(compute[Id]);
		SYNC(compute[Id]);
	}else{
		if(Id ==0){
			SYNC(compute[Id]);
			SYNC(compute[Id+1]);
		}else if(Id==(nRowsCut-1)){
			SYNC(compute[Id]);
			SYNC(compute[Id-1]);
		}else {
			SYNC(compute[Id+1]);
			SYNC(compute[Id-1]);
		}
	}		

	SYNC(compute[Id]);
	
//	if(ts[Id]==0){
//		SYNC(sync);
//		//std::cout <<"Z"<<Id<<std::endl;
//	}

	EXIT_TP();
	
}


void
Stencil2DRowSync::fire(void)
{
//    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(Stencil2DRowDecomposition);
	
	ThreadedProcedure* PTP = getTP();
	PTP->setRef(1);
	
	SIGNAL(signalUp);
    EXIT_TP();
}

