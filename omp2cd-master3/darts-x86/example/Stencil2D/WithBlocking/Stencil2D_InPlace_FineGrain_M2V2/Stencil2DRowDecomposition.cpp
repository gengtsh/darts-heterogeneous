#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "Stencil2DRowDecomposition.h"
#include "Stencil2DKernel.h"
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

void 
Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);

	uint64_t Id	= getID();	

	if((FRAME(tSteps)[Id])==0){
		SYNC(sync);
		EXIT_TP();
	}
	
	RESET(compute[Id]);
	uint64_t nRowsCut = FRAME(nRowsCut);	
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	double *dst=FRAME(dstPtr[Id]);
	double *src=FRAME(srcPtr[Id]);
	uint64_t chunk = FRAME(chunk);
	uint64_t rpos2 = (Id==(nRowsCut-1))? (n_rows-Id*chunk):(chunk+2) ;
	
	//computeInner_stencil2d_v2(dst,src,rpos2,n_cols);
  
    computeInner_stencil2dbkv3(dst,src,rpos2,n_cols,n_rows,n_cols);

	SYNC(syncSwap[Id]);
	
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


	EXIT_TP();
}

void
Stencil2DRowSyncSwap::fire(void)
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	uint64_t Id = getID();	
	//std::cout<<"core num: "<<sched_getcpu()<<" swap: "<<Id <<" timestep: "<<FRAME(tSteps)[Id]<<"\n"<<std::endl;		
	
	//std::cout<<" swap: "<<Id <<" timestep: "<<FRAME(tSteps)[Id]<<"\n"<<std::endl;		
	
	uint64_t *ts = FRAME(tSteps);
	if(ts[Id]==1){
		ts[Id]=ts[Id]-1;
		SYNC(sync);	
		EXIT_TP();
	}

	ts[Id]=ts[Id]-2;
	
	RESET(syncSwap[Id]);
	double *dst=FRAME(dstPtr[Id]);
	double *src=FRAME(srcPtr[Id]);
	//SWAP_PTR(&dst,&src);
	
	uint64_t nRowsCut = FRAME(nRowsCut);	
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	uint64_t chunk = FRAME(chunk);
	uint64_t rpos2 = (Id==(nRowsCut-1))? (n_rows-Id*chunk):chunk+2 ;
	
	//computeInner_stencil2d_v2(src,dst,rpos2,n_cols);
	
    computeInner_stencil2dbkv3(src,dst,rpos2,n_cols,n_rows,n_cols);

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
	
	EXIT_TP();
	
}


void
Stencil2DRowSync::fire(void)
{
    LOAD_FRAME(Stencil2DRowDecomposition);

	ThreadedProcedure* PTP = getTP();
	PTP->setRef(1);
	
	SIGNAL(signalUp);
    EXIT_TP();
}

