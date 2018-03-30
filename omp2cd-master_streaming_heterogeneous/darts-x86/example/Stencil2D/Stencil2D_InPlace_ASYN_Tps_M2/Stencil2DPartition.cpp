#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "Stencil2DPartition.h"
#include "Stencil2DRowDecomposition.h"
#include <math.h>
#include <inttypes.h>

//#include <pthread.h>
//extern pthread_mutex_t mutex;
//#include <sstream>
//#include <iostream>

#include <cassert>

void Stencil2DPartitionChunks::fire(void) 
{
	LOAD_FRAME(Stencil2DPartition);

//	pthread_mutex_lock(&mutex);
//	std::cout<<"PartitionChunks"<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"PartitonChunks"<<"\n";
//	std::cout <<ss.str();
	
	double *initial = FRAME(initial); //matrix pointer initial Matrix[M][n]
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	uint64_t ts=FRAME(timeStep);	
	uint64_t nTpRows = n_rows-2;
	uint64_t nChunks=nTpRows/g_nSU;
	
	for(size_t i=0; i<g_nSU; ++i){
		uint64_t pos = nChunks*i*n_cols;
		uint64_t nTpRowsFinal = ((i==(g_nSU-1))? (nTpRows%g_nSU):0) + nChunks;
		INVOKE(Stencil2DRowDecomposition,initial+pos,nTpRowsFinal+2,n_cols,ts,i,FRAME(checkTP), &FRAME(sync),FRAME(check[i])); 
		
	}


	EXIT_TP();

}



void Stencil2DPartitionCheckTP::fire(void){
	LOAD_FRAME(Stencil2DPartition);
	uint64_t Id = getID();
	RESET(checkTP[Id]);	
	
//	pthread_mutex_lock(&mutex);
//	std::cout<<"CheckTP"<<Id<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"D"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
	
	if(Id==0){
		SYNC(check[Id][0]);
	}else if(Id == g_nSU){
		SYNC(check[Id-1][1]);
	}else{
		SYNC(check[Id-1][1]);
		SYNC(check[Id][0]);
	}

	EXIT_TP();
}


void
Stencil2DPartitionSync::fire(void)
{
    LOAD_FRAME(Stencil2DPartition);

	ThreadedProcedure* PTP = getTP();
	PTP->setRef(1);
	SIGNAL(signalUp);
	//ThreadedProcedure* PTP = getTP();
	//while((PTP->getRef())!=1){
	//	PTP->decRef();
	//}

//	pthread_mutex_lock(&mutex);
//	std::cout<<"PartitionSync!"<<std::endl;
//	pthread_mutex_unlock(&mutex);
    
	EXIT_TP();
}

