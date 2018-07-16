#ifndef DARTS_STENCIL2DPartition_H
#define DARTS_STENCIL2DPartition_H


#include "DARTS.h"
//#include "Stencil2DRowDecomposition.h"
#include <stdint.h>


using namespace darts;

DEF_CODELET(Stencil2DPartitionChunks,0,SHORTWAIT);
DEF_CODELET(Stencil2DPartitionSyncSwap,1,LONGWAIT);

DEF_TP(Stencil2DPartition){
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New; //matrix pointer New matrix[M][N]
	uint64_t timeStep;
	Stencil2DPartitionChunks nChunks;
	Stencil2DPartitionSyncSwap	syncSwap;
	Codelet *signalUP;
	Codelet **compute;	

	Stencil2DPartition( double *inimatrix,const uint64_t m,const uint64_t n,double *newmatrix,uint64_t ts, Codelet *up)
	:Initial(inimatrix)
	,nRows(m)
	,nCols(n)
	,New(newmatrix)
	,timeStep(ts)
	,nChunks(0,1,this,SHORTWAIT)
	,syncSwap(g_nSU,g_nSU,this,LONGWAIT)
	,signalUP(up)
	,compute(new Codelet*[g_nSU])
	{
		add(&nChunks);
		memset(compute,0,sizeof(Codelet*)*g_nSU);
	//	compute=new Codelet*[g_nSU];
	//	for(size_t i=0; i<g_nSU;i++){
	//		compute[i]=0;
	//	}
	}		


};





#endif

