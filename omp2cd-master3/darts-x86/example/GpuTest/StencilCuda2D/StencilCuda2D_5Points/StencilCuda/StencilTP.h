#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H


#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define GPUMETA 0x4

DEF_CODELET(Stencil2D4ptGpuCD,2,LONGWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(StencilTP)
{
    
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t timeStep;
	Stencil2D4ptGpuCD Stencil2D4ptGpu;
	SyncCD	sync;
    Codelet *signalUp;
	uint64_t *tSteps;

	StencilTP(double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t ts,Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	,Stencil2D4ptGpu(0,0,this,GPUMETA)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{
		add(&Stencil2D4ptGpu);
	}
	
	virtual ~StencilTP(){
	}
};

#endif
