#ifndef TESTFGTP_H
#define TESTFGTP_H

#include <cstdint>
#include "DARTS.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES


DEF_CODELET_ITER(ComputeCD,0,SHORTWAIT);
DEF_CODELET_ITER(SwapCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(TestFGTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t timeStep;
    
	ComputeCD *compute;
    SwapCD *swap;
	SyncCD	sync;
    Codelet *signalUp;

	uint64_t nCut;
	uint64_t *tSteps;	
	uint64_t chunk;
	
	double **dstPtr;
	double **srcPtr;
	
	TestFGTP( double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,  uint64_t ts, Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	,compute(new ComputeCD[N_CORES])
	,swap(new SwapCD[N_CORES])
	,sync(N_CORES,N_CORES,this,LONGWAIT)
	,signalUp(up)
    ,nCut(N_CORES)
	,tSteps(new uint64_t[N_CORES])
	{
		for ( size_t i = 0; i < nCut; ++i ) {
			tSteps[i]=ts;
            compute[i] = ComputeCD {0,3,this,SHORTWAIT,i};
			swap[i]= SwapCD{3,3,this,SHORTWAIT,i};
            add( compute + i );
        }


		dstPtr  = new double* [N_CORES];
		srcPtr  = new double* [N_CORES];

		chunk= (nRows-2)/nCut;
		for(size_t i=0;i<nCut;++i){
			size_t pos1=1+chunk*i*nCols+nCols;
			dstPtr[i] = Initial+pos1;
			srcPtr[i] = New+pos1;
		}
	}
	
	virtual ~TestFGTP(){
		delete []compute;
		delete []swap;
		delete []tSteps;
		delete []dstPtr;
		delete []srcPtr;
	}
};

#endif
