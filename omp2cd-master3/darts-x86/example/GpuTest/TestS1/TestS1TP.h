#ifndef TESTFGTP_H
#define TESTFGTP_H

#include <cstdint>
#include "DARTS.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define GPUMETA 0x4

DEF_CODELET_ITER(ComputeCD,0,GPUMETA);
DEF_CODELET_ITER(SwapCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(TestS1TP)
{
	uint64_t timeStep;
    
	size_t nCut;
	ComputeCD *compute;
    SwapCD *swap;
	SyncCD	sync;
    Codelet *signalUp;
	uint64_t *tSteps;

	TestS1TP( uint64_t ts, Codelet *up)
	:timeStep(ts)
    ,nCut(N_CORES-1) //1 GPU core, (N_CORES-1) CPU cores
	,compute(new ComputeCD[nCut])
	,swap(new SwapCD[nCut])
	,sync(nCut,nCut,this,LONGWAIT)
	,signalUp(up)
	,tSteps(new uint64_t[nCut])
	{
		for ( size_t i = 0; i < nCut; ++i ) {
			tSteps[i]=ts;
            compute[i] = ComputeCD {0,3,this,GPUMETA,i};
			swap[i]= SwapCD{3,3,this,SHORTWAIT,i};
            add( compute + i );
        }



	}
	
	virtual ~TestS1TP(){
		delete []compute;
		delete []swap;
		delete []tSteps;
	}
};

#endif
