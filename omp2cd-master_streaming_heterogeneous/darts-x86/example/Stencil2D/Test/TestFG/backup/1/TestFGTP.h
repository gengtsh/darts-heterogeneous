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
	uint64_t timeStep;
    
	ComputeCD *compute;
    SwapCD *swap;
	SyncCD	sync;
    Codelet *signalUp;

	uint64_t nCut;
	uint64_t *tSteps;	

	TestFGTP(uint64_t ts, Codelet *up)
	:timeStep(ts)
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


	}
	virtual ~TestFGTP(){
		delete []compute;
		delete []swap;
		delete []tSteps;
	}
};

#endif
