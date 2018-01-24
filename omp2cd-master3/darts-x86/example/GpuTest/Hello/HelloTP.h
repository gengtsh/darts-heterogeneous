#ifndef TESTFGTP_H
#define TESTFGTP_H

#include <cstdint>
#include "DARTS.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define GPUMETA 0x4

DEF_CODELET(InvokeCudaCD,2,LONGWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(HelloTP)
{
    
	InvokeCudaCD InvokeCuda;
	SyncCD	sync;
    Codelet *signalUp;
	uint64_t *tSteps;

	HelloTP( Codelet *up)
	:InvokeCuda(0,0,this,GPUMETA)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{
		add(&InvokeCuda);
	}
	
	virtual ~HelloTP(){
	}
};

#endif
