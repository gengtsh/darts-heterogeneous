#ifndef TESTRESETTP_H
#define TESTRESETTP_H

#include <stdint.h>
#include "DARTS.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES
using namespace darts;

DEF_CODELET_ITER(TestResetComputerCD,0,SHORTWAIT);
DEF_CODELET(TestResetControllerCD,2,LONGWAIT);

DEF_TP(TestResetTP)
{
	uint64_t nIter;
	uint64_t nComputer;
	TestResetComputerCD *compute;
	TestResetControllerCD controller;
	Codelet *signalUp;

	TestResetTP(uint64_t nIter,uint64_t nComputer,Codelet *up)
	:nIter(nIter)
	,nComputer(nComputer)
	,compute(new TestResetComputerCD[nComputer])
	,controller(nComputer,nComputer,this,LONGWAIT)
	,signalUp(up)
	{
		for ( size_t i = 0; i < nComputer; ++i ) {
            compute[i] = TestResetComputerCD {0,1,this,SHORTWAIT,i};
            add( compute + i );
        }
	}
	
	virtual	~TestResetTP(){delete []compute;}
};

#endif
