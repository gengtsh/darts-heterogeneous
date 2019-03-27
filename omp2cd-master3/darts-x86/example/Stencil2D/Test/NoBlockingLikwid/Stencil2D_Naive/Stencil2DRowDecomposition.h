#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <stdint.h>
#include "DARTS.h"
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

DEF_CODELET_ITER(Stencil2DRowLoop,0,SHORTWAIT);
DEF_CODELET(Stencil2DRowSyncSwap,2,LONGWAIT);

DEF_TP(Stencil2DRowDecomposition)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t timeStep;
    
	Stencil2DRowLoop *compute;
    Stencil2DRowSyncSwap  syncSwap;
    Codelet *signalUp;
	uint32_t nRowsCut;
	
	Stencil2DRowDecomposition(double *inimatrix, 
                              const uint64_t inim, const uint64_t inin,
                              double *newmatrix, uint64_t ts, 
                              Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	,compute(new Stencil2DRowLoop[N_CORES])
	,syncSwap(N_CORES,N_CORES,this,LONGWAIT)
	,signalUp(up)
    ,nRowsCut(N_CORES)
	{
		add(&syncSwap);
		for ( size_t i = 0; i < nRowsCut; ++i ) {
            compute[i] = Stencil2DRowLoop {0,1,this,SHORTWAIT,i};
            add( compute + i );
        }
//		add(&syncSwap);

	}
	~Stencil2DRowDecomposition(){delete []compute;}
};

#endif
