#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <stdint.h>
#include "DARTS.h"
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

using namespace darts;

DEF_CODELET_ITER(Stencil2DRowLoop,0,SHORTWAIT);
DEF_CODELET(Stencil2DRowSync,2,LONGWAIT);
//DEF_CODELET(Stencil2DSyncToMaster,1,LONGWAIT);

DEF_TP(Stencil2DRowDecomposition)
{
	uint32_t nTp;
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
    uint64_t *timeStep;
	Stencil2DRowLoop *compute;
    Stencil2DRowSync  sync;
	//Stencil2DSyncToMaster syncToMaster;
	//bool *reset;	
	//uint64_t *counter;
	uint32_t nRowsCut;
	Codelet*   signalUp;

	Stencil2DRowDecomposition(uint32_t ntp, double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t *ts,Codelet *up, Codelet *cpt=0)
	:nTp(ntp)
	,Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	,compute(static_cast<Stencil2DRowLoop*>(cpt?cpt:new Stencil2DRowLoop[g_nCU]))
	,sync(g_nCU,g_nCU,this,LONGWAIT)
	,signalUp(up)
	{

//		nRowsCut	= computeRowDecomposition(inim,inin);
//		compute		= static_cast<Stencil2DRowLoop*>(cpt?cpt:(new Stencil2DRowLoop[nRowsCut]));
//		sync		= Stencil2DRowSync{nRowsCut,nRowsCut,this,LONGWAIT};
		nRowsCut	=g_nCU;
		add(&sync);
        for ( size_t i = 0; i < nRowsCut; ++i ) {
            compute[i] = Stencil2DRowLoop {0,0,this,SHORTWAIT,i};
            add( compute + i );
        }
//		add(&syncToMaster);
	}


	~Stencil2DRowDecomposition(){if ((*timeStep)==0) delete []compute;}
};

#endif
