#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <stdint.h>
#include <stdlib.h>
#include <darts.h>
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"
#include "DARTS.h"
#include "Stencil2DPartition.h"
using namespace darts;


DEF_CODELET_ITER(Stencil2DRowLoopCopyUp,2,SHORTWAIT);
DEF_CODELET_ITER(Stencil2DRowLoopCopyDown,2,SHORTWAIT);
DEF_CODELET_ITER(Stencil2DRowLoop,4,SHORTWAIT);
DEF_CODELET(Stencil2DRowSync,2,LONGWAIT);

/*
*in Stencil2DRowDecomposition TP:
*
*Stencil2DRowLoop has 4 dependences:up0,down0,up1,down1
*	up0,down0: the uppest and lowest line of sub-array(inner of subarray)
*	up1,down1: the upper and lower line of sub-array(out side of subarray)
*	4 dependences means: the computation can't be execute until copy(up0,down0,up1,down1) finished. 
*
*checkUP,checkDown,ComputeUp,computeDown are related to these 4 dependences.
*
*
*/

DEF_TP(Stencil2DRowDecomposition)
{
	double *initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	uint64_t timeStep;
	uint32_t nTp;

	Stencil2DPartitionCheckTP *checkTP;
	Stencil2DRowLoopCopyUp *copyUp;//codelet copy shared upper line 
	Stencil2DRowLoopCopyDown *copyDown;//codelet copy shared Down line
	Stencil2DRowLoop *compute;
	Stencil2DRowSync  sync;
	
	Codelet*   signalUp;
	
	Stencil2DRowCheck *check;
	double * shareRows;//every nt has an inner matrix which is used to story original data 
    uint64_t *tSteps;
	uint32_t nRowsCut;
	bool *finalize;
	

	Stencil2DRowDecomposition(double *inimatrix,const uint64_t inim,const uint64_t inin,uint64_t ts,uint32_t tp,Stencil2DPartitionCheckTP *pCheckTP,Codelet *up, Stencil2DRowCheck  *rCheck=0)
	:initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,timeStep(ts)
	,nTp(tp)
	,checkTP(pCheckTP)
	,signalUp(up)
//	,check(rCheck?rCheck:new Stencil2DRowCheck[2] )
	,check(rCheck)
	{
	//	nRowsCut	= computeRowDecomposition(inim,inin);
		nRowsCut	= g_nCU;
		copyUp		= new Stencil2DRowLoopCopyUp[nRowsCut];
		copyDown	= new Stencil2DRowLoopCopyDown[nRowsCut];
		compute		= new Stencil2DRowLoop[nRowsCut];
		shareRows	= new double[nRowsCut*2*inin];
		tSteps		= new uint64_t[nRowsCut];
		sync		= Stencil2DRowSync{nRowsCut*3,nRowsCut*3,this,LONGWAIT};
		finalize	= new bool[nRowsCut];	 
//		add(&sync);

		for ( size_t i = 0; i < nRowsCut; ++i ) {
				tSteps[i]=ts;
				finalize[i]=false;
				copyUp[i] = Stencil2DRowLoopCopyUp {0,2,this,SHORTWAIT,i};
				copyDown[i] = Stencil2DRowLoopCopyDown {0,2,this,SHORTWAIT,i};

           		compute[i] = Stencil2DRowLoop {4,4,this,SHORTWAIT,i};
				add (copyUp +i);
				add (copyDown+i);
		}

		check[0] = Stencil2DRowCheck{1,1,this,SHORTWAIT,0};
		check[1] = Stencil2DRowCheck{1,1,this,SHORTWAIT,1};
	}
	static uint32_t computeRowDecomposition(const uint64_t n_rows,const uint64_t n_cols){

        return ((n_rows-2)*(n_cols-2)/TOTAL_TILE_SZ < 1) ? 1 : g_nCU*N_THREADS;
	}

	virtual	~Stencil2DRowDecomposition(){
		delete []copyUp;
		delete []copyDown;
		delete []compute;
		delete []shareRows;
		delete []tSteps;
		
		//std::cout<<"RowFinish,nTP:"<<nTp<<std::endl;	
	}
};

#endif
