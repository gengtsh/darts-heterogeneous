#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <cstdint>
#include <cstdlib>
#include "Stencil2D_main.h"
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES
using namespace darts;

DEF_CODELET_ITER(Stencil2DRowLoopCopyUp,2,SHORTWAIT);
DEF_CODELET_ITER(Stencil2DRowLoopCopyDown,2,SHORTWAIT);
DEF_CODELET_ITER(Stencil2DRowLoop,4,SHORTWAIT);
DEF_CODELET(Stencil2DRowSync,2,LONGWAIT);

DEF_TP(Stencil2DRowDecomposition)
{
//	uint64_t nRowsCut;
	double *initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	uint64_t timeStep;	
	Stencil2DRowLoopCopyUp *copyUp;//codelet copy shared upper line 
	Stencil2DRowLoopCopyDown *copyDown;//codelet copy shared Down line
	Stencil2DRowLoop *compute;
	Stencil2DRowSync  sync;
	Codelet*   signalUp;
	double * shareRows;//every nt has an inner matrix which is used to story original data 
	uint64_t *tSteps;
	uint32_t nRowsCut;
	bool *finalize;

	Stencil2DRowDecomposition(double *inimatrix,const uint64_t inim,const uint64_t inin,uint64_t ts, Codelet *up)
//	:nRowsCut(nt)
	 :initial(inimatrix)
	 ,nRows(inim)
	 ,nCols(inin)
	 ,timeStep(ts)
	 ,signalUp(up)
	{

		//nRowsCut	= computeRowDecomposition(inim,inin);
		setRef(9999);
		nRowsCut	= N_CORES;
		copyUp		= new Stencil2DRowLoopCopyUp[nRowsCut];
		copyDown	= new Stencil2DRowLoopCopyDown[nRowsCut];
		compute		= new Stencil2DRowLoop[nRowsCut];
		shareRows	= new double[nRowsCut*2*inin];
		tSteps		= new uint64_t[nRowsCut];
		sync		= Stencil2DRowSync{nRowsCut,nRowsCut,this,LONGWAIT};
		finalize	= new bool[nRowsCut];	 
		 add(&sync);
        
		for ( size_t i = 0; i < nRowsCut; ++i ) {
				tSteps[i]=ts;
				finalize[i]=false;
				copyUp[i] = Stencil2DRowLoopCopyUp {0,2,this,SHORTWAIT,i};
				copyDown[i] = Stencil2DRowLoopCopyDown {0,2,this,SHORTWAIT,i};

           		compute[i] = Stencil2DRowLoop {4,4,this,SHORTWAIT,i};

				add (copyUp+i);
				add (copyDown+i);
		}
	}
    static uint32_t computeRowDecomposition(const uint64_t n_rows, const uint64_t n_cols) {
//      the total number of pixel in inner matrix 
//      (whole_matrix - left_edge - right_edge - upper_edge - down_edge)
//      N_THREADS are number of threads per cores
        return ((n_rows-2)*(n_cols-2)/TOTAL_TILE_SZ < 1) ? 1 : TOTAL_NUM_CU*N_THREADS;
    }


	~Stencil2DRowDecomposition(){
		delete []copyUp;
		delete []copyDown;
		delete []compute;
		delete []shareRows;
		delete []tSteps;
		delete []finalize;
	}
};

#endif
