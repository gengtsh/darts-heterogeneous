#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <cstdint>
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES
using namespace darts;


DEF_CODELET_ITER(Stencil2DRowLoop,0,SHORTWAIT);
//DEF_CODELET(Stencil2DRowSyncSwap,2,LONGWAIT);
DEF_CODELET_ITER(Stencil2DRowSyncSwap,0,SHORTWAIT);
DEF_CODELET(Stencil2DRowSync,2,LONGWAIT);

DEF_TP(Stencil2DRowDecomposition)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t timeStep;
    
	Stencil2DRowLoop *compute;
    Stencil2DRowSyncSwap *syncSwap;
	Stencil2DRowSync	sync;
    Codelet *signalUp;
	uint32_t nRowsCut;
	uint64_t *tSteps;	
	uint64_t chunk;

//    typedef double (*Array2D)[nCols];
//	Array2D *dstPtr;
//	Array2D *srcPtr;
	double **dstPtr;
	double **srcPtr;

	Stencil2DRowDecomposition(double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t ts, Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	,compute(new Stencil2DRowLoop[N_CORES])
	,syncSwap(new Stencil2DRowSyncSwap[N_CORES])
	,sync(N_CORES,N_CORES,this,LONGWAIT)
	,signalUp(up)
    ,nRowsCut(N_CORES)
	,tSteps(new uint64_t[N_CORES])
	{
//		nRowsCut	= computeRowDecomposition(inim,inin);
//		compute		= new Stencil2DRowLoop[nRowsCut];
//		syncSwap	= Stencil2DRowSyncSwap{nRowsCut,nRowsCut,this,LONGWAIT};
		//std::cout<<"N_CORES"<<N_CORES<<std::endl;
		//add(&sync);
		setRef(9999);
		dstPtr  = new double* [N_CORES];
		srcPtr  = new double* [N_CORES];

		chunk= (nRows-2)/nRowsCut;
		for(size_t i=0;i<nRowsCut;++i){
			size_t pos1=1+chunk*i*nCols+nCols;
			dstPtr[i] = Initial+pos1;
			srcPtr[i] = New+pos1;
		}
		for ( size_t i = 0; i < nRowsCut; ++i ) {
			tSteps[i]=ts;
            compute[i] = Stencil2DRowLoop {0,3,this,SHORTWAIT,i};
			syncSwap[i]= Stencil2DRowSyncSwap{3,3,this,SHORTWAIT,i};
            add( compute + i );
//			add(syncSwap + i );
        }


	}
	~Stencil2DRowDecomposition(){
		delete []compute;
		delete []tSteps;
		delete []syncSwap;
		delete []dstPtr;
		delete []srcPtr;
	}
};

#endif
