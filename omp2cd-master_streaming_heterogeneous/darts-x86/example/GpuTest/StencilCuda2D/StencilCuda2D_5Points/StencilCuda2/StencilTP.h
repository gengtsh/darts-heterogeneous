#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "conf.h"
#include "stencil.h"

#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define GPUMETA 0x4

DEF_CODELET(Stencil2D4ptGpuInitCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuSwapCD,2,LONGWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(StencilTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t ts;
	Stencil2D4ptGpuKernelCD GpuKernel;
	Stencil2D4ptGpuSwapCD GpuSwap;
	SyncCD	sync;
    Codelet *signalUp;
	
	double d_size;
	double *d_dst;
	double *d_src;
	dim3 dimGrid_hack1;
	StencilTP(double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t ts,Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,ts(ts)
	,GpuKernel(0,1,this,GPUMETA)
	,GpuSwap(1,1,this,GPUMETA)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{
		d_size = sizeof(double)*nRows*nCols;
		cudaMalloc( (void **) &d_dst, d_size);
		cudaMalloc( (void **) &d_src, d_size);
		
		dim3 dimGrid_hack1((nCols-HALO*2)/GRID_TILE_X,(nRows-HALO*2)/GRID_TILE_Y);
		
		cudaMemcpy(d_src, Initial, d_size, cudaMemcpyHostToDevice);
		cudaMemcpy(d_dst, New, d_size, cudaMemcpyHostToDevice);
		add(&GpuKernel);
	}
	
	virtual ~StencilTP(){
	
		cudaFree(d_src); 
		cudaFree(d_dst);
	}
};

#endif
