#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "conf.h"
#include "stencil.h"
#include <math.h>
#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define N_CORES TOTAL_NUM_CU
using namespace darts;

#define GPUMETA 0x4

DEF_CODELET(Stencil2D4ptGpuInitCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuSwapCD,2,LONGWAIT);

DEF_CODELET_ITER(Stencil2D4ptCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(Stencil2D4ptCpuSwapCD,0,SHORTWAIT);

DEF_CODELET(Stencil2D4ptSwapCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(StencilTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t ts;
	double GpuRatio = 0.5;
	uint32_t nCPU;
	Stencil2D4ptGpuKernelCD GpuKernel;
	Stencil2D4ptCpuLoopCD *CpuLoop;
	Stencil2D4ptSwapCD Swap;
	SyncCD	sync;
    Codelet *signalUp;
	uint64_t InitRowCpu;
	uint64_t nRowsGpu;
	uint64_t nRowsCpu;

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
//	,nCPU(N_CORES-1)
//	,GpuKernel(0,1,this,GPUMETA)
//	,CpuLoop(new Stencil2D4ptCpuLoopCD[nCPU])
//	,Swap(nCPU+1,nCPU+1,this,SHORTWAIT)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{
//		std::cout<<"invoke TP!"<<std::endl;	
		if(GpuRatio ==0.0){
			nCPU = N_CORES;
			CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
			Swap = Stencil2D4ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
			InitRowCpu = 0; 
			nRowsCpu = nRows ;	
			for(size_t i=0;i<nCPU; ++i){
				CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
				add(CpuLoop + i);
			}
		}else{

			if (GpuRatio == 1.0){
				nCPU = 0;
			}else{
				nCPU = N_CORES-1;
			}
		
			CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
			InitRowCpu = ceil(nRows*GpuRatio)-1; 
			nRowsGpu = InitRowCpu +2;	
			nRowsCpu = nRows - InitRowCpu ;	
			d_size = sizeof(double)*nRowsGpu*nCols;
			cudaMalloc( (void **) &d_dst, d_size);
			cudaMalloc( (void **) &d_src, d_size);
			
			dim3 dimGrid_hack1((nCols-HALO*2)/GRID_TILE_X,(nRows-HALO*2)/GRID_TILE_Y);
			
			cudaMemcpy(d_src, Initial, d_size, cudaMemcpyHostToDevice);
			cudaMemcpy(d_dst, New, d_size, cudaMemcpyHostToDevice);
			
			GpuKernel = Stencil2D4ptGpuKernelCD{0,1,this,GPUMETA};
			add(&GpuKernel);
		
			for(size_t i=0;i<nCPU; ++i){
				CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
				add(CpuLoop + i);
			}
			Swap = Stencil2D4ptSwapCD{nCPU+1,nCPU+1,this,SHORTWAIT}; 

		}

		
	}
	
	virtual ~StencilTP(){
		cudaFree(d_src); 
		cudaFree(d_dst);
		delete []CpuLoop;
	}
};

#endif
