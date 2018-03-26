#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "StencilTP.h"
#include "stencil.h"
#include "StencilCPUKernel.h"

#include <cassert>
//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
#include <iostream>

extern "C"
void 
Stencil2D4ptGpuKernelCD::fire(void) 
{
	LOAD_FRAME(StencilTP);
	RESET(GpuKernel);
	size_t nRowsGpu   = FRAME(nRowsGpu);   
	size_t nCols   = FRAME(nCols);   
//	std::cout<<"Invoke kernel"<<std::endl;	
	dim3 dimGrid_hack1 = FRAME(dimGrid_hack1); 
	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	
	gpu_kernel1(dimGrid_hack1,d_dst,d_src,nRowsGpu,nCols);
	SYNC(Swap);	
	EXIT_TP();
}

void 
Stencil2D4ptCpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint64_t Id = getID();
	RESET(CpuLoop[Id]);	
//	std::cout<<"Invoke CPU["<<Id<<"]\n"<<std::endl;	

	double	*Initial  = FRAME(Initial);
	double	*New      = FRAME(New);
//	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t	nRowsCpu = FRAME(nRowsCpu);	
	uint64_t	nCPU = FRAME(nCPU);
	uint64_t	chunk = (nRowsCpu -2)/nCPU;
	uint64_t	InitRowCpu = FRAME(InitRowCpu);
	size_t		pos1 = InitRowCpu + chunk*Id*nCols+nCols;
	uint64_t rpos2 = (Id==(nCPU-1))? (nRowsCpu-2-chunk*nCPU):chunk ;
	double *dst = Initial+pos1;
	double *src = New+pos1;
	
	computeInner_stencil2d_v2(dst,src,rpos2,nCols);

	SYNC(Swap);

	EXIT_TP();
}

void 
Stencil2D4ptSwapCD::fire(void) 
{
	LOAD_FRAME(StencilTP);
	RESET(Swap);
//	std::cout<<"Invoke swap"<<std::endl;	

	double GpuRatio = FRAME(GpuRatio);

	uint64_t	nCPU = FRAME(nCPU);

	size_t ts = --FRAME(ts);
	if(GpuRatio == 0.0){

		double *dst = FRAME(New);
		double *src = FRAME(Initial);
		
		SWAP_PTR(&dst,&src);
		if(ts!=0){
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}
		}else{
			SYNC(sync);
		}

	}else{

		double *d_dst = FRAME(d_dst);
		double *d_src = FRAME(d_src);
		double *dst = FRAME(New);
		double *src = FRAME(Initial);
		double d_size = FRAME(d_size);
		
	//	uint64_t nRows   = FRAME(nRows);
		uint64_t nCols   = FRAME(nCols);
	//	uint64_t nRowsGpu = FRAME(nRowsGpu);
	//	uint64_t nRwosCpu = FRAME(nRowsCpu);
		
		uint64_t	InitRowCpu = FRAME(InitRowCpu);
		uint64_t	pos1 = InitRowCpu*nCols;
		uint64_t	pos2 = (InitRowCpu+1)*nCols;
		double		line_size = sizeof(double*)*nCols;
		
		cudaDeviceSynchronize();
		//copy GPU the lowest line to CPU
		cudaMemcpy(dst+pos1, d_dst+pos1,line_size , cudaMemcpyDeviceToHost);
		//copy CPU the toppest line to GPU
		cudaMemcpy(d_dst+pos2, dst+pos2, line_size, cudaMemcpyHostToDevice);
		
		double *tmp;	
		if(ts!=0){
			tmp = d_src;
			d_src = d_dst;
			d_dst=tmp;
			SYNC(GpuKernel);
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}
		}else{
			cudaMemcpy(dst, d_dst, d_size, cudaMemcpyDeviceToHost);
			cudaMemcpy(src, d_src, d_size, cudaMemcpyDeviceToHost);
			SYNC(sync);
		}
	}
	EXIT_TP();
}
void
SyncCD::fire(void)
{
  //  std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(StencilTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





