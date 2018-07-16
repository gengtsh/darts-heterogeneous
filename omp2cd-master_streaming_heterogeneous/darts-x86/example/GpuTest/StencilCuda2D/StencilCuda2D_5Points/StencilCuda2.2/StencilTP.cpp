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


extern "C"
void 
Stencil2D4ptGpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint32_t Id = getID();
	RESET(GpuLoop[Id]);	
	std::cout<<"Invoke GpuLoop ["<<Id<<"]"<<std::endl;	
	uint32_t nGPU = FRAME(nGPU);
	double d_size;	
	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);

	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;

	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

	std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	int64_t bk_nRows = gpu_mem_avail_t/(sizeof(double)*2*nCols)-1;
	std::cout<<"bk_nRows: "<<bk_nRows<<std::endl;
	uint64_t d_nRows = (Id!=(nGPU-1))?(bk_nRows):(nRows-Id*bk_nRows);

	d_size = sizeof(double)*d_nRows*nCols;
	std::cout<<"GPU allocate size: "<< std::setprecision(18)  <<d_size*2<<std::endl;
               
	FRAME(d_size_last) = d_size;
	
	cudaError err1,err2;
	err1 =	cudaMalloc( (void **) &d_dst, d_size);

	if(err1!=cudaSuccess){
		std::cout<<"cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	err2 = cudaMalloc( (void **) &d_src, d_size);
	if(err2!=cudaSuccess){
		std::cout<<"cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}


	dim3 dimGrid_hack1((nCols-HALO*2)/GRID_TILE_X,(d_nRows-HALO*2)/GRID_TILE_Y);

	size_t	pos1 = bk_nRows*Id*nCols;

	cudaMemcpy(d_src, h_src+pos1, d_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dst, h_dst+pos1, d_size, cudaMemcpyHostToDevice);

	gpu_kernel1(dimGrid_hack1,d_dst,d_src,d_nRows,nCols);

	cudaDeviceSynchronize();
	//copy GPU the lowest line to CPU
	cudaMemcpy(h_dst+pos1, d_dst+pos1,d_size , cudaMemcpyDeviceToHost);
	//copy CPU the toppest line to GPU
	cudaFree(d_dst);
	cudaFree(d_src);

	std::cout<<"GPU resource release!"<<std::endl;	
	if ((Id+1)<nGPU){	
		SYNC(GpuLoop[Id+1]);
		std::cout<<"Sync GpuLoop["<<Id+1<<"]"<<std::endl;
	}else{
		std::cout<<"Sync Swap"<<std::endl;
		SYNC(Swap);
	}
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
	std::cout<<"Invoke swap"<<std::endl;	

	double GpuRatio = FRAME(GpuRatio);

	uint32_t	nCPU = FRAME(nCPU);
	uint32_t	nGPU = FRAME(nGPU);
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

	}else if(GpuRatio ==1.0){

		double *d_dst = FRAME(d_dst);
		double *d_src = FRAME(d_src);

		double *tmp;	
		
		if(ts!=0){
			tmp = d_src;
			d_src = d_dst;
			d_dst=tmp;
			for (size_t i =0;i<nGPU;++i){
				SYNC(GpuLoop[i]);
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





