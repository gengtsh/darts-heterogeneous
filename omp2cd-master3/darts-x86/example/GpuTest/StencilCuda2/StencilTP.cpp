#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "StencilTP.h"
#include "stencil.h"
#include <cassert>
//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
#include <iostream>

//version1 is just for testing whether swap cd will wait for kernel finishing or not. and version1i has proved that swapcd will be invoke imediately and not wait for kernel finishing, to guarantee the test correct, we need to add cudaDeviceSynchronize() function in swapcd
#define version 1 

extern "C"
void 
Stencil2D4ptGpuKernelCD::fire(void) 
{
	LOAD_FRAME(StencilTP);
	RESET(GpuKernel);
	size_t n_rows   = FRAME(nRows);   // matrix M row
	size_t n_cols   = FRAME(nCols);   // Matrix N column
	std::cout<<"Invoke kernel"<<std::endl;	
	dim3 dimGrid_hack1 = FRAME(dimGrid_hack1); 
	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
#if version == 1
	// version 1	
	gpu_kernel1(dimGrid_hack1,d_dst,d_src,n_rows,n_cols);
	SYNC(GpuSwap);	
#elif version == 2
	//version 2	
	double *dst = FRAME(New);
	double *src = FRAME(Initial);
	double ts = FRAME(ts);
	double d_size = FRAME(d_size);

	gpu_kernel2(dimGrid_hack1,dst, src, d_size, ts, d_dst, d_src, n_rows, n_cols);
	SYNC(GpuSwap);	
#endif
	EXIT_TP();
}


void 
Stencil2D4ptGpuSwapCD::fire(void) 
{
	LOAD_FRAME(StencilTP);
	RESET(GpuSwap);
	std::cout<<"Invoke swap"<<std::endl;	
#if version == 1  	
	//version 1
	size_t ts = --FRAME(ts);
	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	double *tmp;	
	if(ts!=0){
		tmp = d_src;
		d_src = d_dst;
		d_dst=tmp;
		cudaDeviceSynchronize();
		SYNC(GpuKernel);
	}else{
		double *dst = FRAME(New);
		double *src = FRAME(Initial);
		double d_size = FRAME(d_size);
		cudaMemcpy(dst, d_dst, d_size, cudaMemcpyDeviceToHost);
		cudaMemcpy(src, d_src, d_size, cudaMemcpyDeviceToHost);
		SYNC(sync);
	}
#elif version == 2
	//version 2
	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	double *dst = FRAME(New);
	double *src = FRAME(Initial);
	double d_size = FRAME(d_size);

	cudaDeviceSynchronize();
	cudaMemcpy(dst, d_dst, d_size, cudaMemcpyDeviceToHost);
	cudaMemcpy(src, d_src, d_size, cudaMemcpyDeviceToHost);
	SYNC(sync);
#endif
	EXIT_TP();
}
void
SyncCD::fire(void)
{
    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(StencilTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





