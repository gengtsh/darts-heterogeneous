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

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke Gpu kernel!"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(GpuKernel);
	size_t nRowsGpu   = FRAME(nRowsGpu);

	double *d_dst ;
	double *d_src ;

	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);

	double d_size;	

	d_size = sizeof(double)*nRowsGpu*nCols;
	size_t	gpuPos = FRAME(gpuPos);
	
#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	std::cout<<"gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
#endif
	cudaError err1,err2;
	err1 = cudaMalloc( (void **) &d_dst, d_size);
	err2 = cudaMalloc( (void **) &d_src, d_size);
	
	if(err1!=cudaSuccess){
		std::cout<<"cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	FRAME(d_dst) = d_dst;
	FRAME(d_src) = d_src;

	//std::cout<<"GpuKernel: d_src address:"<<d_src<<std::endl;

	dim3 dimGrid_hack1((nCols-HALO*2)/GRID_TILE_X,(nRowsGpu-HALO*2)/GRID_TILE_Y);
	

	cudaMemcpy(d_src, h_src+gpuPos*nCols, d_size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dst, h_dst+gpuPos*nCols, d_size, cudaMemcpyHostToDevice);
	
	__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),true,false);

	gpu_kernel1(dimGrid_hack1,d_dst,d_src,nRowsGpu,nCols);
	
	SYNC(Swap);	
	if(FRAME(nCPU)==0){
		SYNC(Swap);
	}
	//std::cout<<"GpuKernel: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
	EXIT_TP();
}


extern "C"
void 
Stencil2D4ptGpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint32_t Id = getID();
	RESET(GpuLoop[Id]);	
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuLoop ["<<Id<<"]"<<std::endl;
#endif
	uint32_t nGPU = FRAME(nGPU);
	double d_size;	
	double *d_dst ;
	double *d_src ;
	
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t nRowsGpu = FRAME(nRowsGpu);

//	size_t gpu_mem_total_t = 0;
//	size_t gpu_mem_avail_t = 0;

//	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

//	std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	int64_t chunk=nRowsGpu/nGPU;
	uint64_t nRows_bk = (Id==(nGPU-1))?(nRowsGpu-Id*chunk):(chunk+2);
	d_size = sizeof(double)*nRows_bk*nCols;

	std::cout<<"GpuLoop["<<Id<<"], nRows_bk = "<<nRows_bk<<",allocate size: "<< std::setprecision(18)  <<d_size*2<<std::endl;
	
	cudaError err1,err2,err3,err4,err5,err6;
	err1 = cudaMalloc( (void **) &d_dst, d_size);
	err2 = cudaMalloc( (void **) &d_src, d_size);
	
	if(err1!=cudaSuccess){
		std::cout<<"cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}


	FRAME(d_dst) = d_dst;
	FRAME(d_src) = d_src;
	
	//dim3 dimGrid_hack1((nCols-HALO*2)/GRID_TILE_X,(nRows_bk-HALO*2)/GRID_TILE_Y);
		
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRows_bk/GRID_TILE_Y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuLoop["<<Id<<"], gridDimx="<<gridDimx<<",gridDimy="<<gridDimy<<std::endl;
#endif

	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = (gpuPos+chunk*Id)*nCols;

	err3 = cudaMemcpy(d_dst, h_dst+pos1, d_size, cudaMemcpyHostToDevice);
	err4 = cudaMemcpy(d_src, h_src+pos1, d_size, cudaMemcpyHostToDevice);

	if(err3!=cudaSuccess){
		std::cout<<"cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}
	if(err4!=cudaSuccess){
		std::cout<<"cuda memcpyHostToDevice d_src: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
	gpu_kernel4(dimGrid,dimBlock,d_dst,d_src,nRows_bk,nCols);
	
	err5 = cudaDeviceSynchronize();
	if(err5!=cudaSuccess){
		std::cout<<"cuda deviceSynchronize: "<<cudaGetErrorString(err5)<<std::endl;
		exit(-1);
	}
	//copy GPU  to CPU
	
	err6=cudaMemcpy(h_dst+pos1, d_dst,d_size , cudaMemcpyDeviceToHost);

	if(err6!=cudaSuccess){
		std::cout<<"cuda memcpyDeviceToHost: "<<cudaGetErrorString(err6)<<std::endl;
		exit(-1);
	}
	


#ifdef CUDA_DARTS_DEBUG
	std::cout<<"pos1:"<<pos1<<std::endl;
	std::cout<<"d_size:"<<d_size<<std::endl;
	std::cout<<"nRows_bk:"<<nRows_bk<<std::endl;
	std::cout<<"dst:"<<std::endl;
	std::cout<<std::setprecision(3)<<std::endl;
	int tr = (nRows_bk<10)?nRows_bk:10;
	int tc = (nCols<10)?nCols:10;
	for(size_t i=0;i<tr;++i){
		for (size_t j=0;j<tc;++j){
			std::cout<<h_dst[i*nCols+j]<<",";
		}
		std::cout<<"\n";
	}

	std::cout<<"src:"<<std::endl;
	for(size_t i=0;i<tr;++i){
		for (size_t j=0;j<tc;++j){
			std::cout<<h_src[i*nCols+j]<<",";
		}
		std::cout<<"\n";
	}
#endif

	cudaFree(d_dst);
	cudaFree(d_src);


//	std::cout<<"GPU resource release!"<<std::endl;	
	if ((Id+1)<nGPU){	
		SYNC(GpuLoop[Id+1]);
//		std::cout<<"Sync GpuLoop["<<Id+1<<"]"<<std::endl;
	}else{
//		std::cout<<"Sync Swap"<<std::endl;
		SYNC(Swap);
	}
	EXIT_TP();
}


extern "C"
void 
Stencil2D4ptCpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint64_t Id = getID();
	RESET(CpuLoop[Id]);	
	
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke CPU["<<Id<<"]\n"<<std::endl;	
#endif
	double	*Initial  = FRAME(Initial);
	double	*New      = FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t	nRowsCpu = FRAME(nRowsCpu);	
	uint64_t	nCPU = FRAME(nCPU);
	uint64_t	chunk = (nRowsCpu )/nCPU;
	uint64_t	cpuPos = FRAME(cpuPos);
	size_t		pos1 = (cpuPos + chunk*Id)*nCols;
	uint64_t	nRows_bk = (Id==(nCPU-1))? (nRowsCpu-chunk*Id-1):(chunk +1) ;
	double *src = Initial+pos1;
	double *dst = New+pos1;

	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	uint64_t nRowsLeft = FRAME(nRowsLeft);

	if(nRowsLeft!=0){
//		std::cout<<"CPU check GPU!"<<std::endl;
		if(FRAME(GpuRatio)!=0.0 ){
			if(cudaSuccess == cudaStreamQuery(NULL)){
				if((__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),false,true))) {
			
//					std::cout<<"cpuId["<<Id<<"],gpu kernel finish and invode a new gpu kernel"<<std::endl;
			
					uint64_t nRowsGpu = FRAME(nRowsGpu);
					double d_size = sizeof(double)*nRowsGpu*nCols;
					size_t	gpuPos = FRAME(gpuPos);
					
					cudaMemcpy(dst+gpuPos*nCols,d_dst, d_size, cudaMemcpyDeviceToHost);
					cudaFree(d_dst);
					cudaFree(d_src);
					
					size_t gpu_mem_total_t = 0;
					size_t gpu_mem_avail_t = 0;
					cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
//					std::cout<<"gpu memory total: "<<gpu_mem_total_t<<std::endl;
//					std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
					double	req_size = sizeof(double)* nRowsLeft*nCols*2;
					if(req_size < gpu_mem_avail_t){
						FRAME(nRowsGpu) = nRowsLeft;
						FRAME(nRowsLeft)=0;
					}else if(req_size < 2*gpu_mem_avail_t){
						FRAME(nRowsGpu) = nRowsLeft/2;
						FRAME(nRowsLeft) = nRowsLeft - FRAME(nRowsGpu); 
					}else{
						FRAME(nRowsGpu) = floor(gpu_mem_avail_t/(nCols*2*sizeof(double))); 
						FRAME(nRowsLeft) = nRowsLeft-FRAME(nRowsGpu);	
					}
					FRAME(gpuPos)	= nRows-nRowsLeft; 
					INCR(Swap);
					SYNC(GpuKernel);
					
//					std::cout<<"nCPU ,nRowsGpu:"<<FRAME(nRowsGpu)<<std::endl;
//					std::cout<<"nCPU ,nRowsLeft:"<<FRAME(nRowsLeft)<<std::endl;
				}
			}
		}

	}

	computeInner_stencil2d_v2(dst,src,nRows_bk,nCols);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"pos1:"<<pos1<<std::endl;
	std::cout<<"chunk:"<<chunk<<std::endl;
	std::cout<<"nRows_bk:"<<nRows_bk<<std::endl;
	std::cout<<"dst:"<<std::endl;
	std::cout<<std::setprecision(3)<<std::endl;
	int tr = (nRows_bk<10)?nRows_bk:10;
	int tc = (nCols<10)?nCols:10;
	for(size_t i=0;i<tr;++i){
		for (size_t j=0;j<tc;++j){
			std::cout<<dst[i*nCols+j]<<",";
		}
		std::cout<<"\n";
	}
#endif

	if(FRAME(GpuRatio)==0.0){
		SYNC(Swap);
	}else{
		SYNC(CpuSync);
	}

//	std::cout<<"cpu["<<Id<<"]: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
	EXIT_TP();
}
void Stencil2D4ptCpuSyncCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync invoke!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
	uint64_t nRowsLeft = FRAME(nRowsLeft);
	//uint64_t nRows = FRAME(nRows);
	RESET(CpuSync);
	
	if(nRowsLeft==0.0){
		SYNC(Swap);
//		std::cout<<"CpuSync, no RowsLeft"<<std::endl;
//		std::cout<<"CpuSync, Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;

//		size_t gpu_mem_total_t = 0;
//		size_t gpu_mem_avail_t = 0;
//		cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
//		std::cout<<"gpu memory total: "<<gpu_mem_total_t<<std::endl;
//		std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	
	}else{
		
		size_t	gpuPos = FRAME(gpuPos);
		size_t	cpuPos = FRAME(cpuPos);
	
		uint64_t nRowsGpu = FRAME(nRowsGpu);
		uint64_t nRowsCpu = FRAME(nRowsCpu);
		uint64_t nCols = FRAME(nCols);

		double *dst = FRAME(New);
		double *d_dst = FRAME(d_dst);
		double *d_src = FRAME(d_src);
		double cmCpu = FRAME(cmCpu);
		double cmGpu = FRAME(cmGpu);
		
		if((cudaSuccess == cudaStreamQuery(NULL) )&&(__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),false,true))) {
//			std::cout<<"cpuSync, gpu kernel finish and invode a new gpu kernel"<<std::endl;

			size_t gpu_mem_total_t = 0;
			size_t gpu_mem_avail_t = 0;
			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
//			std::cout<<"gpu memory total: "<<gpu_mem_total_t<<std::endl;
//			std::cout<<"gpu memory available: "<<gpu_mem_avail_t<<std::endl;
		
			double d_size = sizeof(double)*nRowsGpu*nCols;
			cudaMemcpy(dst+gpuPos*nCols,d_dst, d_size, cudaMemcpyDeviceToHost);
			cudaFree(d_dst);
			cudaFree(d_src);
			
			double	req_size = sizeof(double)* nRowsLeft*nCols*2;
			if(req_size < gpu_mem_avail_t){
				
				FRAME(gpuPos) = nRowsCpu + nRowsGpu;
				FRAME(nRowsGpu) = nRowsLeft;
				FRAME(nRowsLeft)=0;
				INCR(Swap);
				SYNC(GpuKernel);
				
			}else if(req_size < 2*gpu_mem_avail_t){
				
				FRAME(gpuPos) = nRowsCpu + nRowsGpu;
				FRAME(nRowsGpu) = nRowsLeft/2;
				FRAME(cpuPos) = FRAME(gpuPos) + FRAME(nRowsGpu);
				FRAME(nRowsCpu) = nRowsLeft - FRAME(nRowsGpu); 
				FRAME(nRowsLeft) = 0;
				for(size_t i =0; i<FRAME(nCPU);++i){
					SYNC(CpuLoop[i]);
				}
				INCR(Swap);
				SYNC(GpuKernel);
			}else{
				FRAME(gpuPos) = nRowsCpu + nRowsGpu;
				nRowsGpu = floor(gpu_mem_avail_t/(nCols*2)); 
				FRAME(nRowsGpu) = nRowsGpu;
				FRAME(cpuPos) = FRAME(gpuPos) + FRAME(nRowsGpu);
				nRowsCpu = nRowsGpu*(cmCpu/cmGpu);
				nRowsLeft = nRowsLeft - nRowsGpu - nRowsCpu;
				if (nRowsLeft< (1/cmGpu)*nRowsGpu){
						nRowsCpu = nRowsLeft + nRowsCpu;
						nRowsLeft = 0;
						
				}
								
				FRAME(nRowsCpu) = nRowsCpu;
				FRAME(nRowsLeft) = nRowsLeft;

				for(size_t i =0; i<FRAME(nCPU);++i){
					SYNC(CpuLoop[i]);
				}
				INCR(Swap);
				SYNC(GpuKernel);
			}
		}else{
//			std::cout<<"CpuSync, invoke nCPU"<<std::endl;
			FRAME(cpuPos) = gpuPos + nRowsGpu;
			if(nRowsLeft < nRowsCpu*(cmGpu/cmCpu)){
				nRowsCpu = nRowsLeft;
				nRowsLeft = 0;
			}else{
				nRowsLeft = nRowsLeft - nRowsCpu;
			}
			FRAME(nRowsCpu) = nRowsCpu;
			FRAME(nRowsLeft) = nRowsLeft;

			for(size_t i =0; i<FRAME(nCPU);++i){
				SYNC(CpuLoop[i]);
			}
//			std::cout<<"CpuSync, nRowsLeft:"<<FRAME(nRowsLeft)<<std::endl;
		}


	//	std::cout<<"CpuSync, Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
	}
	EXIT_TP();

}


void 
Stencil2D4ptSwapCD::fire(void) 
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke swap"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(Swap);

	double GpuRatio = FRAME(GpuRatio);

	uint32_t	nCPU = FRAME(nCPU);
	uint32_t	nGPU = FRAME(nGPU);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"ts: "<<FRAME(ts)<<std::endl;
#endif	
	size_t ts = --FRAME(ts);
	if(GpuRatio == 0.0){

		double *dst = FRAME(New);
		double *src = FRAME(Initial);
		
		if(ts!=0){
			SWAP_PTR(&dst,&src);
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}
		}else{
			SYNC(sync);
		}

	}else if(GpuRatio ==1.0){

		double *dst = FRAME(New);
		double *src = FRAME(Initial);
		
		if(ts!=0){
			
			SWAP_PTR(&dst,&src);
			for (size_t i =0;i<nGPU;++i){
				SYNC(GpuLoop[i]);
			}
		}else{
			SYNC(sync);
		}


	}else{

		double *d_dst = FRAME(d_dst);
		double *d_src = FRAME(d_src);
	//	std::cout<<"swap: d_src address:"<<d_src<<std::endl;
		double *dst = FRAME(New);
		double *src = FRAME(Initial);

		uint64_t	gpuPos = FRAME(gpuPos);
		uint64_t	cpuPos = FRAME(cpuPos);
		
		uint64_t nRows   = FRAME(nRows);
		uint64_t nCols   = FRAME(nCols);
		uint64_t nRowsGpu = FRAME(nRowsGpu);
		uint64_t nRowsCpu = FRAME(nRowsCpu);
		
//		uint64_t	pos1 = (cpuPos+nRowsCpu)*nCols;
//		uint64_t	pos2 = (cpuPos+1)*nCols;
//		double		line_size = sizeof(double*)*nCols;
//		
//		cudaDeviceSynchronize();
//		//copy GPU the lowest line to CPU
//		cudaMemcpy(dst+pos1, d_dst+pos1,line_size , cudaMemcpyDeviceToHost);
//		//copy CPU the toppest line to GPU
//		cudaMemcpy(d_dst+pos2, dst, line_size, cudaMemcpyHostToDevice);
	
		
		double d_size = sizeof(double)*nRowsGpu*nCols;

		cudaDeviceSynchronize();
		

		cudaError err1,err2,err3;
		cudaMemcpy(dst+gpuPos*nCols, d_dst,d_size , cudaMemcpyDeviceToHost);
	
		err1 = cudaFree(d_dst);
		err2 = cudaFree(d_src);

		if(err1!=cudaSuccess){
			std::cout<<"swap: cudaFree d_dst: "<<cudaGetErrorString(err1)<<std::endl;
			exit(-1);
		}

		if(err2!=cudaSuccess){
			std::cout<<"swap: cudaFree d_src: "<<cudaGetErrorString(err2)<<std::endl;
			exit(-1);
		}
		size_t gpu_mem_total_t = 0;
		size_t gpu_mem_avail_t = 0;

		err3 = cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

		if(err3!=cudaSuccess){
			std::cout<<"swap: cuda get mem info: "<<cudaGetErrorString(err3)<<std::endl;
			exit(-1);
		}
//		std::cout<<"swap: gpu memory total: "<<gpu_mem_total_t<<std::endl;
//		std::cout<<"swap: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
		
		double *tmp;	
		if(ts!=0){
		//	tmp = d_src;
		//	d_src = d_dst;
		//	d_dst=tmp;
			SWAP_PTR(&dst,&src);
			SYNC(GpuKernel);
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}

		}else{
			SYNC(sync);
		}
	}
	EXIT_TP();
}
void
SyncCD::fire(void)
{
  
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Sync!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





