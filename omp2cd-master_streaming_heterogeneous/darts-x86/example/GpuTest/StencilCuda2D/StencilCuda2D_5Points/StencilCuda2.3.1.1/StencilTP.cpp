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
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernel: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernel: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernel: require memory size:"<<d_size*2<<std::endl;
	std::cout<<"GpuKernel: Gpupos:"<<gpuPos<<std::endl;
	std::cout<<"GpuKernel: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernel: nRowsGpu:"<<nRowsGpu<<std::endl;

#endif
	cudaError err1,err2,err3,err4;
	err1 = cudaMalloc( (void **) &d_dst, d_size);
	err2 = cudaMalloc( (void **) &d_src, d_size);
	
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernel: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuKernel: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	FRAME(d_dst) = d_dst;
	FRAME(d_src) = d_src;

	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRowsGpu/GRID_TILE_Y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	err3 = cudaMemcpy(d_src, h_src+gpuPos*nCols, d_size, cudaMemcpyHostToDevice);
	err4 = cudaMemcpy(d_dst, h_dst+gpuPos*nCols, d_size, cudaMemcpyHostToDevice);

	if(err3 != cudaSuccess){
		std::cout<<"GpuKernel: cuda memcpy HostToDevice d_src: "<<cudaGetErrorString(err3)<<std::endl;
	}

	if(err4 != cudaSuccess){
		std::cout<<"GpuKernel: cuda memcpy dst HostToDevice d_dst: "<<cudaGetErrorString(err4)<<std::endl;
	}

	__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),true,false);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"FRAME(CpuIvGpu): "<<FRAME(CpuIvGpu)<<std::endl;
#endif
	gpu_kernel4(dimGrid,dimBlock,d_dst,d_src,nRowsGpu,nCols);

	SYNC(Swap);	

	if(FRAME(nCPU)==0){
		SYNC(Swap);
	}
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernel: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif

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
	double *d_sharedCols ;
	double *d_sharedRows ;
	
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t nRowsGpu = FRAME(nRowsGpu);

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

	gpu_mem_valid_t = gpu_mem_avail_t -XMB;
	
	std::cout<<"GpuLoop: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
#endif
	int64_t chunk=nRowsGpu/nGPU;
	uint64_t nRows_bk = (Id==(nGPU-1))?(nRowsGpu-Id*chunk):(chunk+2);
	d_size = sizeof(double)*nRows_bk*nCols;
	
	int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRows_bk/tile_y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = (gpuPos+chunk*Id)*nCols;

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuLoop["<<Id<<"]: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: nRows_bk:"<<nRows_bk<<std::endl;

	std::cout<<"GpuLoop["<<Id<<"]: blockDimx:"<<blockDimx<<std::endl;
	
	std::cout<<"GpuLoop["<<Id<<"]: gridDimx="<<gridDimx<<",gridDimy="<<gridDimy<<std::endl;
#endif
	
	int64_t d_size_sharedCols = sizeof(double)*nRows_bk*gridDimx*2;
	int64_t d_size_sharedRows = sizeof(double)*nCols*gridDimy*2;

	cudaError err1,err2,err3,err4,err5,err6,err7,err8,err9,err10, err11;
	err1 = cudaMalloc( (void **) &d_dst, d_size);
//	err2 = cudaMalloc( (void **) &d_src, d_size);

	err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
	err9 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuLoop: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuLoop: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	if(err9!=cudaSuccess){
		std::cout<<"GpuLoop: cuda malloc9: "<<cudaGetErrorString(err9)<<std::endl;
		exit(-1);
	}

#endif

	FRAME(d_dst) = d_dst;
	FRAME(d_src) = d_src;
		
	err3 = cudaMemcpy(d_dst, h_dst+pos1, d_size, cudaMemcpyHostToDevice);
//	err4 = cudaMemcpy(d_src, h_src+pos1, d_size, cudaMemcpyHostToDevice);


#ifdef CUDA_ERROR_CHECKING
	if(err3!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

#endif

	//gpu_kernel4(dimGrid,dimBlock,d_dst,d_src,nRows_bk,nCols);

	int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
	int blockDimy_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	int gridDimy_rows = std::ceil(1.0*nRows_bk/tile_y);


	int blockDimx_cols = 1 ;
	int blockDimy_cols = (nRows_bk>NUM_THREADS)?NUM_THREADS:nRows_bk;
	int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows_bk/blockDimy_cols);

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuLoop["<<Id<<"]: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: grimDimy_rows:"<<gridDimy_rows<<std::endl;

	std::cout<<"GpuLoop["<<Id<<"]: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	std::cout<<"GpuLoop["<<Id<<"]: grimDimy_cols:"<<gridDimy_cols<<std::endl;

#endif
	
	dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);

	dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
	
	gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,nRows_bk, nCols);

	gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,nRows_bk, nCols);

	int nCmpRep = FRAME(nCmpRep);
	for (size_t j=0;j<nCmpRep;++j){
		gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,nRows_bk,nCols);
	}
	
	err5 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err5!=cudaSuccess){
		std::cout<<"GpuLoop: cuda deviceSynchronize: "<<cudaGetErrorString(err5)<<std::endl;
		exit(-1);
	}
#endif
	//copy from GPU  to CPU
	err6=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
	if(err6!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err6)<<std::endl;
		exit(-1);
	}
#endif

#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"dst:"<<std::endl;
//	std::cout<<std::setprecision(3)<<std::endl;
//	int tr = (nRows_bk<10)?nRows_bk:10;
//	int tc = (nCols<10)?nCols:10;
//	for(size_t i=0;i<tr;++i){
//		for (size_t j=0;j<tc;++j){
//			std::cout<<h_dst[i*nCols+j]<<",";
//		}
//		std::cout<<"\n";
//	}
//
//	std::cout<<"src:"<<std::endl;
//	for(size_t i=0;i<tr;++i){
//		for (size_t j=0;j<tc;++j){
//			std::cout<<h_src[i*nCols+j]<<",";
//		}
//		std::cout<<"\n";
//	}
#endif

	err7 = cudaFree(d_dst);
//	err8 = cudaFree(d_src);
	err8 = cudaFree(d_sharedCols);
	err11 = cudaFree(d_sharedRows);
	
	if(err7!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_dst: "<<cudaGetErrorString(err7)<<std::endl;
		exit(-1);
	}

	if(err8!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_src: "<<cudaGetErrorString(err8)<<std::endl;
		exit(-1);
	}

	if(err11!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_src: "<<cudaGetErrorString(err11)<<std::endl;
		exit(-1);
	}

	if ((Id+1)<nGPU){	
		SYNC(GpuLoop[Id+1]);
	}else{
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
	std::cout<<"Invoke CPU["<<Id<<"]"<<std::endl;	
#endif
	double	*h_src  = FRAME(Initial);
	double	*h_dst      = FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t	nRowsCpu = FRAME(nRowsCpu);	
	uint64_t	nCPU = FRAME(nCPU);
	uint64_t	chunk = (nRowsCpu )/nCPU;
	uint64_t	cpuPos = FRAME(cpuPos);
	size_t		pos1 = (cpuPos + chunk*Id)*nCols;
	uint64_t	nRows_bk = (Id==(nCPU-1))? (nRowsCpu-chunk*Id-1):(chunk +1) ;
	double *src = h_src+pos1;
	double *dst = h_dst+pos1;

	double *d_dst = FRAME(d_dst);
	double *d_src = FRAME(d_src);
	uint64_t nRowsLeft = FRAME(nRowsLeft);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuLoop["<<Id<<"]: cpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nRowsCpu:"<<nRowsCpu<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: chunk:"<<chunk<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nRows_bk:"<<nRows_bk<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nRowsLeft:"<<nRowsLeft<<std::endl;
#endif


	if(nRowsLeft!=0){

#ifdef CUDA_DARTS_DEBUG
		std::cout<<"CPULoop["<<Id<<"]: check GPU!"<<std::endl;
#endif
		if(FRAME(GpuRatio)!=0.0 ){
			if(cudaSuccess == cudaStreamQuery(NULL)){
				if((__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),false,true))) {
			
#ifdef CUDA_DARTS_DEBUG
					std::cout<<"CPULoop["<<Id<<"]: gpu kernel finish and invode a new gpu kernel"<<std::endl;
#endif			
					uint64_t nRowsGpu = FRAME(nRowsGpu);
					double d_size = sizeof(double)*(nRowsGpu-2)*nCols;
					size_t	gpuPos = FRAME(gpuPos);
					
					cudaError err1,err2,err3;
					err1 = cudaMemcpy(h_dst+(gpuPos+1)*nCols,d_dst+nCols, d_size, cudaMemcpyDeviceToHost);
					err2 = cudaFree(d_dst);
					err3 = cudaFree(d_src);
				

					if(err1!=cudaSuccess){
						std::cout<<"CpuLoop: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err1)<<std::endl;
						exit(-1);
					}

					if(err2!=cudaSuccess){
						std::cout<<"CpuLoop: cuda mem free d_dst: "<<cudaGetErrorString(err2)<<std::endl;
						exit(-1);
					}

					if(err3!=cudaSuccess){
						std::cout<<"CpuLoop: cuda mem free d_src: "<<cudaGetErrorString(err2)<<std::endl;
						exit(-1);
					}

#ifdef CUDA_DARTS_DEBUG
			//		std::cout<<"CPULoop["<<Id<<"]: last gpu kernel dst:"<<std::endl;
			//		std::cout<<std::setprecision(3)<<std::endl;
			//		int tr = 10;
			//		int tc = 10;
			//		double *lastdst=h_dst+gpuPos*nCols;
			//		for(size_t i=0;i<tr;++i){
			//			for (size_t j=0;j<tc;++j){
			//				std::cout<<lastdst[i*FRAME(nCols)+j]<<",";
			//			}
			//			std::cout<<"\n";
			//		}

#endif

					size_t gpu_mem_total_t = 0;
					size_t gpu_mem_avail_t = 0;
				
					size_t gpu_mem_valid_t = 0;
					cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
					gpu_mem_valid_t = gpu_mem_avail_t -XMB;

					double	req_size = sizeof(double)* nRowsLeft*nCols*2;
#ifdef CUDA_DARTS_DEBUG
					std::cout<<"CpuLoop["<<Id<<"]: gpu memory total: "<<gpu_mem_total_t<<std::endl;
					std::cout<<"CpuLoop["<<Id<<"]: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
					std::cout<<"CpuLoop["<<Id<<"]: require gpu memory:"<<req_size<<std::endl;
#endif
					if(req_size < gpu_mem_valid_t){
						FRAME(nRowsGpu) = nRowsLeft;
						FRAME(nRowsLeft)=0;
						INCR(Swap);
					}else if(req_size < 2*gpu_mem_valid_t){
						FRAME(nRowsGpu) = nRowsLeft/2;
						FRAME(nRowsLeft) = nRowsLeft - FRAME(nRowsGpu)+2; 
						INCR(Swap);
					}else{
						FRAME(nRowsGpu) = floor(gpu_mem_valid_t/(nCols*2*sizeof(double))); 
						FRAME(nRowsLeft) = nRowsLeft-FRAME(nRowsGpu)+2;	
						INCR(Swap);
					}
					FRAME(gpuPos)	= nRows-nRowsLeft; 
					SYNC(GpuKernel);
					
#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuLoop["<<Id<<"]: invoke new GpuKernel: gpuPos: "<<FRAME(gpuPos)<<std::endl;
			std::cout<<"cpuLoop["<<Id<<"]: invoke new GpuKernel: gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
			std::cout<<"cpuLoop["<<Id<<"]: invoke new GpuKernel: nRowsLeft: "<<FRAME(nRowsLeft)<<std::endl;
			std::cout<<"cpuLoop["<<Id<<"]: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif
				}
			}
		}

	}
	int nCmpRep = FRAME(nCmpRep);
	for (size_t j=0;j<nCmpRep;++j){
		computeInner_stencil2d_v2(dst,src,nRows_bk,nCols);
	}
#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"dst:"<<std::endl;
//	std::cout<<std::setprecision(3)<<std::endl;
//	int tr = (nRows_bk<10)?nRows_bk:10;
//	int tc = (nCols<10)?nCols:10;
//	for(size_t i=0;i<tr;++i){
//		for (size_t j=0;j<tc;++j){
//			std::cout<<dst[i*nCols+j]<<",";
//		}
//		std::cout<<"\n";
//	}
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
	
	}else{
		
		size_t	gpuPos = FRAME(gpuPos);
		size_t	cpuPos = FRAME(cpuPos);
	
		uint64_t nRowsGpu = FRAME(nRowsGpu);
		uint64_t nRowsCpu = FRAME(nRowsCpu);
		uint64_t nCols = FRAME(nCols);
		uint64_t nRows = FRAME(nRows);

		double *dst = FRAME(New);
		double *d_dst = FRAME(d_dst);
		double *d_src = FRAME(d_src);
		double cmCpu = FRAME(cmCpu);
		double cmGpu = FRAME(cmGpu);
		
		if((cudaSuccess == cudaStreamQuery(NULL) )&&(__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),false,true))) {

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync: gpu kernel finish and invode a new gpu kernel"<<std::endl;
#endif

			double d_size = sizeof(double)*(nRowsGpu-2)*nCols;
			cudaError err1,err2,err3;
			err1 = cudaMemcpy(dst+(gpuPos+1)*nCols,d_dst+nCols, d_size, cudaMemcpyDeviceToHost);
			err2 = cudaFree(d_dst);
			err3 = cudaFree(d_src);
		
			if(err1!=cudaSuccess){
				std::cout<<"cpu sync: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err1)<<std::endl;
				exit(-1);
			}

			if(err2!=cudaSuccess){
				std::cout<<"CpuSync: cuda mem free d_dst: "<<cudaGetErrorString(err2)<<std::endl;
				exit(-1);
			}

			if(err3!=cudaSuccess){
				std::cout<<"CpuSync: cuda mem free d_src: "<<cudaGetErrorString(err2)<<std::endl;
				exit(-1);
			}

			size_t gpu_mem_total_t = 0;
			size_t gpu_mem_avail_t = 0;
			size_t gpu_mem_valid_t = 0;
			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
			gpu_mem_valid_t = gpu_mem_avail_t -XMB;

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"CpuSync: gpu memory total: "<<gpu_mem_total_t<<std::endl;
			std::cout<<"CpuSync: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
#endif		
			
			double	req_size = sizeof(double)* nRowsLeft*nCols*2;
			if(req_size < gpu_mem_valid_t){
				
				FRAME(gpuPos) = nRows-nRowsLeft;
				FRAME(nRowsGpu) = nRowsLeft;
				FRAME(nRowsLeft)=0;
				//INCR(Swap);
				SYNC(GpuKernel);
			
#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync: invoke new GpuKernel: gpuPos: "<<FRAME(gpuPos)<<std::endl;
			std::cout<<"cpuSync: invoke new GpuKernel: gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
			std::cout<<"cpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif

			}else if(req_size < 2*gpu_mem_valid_t){
				
				FRAME(gpuPos) = nRows-nRowsLeft;
				FRAME(nRowsGpu) = nRowsLeft/2;
				FRAME(cpuPos) = FRAME(gpuPos) + FRAME(nRowsGpu)-2;
				FRAME(nRowsCpu) = nRowsLeft - FRAME(nRowsGpu)+2; 
				FRAME(nRowsLeft) = 0;
				for(size_t i =0; i<FRAME(nCPU);++i){
					SYNC(CpuLoop[i]);
				}
				INCR(Swap);
				SYNC(GpuKernel);


#ifdef CUDA_DARTS_DEBUG
				std::cout<<"cpuSync: invoke new GpuKernel: gpuPos: "<<FRAME(gpuPos)<<std::endl;
				std::cout<<"cpuSync: invoke new CpuLoop: cpuPos: "<<FRAME(cpuPos)<<std::endl;
				std::cout<<"cpuSync: invoke new GpuKernel: gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
				std::cout<<"cpuSync: invoke new CpuLoop: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
				std::cout<<"cpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif


			}else{
				FRAME(gpuPos) = nRows-nRowsLeft;
				nRowsGpu = floor(gpu_mem_valid_t/(nCols*2*sizeof(double))); 
				FRAME(nRowsGpu) = nRowsGpu;
				FRAME(cpuPos) = FRAME(gpuPos) + FRAME(nRowsGpu)-2;
				nRowsCpu = nRowsGpu*(cmCpu/cmGpu);
				nRowsLeft = nRowsLeft - nRowsGpu - nRowsCpu+4;
				if (nRowsLeft< (1/cmGpu)*nRowsGpu){
						nRowsCpu = nRows-FRAME(cpuPos);
						nRowsLeft = 0;
				}
								
				FRAME(nRowsCpu) = nRowsCpu;
				FRAME(nRowsLeft) = nRowsLeft;

				for(size_t i =0; i<FRAME(nCPU);++i){
					SYNC(CpuLoop[i]);
				}
				INCR(Swap);
				SYNC(GpuKernel);

#ifdef CUDA_DARTS_DEBUG
				std::cout<<"cpuSync new GpuKernel: gpuPos: "<<FRAME(gpuPos)<<std::endl;
				std::cout<<"cpuSync new CpuLoop: cpuPos: "<<FRAME(cpuPos)<<std::endl;
				std::cout<<"cpuSync new GpuKernel: gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
				std::cout<<"cpuSync new CpuLoop: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
				std::cout<<"cpuSync new CpuLoop: reset rowsLeft: "<<FRAME(nRowsLeft)<<std::endl;
				std::cout<<"cpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;

#endif
			}


		}else{

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"CpuSync: invoke new CpuLoop!"<<std::endl;
#endif
			FRAME(cpuPos) = nRows-nRowsLeft;
			if(nRowsLeft < nRowsCpu*(cmGpu/cmCpu)){
				nRowsCpu = nRowsLeft;
				nRowsLeft = 0;
			}else{
				nRowsLeft = nRowsLeft - nRowsCpu+2;
			}
			FRAME(nRowsCpu) = nRowsCpu;
			FRAME(nRowsLeft) = nRowsLeft;

			for(size_t i =0; i<FRAME(nCPU);++i){
				SYNC(CpuLoop[i]);
			}
//			std::cout<<"CpuSync, nRowsLeft:"<<FRAME(nRowsLeft)<<std::endl;

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync new CpuLoop: cpuPos: "<<FRAME(cpuPos)<<std::endl;
			std::cout<<"cpuSync new CpuLoop: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
#endif
		
		}

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernel: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif
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

		SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
		if(ts!=0){
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}
		}else{
			SYNC(sync);
		}

	}else if(GpuRatio ==1.0){

		SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
		if(ts!=0){
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

			uint64_t	gpuPos = FRAME(gpuPos);
			uint64_t	cpuPos = FRAME(cpuPos);
			
			uint64_t nRows   = FRAME(nRows);
			uint64_t nCols   = FRAME(nCols);
			uint64_t nRowsGpu = FRAME(nRowsGpu);
			uint64_t nRowsCpu = FRAME(nRowsCpu);
			
			double d_size = sizeof(double)*(nRowsGpu-2)*nCols;
			
			cudaError err0, err1,err2,err3,err4;
			err0 = cudaDeviceSynchronize();

			err1 = cudaMemcpy(dst+(gpuPos+1)*nCols, d_dst+nCols,d_size , cudaMemcpyDeviceToHost);
		
			err2	= cudaFree(d_dst);
			err3	= cudaFree(d_src);

			if(err0!=cudaSuccess){
				std::cout<<"swap: cuda deviceSynchronize: "<<cudaGetErrorString(err0)<<std::endl;
				exit(-1);
			}
			if(err1!=cudaSuccess){
				std::cout<<"swap: cuda Memcpy DeviceToHost: "<<cudaGetErrorString(err1)<<std::endl;
				exit(-1);
			}
			if(err2!=cudaSuccess){
				std::cout<<"swap: cudaFree d_dst: "<<cudaGetErrorString(err2)<<std::endl;
				exit(-1);
			}

			if(err3!=cudaSuccess){
				std::cout<<"swap: cudaFree d_src: "<<cudaGetErrorString(err3)<<std::endl;
				exit(-1);
			}

		if(ts!=0){
			FRAME(gpuPos)	=FRAME(gpuPosInit);
			FRAME(cpuPos)	=FRAME(cpuPosInit);
			FRAME(nRowsGpu)	=FRAME(nRowsGpuInit);
			FRAME(nRowsCpu) =FRAME(nRowsCpuInit);
			FRAME(nRowsLeft)=FRAME(nRowsLeftInit);
			if(FRAME(nRowsLeft)!=0){
				FRAME(CpuIvGpu) = true;
			}
#ifdef CUDA_DARTS_DEBUG
			std::cout<<"swap: reset gpuPos: "<<FRAME(gpuPos)<<std::endl;
			std::cout<<"swap: reset cpuPos: "<<FRAME(cpuPos)<<std::endl;
			std::cout<<"swap: reset gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
			std::cout<<"swap: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
			std::cout<<"swap: reset cpu rowsLeft: "<<FRAME(nRowsLeft)<<std::endl;
#endif


#ifdef CUDA_DARTS_DEBUG
			size_t gpu_mem_total_t = 0;
			size_t gpu_mem_avail_t = 0;
			size_t gpu_mem_valid_t = 0;
			err4 = cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

			gpu_mem_valid_t = gpu_mem_avail_t -XMB;
			
			if(err4!=cudaSuccess){
				std::cout<<"swap: cuda get mem info: "<<cudaGetErrorString(err4)<<std::endl;
				exit(-1);
			}

			std::cout<<"swap: gpu memory total: "<<gpu_mem_total_t<<std::endl;
			std::cout<<"swap: gpu memory available: "<<gpu_mem_avail_t<<std::endl;

#endif		
		
			SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
			SYNC(GpuKernel);
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}

		}else{

			SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
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
#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"Sync: dst:"<<std::endl;
//	std::cout<<std::setprecision(3)<<std::endl;
//	int tr = 10;
//	int tc = 10;
//	double *lastdst=FRAME(Initial);
//	for(size_t i=0;i<tr;++i){
//		for (size_t j=0;j<tc;++j){
//			std::cout<<lastdst[i*FRAME(nCols)+j]<<",";
//		}
//		std::cout<<"\n";
//	}

#endif

	SIGNAL(signalUp);
    EXIT_TP();
}





