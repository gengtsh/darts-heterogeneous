#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "StencilTP.h"
#include "stencil.h"
#include "StencilCPUKernel.h"

#include <cassert>
#include <pthread.h>
pthread_mutex_t mutex;
pthread_mutex_t mutex2;
//#include <sstream>
#include <iostream>


////	pthread_mutex_lock(&mutex);
////	pthread_mutex_unlock(&mutex);

extern "C"
void
Stencil2D4ptGpuKernelWithAllTimeStepsCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelWithAllTimeSteps!"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(GpuKernelWithAllTimeSteps);

	double d_size;	
	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;


	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t nRowsGpu = FRAME(nRowsGpu);

	d_size = sizeof(double)*nRowsGpu*nCols;
	
	int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRowsGpu/tile_y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = gpuPos*nCols;

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuKernelWithAllTimeSteps: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nRowsGpu:"<<nRowsGpu<<std::endl;

	std::cout<<"GpuKernelWithAllTimeSteps: blockDimx:"<<blockDimx<<std::endl;
	
	std::cout<<"GpuKernelWithAllTimeSteps: gridDimx="<<gridDimx<<",gridDimy="<<gridDimy<<std::endl;
#endif
	
	int64_t d_size_sharedCols = sizeof(double)*nRowsGpu*gridDimx*2;
	int64_t d_size_sharedRows = sizeof(double)*nCols*gridDimy*2;

	cudaError err1,err2,err3,err4,err5,err6,err7,err8,err9,err10, err11;
	err1 = cudaMalloc( (void **) &d_dst, d_size);

	err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
	err9 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	if(err9!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc9: "<<cudaGetErrorString(err9)<<std::endl;
		exit(-1);
	}

#endif

	FRAME(d_dst) = d_dst;
		
	err3 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);


#ifdef CUDA_ERROR_CHECKING
	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

#endif

	int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
	int blockDimy_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	int gridDimy_rows = std::ceil(1.0*nRowsGpu/tile_y);


	int blockDimx_cols = 1 ;
	int blockDimy_cols = (nRowsGpu>NUM_THREADS)?NUM_THREADS:nRows;
	int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRowsGpu/blockDimy_cols);

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuKernelWithAllTimeSteps: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimy_rows:"<<gridDimy_rows<<std::endl;

	std::cout<<"GpuKernelWithAllTimeSteps: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimy_cols:"<<gridDimy_cols<<std::endl;

#endif
	
	dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);

	dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);


	size_t ts = FRAME(ts);
	while(ts-- >0){
		gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,nRowsGpu, nCols);
		gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,nRowsGpu, nCols);
		gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,nRowsGpu,nCols);
	}
	
	err5 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err5!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda deviceSynchronize: "<<cudaGetErrorString(err5)<<std::endl;
		exit(-1);
	}
#endif
        
#ifdef VERIFICATION
        if(FRAME(ts)%2==0){
	    SWAP_PTR(&h_dst ,&h_src);
        }
#endif

	//copy from GPU  to CPU
        err6=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
	if(err6!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err6)<<std::endl;
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
	err8 = cudaFree(d_sharedCols);
	err11 = cudaFree(d_sharedRows);
	
#ifdef CUDA_ERROR_CHECKING
	if(err7!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_dst: "<<cudaGetErrorString(err7)<<std::endl;
		exit(-1);
	}

	if(err8!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err8)<<std::endl;
		exit(-1);
	}

	if(err11!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err11)<<std::endl;
		exit(-1);
	}
#endif

	SYNC(sync);

	EXIT_TP();

}

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

	double *d_sharedCols ;
	double *d_sharedRows ;
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);

	double d_size;	

	d_size = sizeof(double)*nRowsGpu*nCols;
	
	int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRowsGpu/tile_y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = gpuPos*nCols;

	int64_t d_size_sharedCols = sizeof(double)*nRowsGpu*gridDimx*2;
	int64_t d_size_sharedRows = sizeof(double)*nCols*gridDimy*2;
#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernel: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernel: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernel: require memory size:"<<d_size + d_size_sharedCols + d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernel: Gpupos:"<<gpuPos<<std::endl;
	std::cout<<"GpuKernel: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernel: nRowsGpu:"<<nRowsGpu<<std::endl;

#endif

	cudaError err1,err2,err3,err4,err5,err6,err7,err8,err9,err10, err11;

    err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
    err9 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
#ifdef CUDA_ERROR_CHECKING
    if(err2!=cudaSuccess){
    	std::cout<<"GpuKernel: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
    	exit(-1);
    }
    if(err9!=cudaSuccess){
    	std::cout<<"GpuKernel: cuda malloc9: "<<cudaGetErrorString(err9)<<std::endl;
    	exit(-1);
    }
#endif
    
    if(FRAME(invokeStreams)==false){

        if(FRAME(ts)==FRAME(tsInit) ){
            err1 = cudaMalloc( (void **) &d_dst, d_size);
            FRAME(d_dst) = d_dst;
            err3 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);

#ifdef CUDA_ERROR_CHECKING
            if(err1!=cudaSuccess){
                std::cout<<"GpuKernel: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
    	        exit(-1);
            }
            if(err3!=cudaSuccess){
    	        std::cout<<"GpuKernel: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err3)<<std::endl;
    	        exit(-1);
            }
#endif
        }
    
    	d_dst = FRAME(d_dst);
    	int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
    	int blockDimy_rows = 1;
    	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
    	int gridDimy_rows = std::ceil(1.0*nRowsGpu/tile_y);
    
    
    	int blockDimx_cols = 1 ;
    	int blockDimy_cols = (nRowsGpu>NUM_THREADS)?NUM_THREADS:nRows;
    	int gridDimx_cols = gridDimx;
    	int gridDimy_cols = std::ceil(1.0*nRowsGpu/blockDimy_cols);
    
    
    	dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
    	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    
    	dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
    	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
    
    
    	gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,nRowsGpu, nCols);
    	gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,nRowsGpu, nCols);
    	gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,nRowsGpu,nCols);
	
    }else{
        
#ifdef CUDA_DARTS_DEBUG
    	std::cout<<"GpuKernel multiple streams: FRAME(invokeStreams): "<<FRAME(invokeStreams)<<std::endl;
#endif
        int nStream = FRAME(nStream);
        d_size = d_size+sizeof(double)*nStream*nCols*2; 
        
#ifdef CUDA_DARTS_DEBUG		
			std::cout<<"GpuKernel multiple streams: required memory size: "<<d_size+d_size_sharedCols+ d_size_sharedRows<<",stream computation size:"<<d_size<<std::endl;
#endif
        err1 = cudaMalloc( (void **) &d_dst, d_size);
#ifdef CUDA_ERROR_CHECKING
        if(err1!=cudaSuccess){
            std::cout<<"GpuKernel multiple streams: cuda malloc d_dst: "<<cudaGetErrorString(err1)<<std::endl;
            exit(-1);
        }
#endif
        
        FRAME(d_dst) = d_dst;
         
    	int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
    	int blockDimy_rows = 1;
    	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
    	int gridDimy_rows; 
        //int gridDimy_rows = std::ceil(1.0*nRowsGpu/tile_y);
    
    
    	int blockDimx_cols = 1 ;
    	int blockDimy_cols;
        //int blockDimy_cols = (nRowsGpu>NUM_THREADS)?NUM_THREADS:nRowsGpu;
    	int gridDimx_cols = gridDimx;
        int gridDimy_cols; 
    	//int gridDimy_cols = std::ceil(1.0*nRowsGpu/blockDimy_cols);
    	
        //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
    	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    
    	//dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
    	//dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
    
        uint64_t d_size_stream; 
        uint64_t nRowsStream;
        int nTile_y = (nRowsGpu/(tile_y*nStream));
        int chunk=nTile_y*tile_y;
        int chunk2 = chunk+2;
        uint64_t h_pos;
        uint64_t d_pos;
        uint64_t pos2;
        uint64_t pos3=0;
        for(size_t i =0; i<nStream;++i){  

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernel multiple streams: kernel5 stream error: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
#ifdef CUDA_DARTS_DEBUG
        	std::cout<<"GpuKernel multiple streams: invoke stream["<<i<<"]"<<std::endl;
#endif
            nRowsStream = (i==(nStream-1))?(nRowsGpu-i*chunk):(chunk+2);
            d_size_stream= sizeof(double)*nCols*nRowsStream;
            pos2 = i*chunk*nCols;
            h_pos = pos1+pos2;
            d_pos = i*chunk2*nCols;
#ifdef CUDA_DARTS_DEBUG
        	std::cout<<"GpuKernel multiple streams: stream["<<i<<"]:size:"<<d_size_stream<<",d_pos:"<<d_pos/nCols<<",h_pos:"<<h_pos/nCols<<",nRowsStream:"<<nRowsStream<<std::endl;
#endif
            err5 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,FRAME(stream[i]));
#ifdef CUDA_ERROR_CHECKING
            if(err5!=cudaSuccess){
                std::cout<<"GpuKernel multiple streams: cuda MemcpyAsync 5 from host to device: "<<cudaGetErrorString(err5)<<std::endl;
                exit(-1);
            }
#endif

            gridDimy_rows = std::ceil(1.0*nRowsStream/tile_y);
            dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
            gpu_kernel5_stream_cp_rows(FRAME(stream[i]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+i*nTile_y*2*nCols, tile_y,nRowsGpu, nCols);
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernel multiple streams: kernel5 cuda cp rows: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
            blockDimy_cols = (nRowsStream>NUM_THREADS)?NUM_THREADS:nRowsStream;
            gridDimy_cols = std::ceil(1.0*nRowsStream/blockDimy_cols);
    	    dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
    	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);

#ifdef CUDA_DARTS_DEBUG
        	std::cout<<"GpuKernel multiple streams cp cols: stream["<<i<<"]:cols begin: "<<i*chunk<<",blockDimy_cols:"<<blockDimy_cols<<", gridDimy_cols: "<<gridDimy_cols<<std::endl;
#endif
            gpu_kernel5_stream_cp_cols(FRAME(stream[i]),dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+i*chunk, d_sharedRows, tile_x,tile_x,nRowsStream, nCols);
    	     
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernel multiple streams: cuda kernel5 cp cols: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
	        int gridDimy_stream = std::ceil(1.0*nRowsStream/tile_y);
	        dim3 dimGrid_stream(gridDimx,gridDimy_stream);

#ifdef CUDA_DARTS_DEBUG
        	std::cout<<"GpuKernel multiple streams computation: stream["<<i<<"]:device begin: "<<d_pos/nCols<<",host begin:"<<h_pos/nCols<<",nRowsStream:"<<nRowsStream<<std::endl;
#endif
            gpu_kernel5_stream(FRAME(stream[i]) ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+i*chunk,d_sharedRows+i*nTile_y*2*nCols,tile_y,nRowsStream,nCols);

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernel multiple streams: cuda kernel5 compute: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
	        err6=cudaMemcpyAsync(h_dst+h_pos+nCols, d_dst+d_pos+nCols,d_size_stream-(nCols)*2*sizeof(double), cudaMemcpyDeviceToHost,FRAME(stream[i]));

#ifdef CUDA_ERROR_CHECKING
        
            if(err6!=cudaSuccess){
                std::cout<<"GpuKernel multiple streams: cuda Memcpy Async 6  from Device to Host : "<<cudaGetErrorString(err6)<<std::endl;
                exit(-1);
            }
#endif
        
        }
    }
	err5 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err5!=cudaSuccess){
		std::cout<<"GpuKernel: cuda deviceSynchronize: "<<cudaGetErrorString(err5)<<std::endl;
		exit(-1);
	}
#endif
       
	err8 = cudaFree(d_sharedCols);
	err11 = cudaFree(d_sharedRows);
	if(FRAME(invokeStreams)==true){
	    err2 = cudaFree(d_dst);
#ifdef CUDA_ERROR_CHECKING
	    if(err2!=cudaSuccess){
		    std::cout<<"GpuKernel: cuda memcpy free d_dst: "<<cudaGetErrorString(err2)<<std::endl;
	    	exit(-1);
	    }
#endif 
    }
#ifdef CUDA_ERROR_CHECKING
	if(err8!=cudaSuccess){
		std::cout<<"GpuKernel: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err8)<<std::endl;
		exit(-1);
	}

	if(err11!=cudaSuccess){
		std::cout<<"GpuKernel: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err11)<<std::endl;
		exit(-1);
	}
#endif

    	__sync_bool_compare_and_swap(&FRAME(CpuIvGpu),true,false);
    
#ifdef CUDA_DARTS_DEBUG
    	std::cout<<"FRAME(CpuIvGpu): "<<FRAME(CpuIvGpu)<<std::endl;
#endif
        
        SYNC(Swap);	

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernel: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif

	EXIT_TP();
}


extern "C"
void
Stencil2D4ptGpuKernelPureGpuWithStreamsCD::fire(void)
{

	LOAD_FRAME(StencilTP);
    RESET(GpuKernelPureGpuWithStreams);
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelPureGpuWithStreams"<<std::endl;
#endif

	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    uint32_t nGPU = FRAME(nGPU);
    uint64_t gpuPos = FRAME(gpuPos); 
    int nStream = FRAME(nStream);
    int vnStream = nStream*nGPU;
    
    int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
    
    int nTile_y = nRows/(tile_y * vnStream);

    int chunk = nTile_y*tile_y;
    int chunk2= chunk+2;
    
    int nRowsGpuBlock = nStream*chunk2 + nRows-nGPU*nStream*chunk;
    int64_t nRowsGpuStream;
	int64_t d_size_stream;
   
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRowsGpuBlock/tile_y); 
    
	dim3 dimBlock(blockDimx,blockDimy);
//	dim3 dimGrid(gridDimx,gridDimy);
    
    int64_t d_size = sizeof(double)*nRowsGpuBlock*nCols;
    int64_t d_size_sharedCols = sizeof(double) * nRowsGpuBlock*gridDimx*2 ;
    int64_t d_size_sharedRows = sizeof(double) * nCols* gridDimy*2;

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernelWithStream : require memory size:"<<d_size + d_size_sharedCols + d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: Gpupos:"<<FRAME(gpuPos)<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: nRowsGpuBlock:"<<nRowsGpuBlock<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams: chunk:"<<chunk<<std::endl;
#endif
    
    cudaError err1,err2,err3,err4;
	
    err1 = cudaMalloc( (void **) &d_dst, d_size);
    err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
    err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
    
#ifdef CUDA_ERROR_CHECKING
    if(err1!=cudaSuccess){
        std::cout<<"GpuKernelPureGpuWithStreams: cuda malloc d_dst: "<<cudaGetErrorString(err1)<<std::endl;
        exit(-1);
    }
    if(err2!=cudaSuccess){
        std::cout<<"GpuKernelPureGpuWithStreams: cuda mallock d_sharedRows: "<<cudaGetErrorString(err2)<<std::endl;
        exit(-1);
    }
    
    if(err3!=cudaSuccess){
        std::cout<<"GpuKernelPureGpuWithStreams: cuda mallock d_sharedCols: "<<cudaGetErrorString(err3)<<std::endl;
        exit(-1);
    }
#endif
    
    int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
    int blockDimy_rows = 1;
    int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
    int gridDimy_rows; 
    
    int blockDimx_cols = 1 ;
    int blockDimy_cols;
    int gridDimx_cols = gridDimx;  
    int gridDimy_cols; 
    
    dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);

    uint64_t h_pos;
    uint64_t d_pos;
    uint64_t pos0 = gpuPos*nCols;
    
    size_t ts = FRAME(ts);
    while(ts-- >0){
        for (size_t i = 0; i<nGPU; ++i){
            for (size_t j =0; j<nStream;++j){
                int ps = i*nStream+j;
                nRowsGpuStream = ((i==(nGPU-1))&&(j==(nStream-1)))? (nRows-ps*chunk) :chunk2;
                h_pos = pos0+ps*chunk*nCols;
                d_pos = j*chunk2*nCols; 

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelPureGpuWithStreams multiple streams: kernel5 stream error: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif

#ifdef CUDA_DARTS_DEBUG
        	    std::cout<<"GpuKernelPureGpuWithStreams multiple streams: invoke stream["<<j<<"]"<<",ps:"<<ps<<std::endl;
#endif

                d_size_stream = sizeof(double)*nCols*nRowsGpuStream;
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"GpuKernelPureGpuWithStreams multiple streams: stream["<<j<<"]"<<",ps:"<<ps<<", size:"<<d_size_stream<<",d_pos:"<<d_pos/nCols<<",h_pos:"<<h_pos/nCols<<",nRowsGpuStream:"<<nRowsGpuStream<<std::endl;
#endif
                err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
                if(err1!=cudaSuccess){
                    std::cout<<"GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                    exit(-1);
                }
#endif
            
                gridDimy_rows = std::ceil(1.0*nRowsGpuStream/tile_y);
                dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
                gpu_kernel5_stream_cp_rows(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+j*nTile_y*2*nCols, tile_y,nRowsGpuStream, nCols);
            
#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp rows: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
                blockDimy_cols = (nRowsGpuStream>NUM_THREADS)?NUM_THREADS:nRowsGpuStream;
                gridDimy_cols = std::ceil(1.0*nRowsGpuStream/blockDimy_cols);
        	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
                dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
                int addrCol = j*chunk2*2*gridDimx_cols; 
                gpu_kernel5_stream_cp_cols(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+addrCol, d_sharedRows, tile_x,tile_x,nRowsGpuStream, nCols);

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp cols: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
	            int gridDimy_stream = std::ceil(1.0*nRowsGpuStream/tile_y);
	            dim3 dimGrid_stream(gridDimx,gridDimy_stream);
                gpu_kernel5_stream(FRAME(stream[j]) ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+addrCol,d_sharedRows+j*nTile_y*2*nCols,tile_y,nRowsGpuStream,nCols);
            
#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda computation: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
	            err3=cudaMemcpyAsync(h_dst+h_pos+nCols, d_dst+d_pos+nCols,d_size_stream-(nCols)*2*sizeof(double), cudaMemcpyDeviceToHost,FRAME(stream[j]));

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 Asyn Memory copy from device to host: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            }
        }
   
	    err4 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
	    if(err4!=cudaSuccess){
		    std::cout<<"GpuKernelPureGpuWithStreams: cuda deviceSynchronize: "<<cudaGetErrorString(err4)<<std::endl;
		    exit(-1);
	    }
#endif
	    SWAP_PTR(&h_dst ,&h_src);
    }
	
    
	err4 = cudaDeviceSynchronize();
    err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedRows);
    err3 = cudaFree(d_sharedCols);

#ifdef CUDA_ERROR_CHECKING
	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams: cuda deviceSynchronize: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}

    if(err1!=cudaSuccess){
	    std::cout<<"GpuKernelPureGpuWithStreams: cuda memcpy free d_dst: "<<cudaGetErrorString(err1)<<std::endl;
    	exit(-1);
    }

    if(err2!=cudaSuccess){
	    std::cout<<"GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err2)<<std::endl;
    	exit(-1);
    }

    if(err3!=cudaSuccess){
	    std::cout<<"GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err3)<<std::endl;
    	exit(-1);
    }
#endif
    if(FRAME(tsInit%2)){
	    SWAP_PTR(&h_dst ,&h_src);
    }

    SYNC(sync);
	EXIT_TP();
}



extern "C"
void
Stencil2D4ptGpuKernelHybridWithStreamsCD::fire(void)
{

	LOAD_FRAME(StencilTP);
    RESET(GpuKernelHybridWithStreams);
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelHybridWithStreams"<<std::endl;
#endif

	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    uint32_t nGPU = FRAME(nGPU);
    uint64_t gpuPos = FRAME(gpuPos);
    uint64_t nRowsGpu = FRAME(nRowsGpu);
    int nStream = FRAME(nStream);
    int vnStream = nStream*nGPU;
    
    int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
    
    int nTile_y = nRowsGpu/(tile_y * vnStream);

    int chunk = nTile_y*tile_y;
    int chunk2= chunk+2;
    
    int nRowsGpuBlock = nStream*chunk2 + nRowsGpu-nGPU*nStream*chunk;
    int64_t nRowsGpuStream;
	int64_t d_size_stream;
   
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*nRowsGpuBlock/tile_y); 
    
	dim3 dimBlock(blockDimx,blockDimy);
//	dim3 dimGrid(gridDimx,gridDimy);
    
    int64_t d_size = sizeof(double)*(nRowsGpuBlock)*nCols;
    int64_t d_size_sharedCols = sizeof(double) * (nRowsGpuBlock+3)*gridDimx*2 ;
    int64_t d_size_sharedRows = sizeof(double) * (nCols +NUM_THREADS )* gridDimy*2;

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: require memory size:"<<d_size + d_size_sharedCols + d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: Gpupos:"<<FRAME(gpuPos)<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: nRowsGpuBlock:"<<nRowsGpuBlock<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: chunk:"<<chunk<<std::endl;
#endif
    
    cudaError err1,err2,err3,err4;
	
    err1 = cudaMalloc( (void **) &d_dst, d_size);
    err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
    err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
    
#ifdef CUDA_ERROR_CHECKING
    if(err1!=cudaSuccess){
        std::cout<<"GpuKernelHybridWithStreams: cuda malloc d_dst: "<<cudaGetErrorString(err1)<<std::endl;
        exit(-1);
    }
    if(err2!=cudaSuccess){
        std::cout<<"GpuKernelHybridWithStreams: cuda mallock d_sharedRows: "<<cudaGetErrorString(err2)<<std::endl;
        exit(-1);
    }
    
    if(err3!=cudaSuccess){
        std::cout<<"GpuKernelHybridWithStreams: cuda mallock d_sharedCols: "<<cudaGetErrorString(err3)<<std::endl;
        exit(-1);
    }
#endif
    
    int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
    int blockDimy_rows = 1;
    int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
    int gridDimy_rows; 
    
    int blockDimx_cols = 1 ;
    int blockDimy_cols;
    int gridDimx_cols = gridDimx;  
    int gridDimy_cols; 
    
    dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);

    uint64_t h_pos;
    uint64_t d_pos;
    uint64_t pos0 = gpuPos*nCols;
    
    for (size_t i = 0; i<nGPU; ++i){
        for (size_t j =0; j<nStream;++j){
            int ps = i*nStream+j;
            nRowsGpuStream = ((i==(nGPU-1))&&(j==(nStream-1)))? (nRowsGpu-ps*chunk) :chunk2;
            h_pos = pos0+ps*chunk*nCols;
            d_pos = j*chunk2*nCols; 

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelHybridWithStreams multiple streams: kernel5 stream error: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif

#ifdef CUDA_DARTS_DEBUG
    	    std::cout<<"GpuKernelHybridWithStreams multiple streams: invoke stream["<<j<<"]"<<",ps:"<<ps<<std::endl;
#endif

            d_size_stream = sizeof(double)*nCols*nRowsGpuStream;
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams multiple streams: stream["<<j<<"]"<<",ps:"<<ps<<", size:"<<d_size_stream<<",d_pos:"<<d_pos/nCols<<",h_pos:"<<h_pos/nCols<<",nRowsGpuStream:"<<nRowsGpuStream<<std::endl;
#endif
            err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
            if(err1!=cudaSuccess){
                std::cout<<"GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                exit(-1);
            }
#endif
        
            gridDimy_rows = std::ceil(1.0*nRowsGpuStream/tile_y);
            dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
            gpu_kernel5_stream_cp_rows(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+j*nTile_y*2*nCols, tile_y,nRowsGpuStream, nCols);
        
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp rows: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        
            blockDimy_cols = (nRowsGpuStream>NUM_THREADS)?NUM_THREADS:nRowsGpuStream;
            gridDimy_cols = std::ceil(1.0*nRowsGpuStream/blockDimy_cols);
    	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
            dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
            int addrCol = j*chunk2*2*gridDimx_cols; 
            gpu_kernel5_stream_cp_cols(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+addrCol, d_sharedRows, tile_x,tile_x,nRowsGpuStream, nCols);

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp cols: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        
	        int gridDimy_stream = std::ceil(1.0*nRowsGpuStream/tile_y);
	        dim3 dimGrid_stream(gridDimx,gridDimy_stream);
            gpu_kernel5_stream(FRAME(stream[j]) ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+addrCol,d_sharedRows+j*nTile_y*2*nCols,tile_y,nRowsGpuStream,nCols);
        
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda computation: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        
	        err3=cudaMemcpyAsync(h_dst+h_pos+nCols, d_dst+d_pos+nCols,d_size_stream-(nCols)*2*sizeof(double), cudaMemcpyDeviceToHost,FRAME(stream[j]));

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 Asyn Memory copy from device to host: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        }
    }
	
    
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 finish and get last error! "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
	err4 = cudaDeviceSynchronize();
    err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedRows);
    err3 = cudaFree(d_sharedCols);

#ifdef CUDA_ERROR_CHECKING
	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams: cuda deviceSynchronize: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}

    if(err1!=cudaSuccess){
	    std::cout<<"GpuKernelHybridWithStreams: cuda memcpy free d_dst: "<<cudaGetErrorString(err1)<<std::endl;
    	exit(-1);
    }

    if(err2!=cudaSuccess){
	    std::cout<<"GpuKernelHybridWithStreams: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err2)<<std::endl;
    	exit(-1);
    }

    if(err3!=cudaSuccess){
	    std::cout<<"GpuKernelHybridWithStreams: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err3)<<std::endl;
    	exit(-1);
    }
#endif
    
    __sync_bool_compare_and_swap(&FRAME(GpuFinish),false,true);
    uint64_t nRowsLeft = FRAME(nRowsLeft);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: GPU finish! and nRowsLeft: "<<nRowsLeft<<std::endl;
#endif
    ++FRAME(gpuCnt); 
    if(nRowsLeft == 0){
        SYNC(Swap);
    }else{

        if(__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,false)&&(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false) )){ 
            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));
            double gpuStepR = FRAME(gpuStepR);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: invoke new GPU"<<std::endl;
#endif
            
            pthread_mutex_lock(&mutex);
            nRowsLeft = FRAME(nRowsLeft);
            if(nRowsLeft==0){
                SYNC(Swap);
            }else{
             
                if(nRowsLeft <= FRAME(lastCnt)){
                    FRAME(nRowsGpu) = nRowsLeft;
                    FRAME(nRowsLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*nRowsGpu;
                        t1 = (t2<FRAME(nRowsGpuMax))?t2:FRAME(nRowsGpuMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = nRowsGpu;
                    }else{
                        t1 = (1+gpuStepR)*nRowsGpu;
                    }

                    if(nRowsLeft<t1){
                        FRAME(nRowsGpu) = nRowsLeft*rt;
                    }else{
                        FRAME(nRowsGpu) = t1;
                    }
               
                    FRAME(nRowsLeft) = nRowsLeft - FRAME(nRowsGpu)+2;
                }
                
                FRAME(nGPU)=2;
                FRAME(gpuPos) = nRows-nRowsLeft ;
                __sync_synchronize();

                SYNC(GpuKernelHybridWithStreams);
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, gpuPos: "<<FRAME(gpuPos)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, nRowsGpu: "<<FRAME(nRowsGpu)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, nRowsLeft: "<<FRAME(nRowsLeft)<<std::endl;

#endif
        }
    }

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernelHybridWithStreams: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
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
		
	err3 = cudaMemcpy(d_dst, h_dst+pos1, d_size, cudaMemcpyHostToDevice);


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

	gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,nRows_bk,nCols);
	
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
	err8 = cudaFree(d_sharedCols);
	err11 = cudaFree(d_sharedRows);
	
#ifdef CUDA_ERROR_CHECKING
	if(err7!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_dst: "<<cudaGetErrorString(err7)<<std::endl;
		exit(-1);
	}

	if(err8!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err8)<<std::endl;
		exit(-1);
	}

	if(err11!=cudaSuccess){
		std::cout<<"GpuLoop: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err11)<<std::endl;
		exit(-1);
	}
#endif
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
	//const uint64_t nRows   = FRAME(nRows);
	//const uint64_t nCols   = FRAME(nCols);
	uint64_t	nRowsCpu = FRAME(nRowsCpu);	
	//uint64_t	nCPU = FRAME(nCPU);
	//uint64_t	chunk = (nRowsCpu )/nCPU;
	//uint64_t	cpuPos = FRAME(cpuPos);
	//size_t		pos1 = (cpuPos + chunk*Id)*nCols;
	//uint64_t	nRows_bk = (Id==(nCPU-1))? (nRowsCpu-chunk*Id-1):(chunk +1) ;
	//double *src = h_src+pos1;
	//double *dst = h_dst+pos1;

	//double *d_dst = FRAME(d_dst);
	//uint64_t nRowsLeft = FRAME(nRowsLeft);


    int *arrnEdge       = FRAME(arrnEdge);
    int *arrnCpuPos     = FRAME(arrnCpuPos);
    int *arrnCpuEdge    = FRAME(arrnCpuEdge);
    int *arrnCpuBlock   = FRAME(arrnCpuBlock);
    int *arrnCpuGrid    = FRAME(arrnCpuGrid);

	const int nCols		= FRAME(nCols);
	const int nRows		= FRAME(nRows);

	int cpuCols		= arrnCpuEdge[0];
	int cpuRows		= arrnCpuEdge[1];

    int cpuBlockDimx	= arrnCpuBlock[0];
    int cpuBlockDimy	= arrnCpuBlock[1];
    int cpuGridDimx		= arrnCpuGrid[0];
    int cpuGridDimy		= arrnCpuGrid[1];

    int cpuPosx         = arrnCpuPos[0]; 
    int cpuPosy         = arrnCpuPos[1];

    int  arrnCpuEdge2[DIM];

    calcEdge(arrnCpuEdge2,arrnCpuPos,arrnCpuEdge, arrnEdge, 1,0,0,2);
	
    int cpuCols2		= arrnCpuEdge2[0];
	int cpuRows2		= arrnCpuEdge2[1];
    
    int cpuPos          = cpuPosx+cpuPosy*nCols;
	
	int posy = Id/cpuGridDimx;
	int posx = Id - posy*cpuGridDimx;

	uint64_t pos = posy*cpuBlockDimy*nCols + posx*cpuBlockDimx;
	double *src = h_src+cpuPos+pos;
	double *dst = h_dst+cpuPos+pos;
	
	int nRowsChunk		= (((posy+1)*cpuBlockDimy)>=cpuRows2)?(cpuRows2-posy*cpuBlockDimy-1) : (cpuBlockDimy+1) ;
	int nColsChunk		= (((posx+1)*cpuBlockDimx)>=cpuCols2)?(cpuCols2-posx*cpuBlockDimx-1):(cpuBlockDimx+1) ;
#ifdef CUDA_DARTS_DEBUG

	pthread_mutex_lock(&mutex2);
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPos:"<<cpuPos<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPosx:"<<cpuPosx<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPosy:"<<cpuPosy<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posx:"<<posx<<std::endl;
    std::cout<<"CpuLoop37["<<Id<<"]: posy:"<<posy<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nRowsChunk:"<<nRowsChunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nColsChunk:"<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif

	computeBlock_stencil25(dst,src,nRows,nCols,nRowsChunk,nColsChunk);


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

extern "C"
void Stencil2D4ptCpuSyncCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync invoke!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
	uint64_t nRowsLeft = FRAME(nRowsLeft);
	//uint64_t nRows = FRAME(nRows);
	RESET(CpuSync);

    ++FRAME(cpuCnt);

    if(nRowsLeft==0.0){
		SYNC(Swap);
	
	}else{
		__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,true);
	
        size_t	gpuPos = FRAME(gpuPos);
		size_t	cpuPos = FRAME(cpuPos);
	
		uint64_t nRowsGpu = FRAME(nRowsGpu);
		uint64_t nRowsCpu = FRAME(nRowsCpu);
		uint64_t nCols = FRAME(nCols);
		uint64_t nRows = FRAME(nRows);

		double *dst = FRAME(New);
		double *d_dst = FRAME(d_dst);
		double cmCpu = FRAME(cmCpu);
		double cmGpu = FRAME(cmGpu);
        
        double cpuStepR = FRAME(cpuStepR);
        double gpuStepR = FRAME(gpuStepR);
       
        
        int *arrnCpuPos  = FRAME(arrnCpuPos);
        int *arrnCpuEdge = FRAME(arrnCpuEdge); 
        int *arrnCpuTile = FRAME(arrnCpuTile); 
        int *arrnCpuGridTileBase = FRAME(arrnCpuGridTileBase);
        int *arrnCpuBlock = FRAME(arrnCpuBlock);
        int *arrnCpuGrid  = FRAME(arrnCpuGrid);



        if(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false)) {

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync: gpu kernel finish and invode a new gpu kernel"<<std::endl;
#endif

            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            nRowsLeft = FRAME(nRowsLeft);
            if(nRowsLeft==0.0){
                SYNC(Swap);
            }else{
                
                if(nRowsLeft <= FRAME(lastCnt)){
                    FRAME(nRowsGpu) = nRowsLeft;
                    FRAME(nRowsLeft) = 0;
                }else{
                    uint64_t t1;
                  
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*nRowsGpu;
                        t1 = (t2<FRAME(nRowsGpuMax))?t2:FRAME(nRowsGpuMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = nRowsGpu;
                    }else{
                        t1 = (1+gpuStepR)*nRowsGpu;
                    }
                    
                    if(nRowsLeft<t1){
                        FRAME(nRowsGpu) = nRowsLeft*rt;
                    }else{
                        FRAME(nRowsGpu) = t1;
                    }
               
                    FRAME(nRowsLeft) = nRowsLeft - FRAME(nRowsGpu)+2;
                }
                FRAME(nGPU)=2;
                FRAME(gpuPos) = nRows-nRowsLeft;
                __sync_synchronize();
                SYNC(GpuKernelHybridWithStreams);
            }

#ifdef CUDA_DARTS_DEBUG
			size_t gpu_mem_total_t = 0;
			size_t gpu_mem_avail_t = 0;
			size_t gpu_mem_valid_t = 0;
			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
			gpu_mem_valid_t = gpu_mem_avail_t -XMB;

			std::cout<<"CpuSync: gpu memory total: "<<gpu_mem_total_t<<std::endl;
			std::cout<<"CpuSync: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
			std::cout<<"CpuSync: reqire gpu memory : "<<sizeof(double)*nCols*FRAME(nRowsGpu)<<std::endl;
#endif		


#ifdef CUDA_DARTS_DEBUG
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: gpuPos: "<<FRAME(gpuPos)<<std::endl;
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: nRowsGpu: "<<FRAME(nRowsGpu)<<std::endl;
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: reset rowsLeft: "<<FRAME(nRowsLeft)<<std::endl;
				std::cout<<"CpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;

#endif

		}
		if(FRAME(nRowsLeft)==0){
                SYNC(Swap);
		}else{

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"CpuSync: invoke new CpuLoop!"<<std::endl;
#endif
	
            double rt = FRAME(cmCpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            pthread_mutex_lock(&mutex);
            nRowsLeft = FRAME(nRowsLeft);
            if(nRowsLeft==0){
                SYNC(Swap);
            }else{
                
                if(nRowsLeft <= FRAME(lastCnt)){
                    FRAME(nRowsCpu) = nRowsLeft;
                    FRAME(nRowsLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        t1 = (1-cpuStepR)*nRowsCpu;
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)){
                        t1 = nRowsCpu;
                    }else{
                        t1 = (1+cpuStepR)*nRowsCpu;
                    }

                    if(nRowsLeft<t1){
                        FRAME(nRowsCpu) = nRowsLeft*rt;
                    }else{
                        FRAME(nRowsCpu) = t1;
                    }
              


                    FRAME(nRowsLeft) = nRowsLeft - FRAME(nRowsCpu)+2;
                }

                arrnCpuEdge[0] = FRAME(nCols);
                arrnCpuEdge[2] = FRAME(nRowsCpu);
                //============================= invoke new cpu codelet ==============================//
                chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
                chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
                calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
		        int nCPU = arrnCpuGrid[0]*arrnCpuGrid[1];
                FRAME(nCPU) = nCPU;
                SETSYNC(CpuSync,nCPU,nCPU);

				FRAME(cpuPos) = nRows-nRowsLeft;
                arrnCpuPos[0]=0;
                arrnCpuPos[1] = FRAME(cpuPos) ;
                __sync_synchronize();
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"CpuSync: invoked  new CPU, cpuPos: "<<FRAME(cpuPos)<<std::endl;

            std::cout<<"CpuSync: invoked  new CPU, nRowsCpu: "<<FRAME(nRowsCpu)<<std::endl;

            std::cout<<"CpuSync: invoked  new CPU, nRowsLeft: "<<FRAME(nRowsLeft)<<std::endl;

#endif
			__sync_bool_compare_and_swap(&FRAME(CpuFinish),true,false);
			for(size_t i =0; i<FRAME(nCPU);++i){
				SYNC(CpuLoop[i]);
			}

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync new CpuLoop: cpuPos: "<<FRAME(cpuPos)<<std::endl;
			std::cout<<"cpuSync new CpuLoop: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
#endif
		
		}

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
#endif
	}
	EXIT_TP();

}


extern "C"
void 
Stencil2D4ptSwapCD::fire(void) 
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke swap"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(Swap);

	double GpuRatio = FRAME(GpuRatio);

	uint32_t	nGPU = FRAME(nGPU);
    int nCPU;
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Swap: ts: "<<FRAME(ts)<<std::endl;
#endif	


	size_t ts = --FRAME(ts);
	if(GpuRatio == 0.0){

		SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
		if(ts!=0){
            nCPU = FRAME(nCPU);
			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop[i]);
			}
		}else{
			SYNC(sync);
		}
	}else{

		double *d_dst = FRAME(d_dst);
		double *dst = FRAME(New);
		double *src = FRAME(Initial);

		uint64_t	gpuPos = FRAME(gpuPos);
		uint64_t	cpuPos = FRAME(cpuPos);
		
		uint64_t nRows   = FRAME(nRows);
		uint64_t nCols   = FRAME(nCols);
		uint64_t nRowsGpu = FRAME(nRowsGpu);
		uint64_t nRowsCpu = FRAME(nRowsCpu);
		
		
		cudaError err0;
		err0 = cudaDeviceSynchronize();
                
#ifdef CUDA_ERROR_CHECKING
		if(err0!=cudaSuccess){
			std::cout<<"swap: cuda deviceSynchronize: "<<cudaGetErrorString(err0)<<std::endl;
			exit(-1);
		}
#endif
		
        if(ts!=0){
                    
            if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                FRAME(nRowsGpuInit) = FRAME(nRowsGpuInit)*(1-FRAME(gpuStepR));
                FRAME(nRowsLeftInit) = FRAME(nRows) - FRAME(nRowsGpuInit)-FRAME(nRowsCpuInit) +4;
                FRAME(cpuPosInit) = FRAME(nRowsGpuInit)-2;
            } 
            FRAME(gpuPos)	=FRAME(gpuPosInit);
		    FRAME(cpuPos)	=FRAME(cpuPosInit);
		    FRAME(nRowsGpu)	=FRAME(nRowsGpuInit);
		    FRAME(nRowsCpu) =FRAME(nRowsCpuInit);
		    FRAME(nRowsLeft)=FRAME(nRowsLeftInit);
			FRAME(CpuFinish) = false;
			FRAME(GpuFinish) = false;
			FRAME(nCPU) = FRAME(nCPUInit);
			FRAME(nGPU) = FRAME(nGPUInit);
            FRAME(gpuCnt) = 0;
            FRAME(cpuCnt) = 0;
#ifdef CUDA_DARTS_DEBUG
	        std::cout<<"swap: reset gpuPos: "<<FRAME(gpuPos)<<std::endl;
		    std::cout<<"swap: reset cpuPos: "<<FRAME(cpuPos)<<std::endl;
		    std::cout<<"swap: reset gpu rows: "<<FRAME(nRowsGpu)<<std::endl;
		    std::cout<<"swap: reset cpu rows: "<<FRAME(nRowsCpu)<<std::endl;
		    std::cout<<"swap: reset cpu rowsLeft: "<<FRAME(nRowsLeft)<<std::endl;

		    std::cout<<"swap: reset CpuFinsh: "<<FRAME(CpuFinish)<<std::endl;

		    std::cout<<"swap: reset GpuFinsh: "<<FRAME(GpuFinish)<<std::endl;

#endif

            int *arrnEdge               = FRAME(arrnEdge);
            int *arrnCpuPos             = FRAME(arrnCpuPos);
            int *arrnCpuEdge            = FRAME(arrnCpuEdge);
            int *arrnCpuTile            = FRAME(arrnCpuTile); 
            int *arrnCpuGridTileBase    = FRAME(arrnCpuGridTileBase);
            int *arrnCpuBlock           = FRAME(arrnCpuBlock);
            int *arrnCpuGrid            = FRAME(arrnCpuGrid);

            arrnCpuEdge[0] = FRAME(nCols);
            arrnCpuEdge[2] = FRAME(nRowsCpu);

            arrnCpuPos[0]=0;
            arrnCpuPos[1] = FRAME(cpuPos) ;

            chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
            chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
            calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
		    nCPU = arrnCpuGrid[0]*arrnCpuGrid[1];
            FRAME(nCPU) = nCPU;
            SETSYNC(CpuSync,nCPU,nCPU);

		    SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );

            SYNC(GpuKernelHybridWithStreams);
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
	std::cout<<"invoke Sync!"<<std::endl;
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





