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
	uint64_t gpuWL = FRAME(gpuWL);

	d_size = sizeof(double)*gpuWL*nCols;
	
	int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*gpuWL/tile_y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = gpuPos*nCols;

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuKernelWithAllTimeSteps: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: gpuWL:"<<gpuWL<<std::endl;

	std::cout<<"GpuKernelWithAllTimeSteps: blockDimx:"<<blockDimx<<std::endl;
	
	std::cout<<"GpuKernelWithAllTimeSteps: gridDimx="<<gridDimx<<",gridDimy="<<gridDimy<<std::endl;
#endif
	
	int64_t d_size_sharedCols = sizeof(double)*gpuWL*gridDimx*2;
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
	int gridDimy_rows = std::ceil(1.0*gpuWL/tile_y);


	int blockDimx_cols = 1 ;
	int blockDimy_cols = (gpuWL>NUM_THREADS)?NUM_THREADS:nRows;
	int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*gpuWL/blockDimy_cols);

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
		gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,gpuWL, nCols);
		gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,gpuWL, nCols);
		gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,gpuWL,nCols);
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
	size_t gpuWL   = FRAME(gpuWL);

	double *d_dst ;

	double *d_sharedCols ;
	double *d_sharedRows ;
	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);

	double d_size;	

	d_size = sizeof(double)*gpuWL*nCols;
	
	int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*gpuWL/tile_y); //GRID_TILE_Y=10, it needs to change.
	dim3 dimGrid(gridDimx,gridDimy);
	dim3 dimBlock(blockDimx,blockDimy);
	
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = gpuPos*nCols;

	int64_t d_size_sharedCols = sizeof(double)*gpuWL*gridDimx*2;
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
	std::cout<<"GpuKernel: gpuWL:"<<gpuWL<<std::endl;

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
    	int gridDimy_rows = std::ceil(1.0*gpuWL/tile_y);
    
    
    	int blockDimx_cols = 1 ;
    	int blockDimy_cols = (gpuWL>NUM_THREADS)?NUM_THREADS:nRows;
    	int gridDimx_cols = gridDimx;
    	int gridDimy_cols = std::ceil(1.0*gpuWL/blockDimy_cols);
    
    
    	dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
    	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    
    	dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
    	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
    
    
    	gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,gpuWL, nCols);
    	gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,gpuWL, nCols);
    	gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,gpuWL,nCols);
	
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
        //int gridDimy_rows = std::ceil(1.0*gpuWL/tile_y);
    
    
    	int blockDimx_cols = 1 ;
    	int blockDimy_cols;
        //int blockDimy_cols = (gpuWL>NUM_THREADS)?NUM_THREADS:gpuWL;
    	int gridDimx_cols = gridDimx;
        int gridDimy_cols; 
    	//int gridDimy_cols = std::ceil(1.0*gpuWL/blockDimy_cols);
    	
        //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
    	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
    
    	//dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
    	//dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
    
        uint64_t d_size_stream; 
        uint64_t nRowsStream;
        int nTile_y = (gpuWL/(tile_y*nStream));
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
            nRowsStream = (i==(nStream-1))?(gpuWL-i*chunk):(chunk+2);
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
            gpu_kernel5_stream_cp_rows(FRAME(stream[i]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+i*nTile_y*2*nCols, tile_y,gpuWL, nCols);
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
    
    int gpuWLBlock = nStream*chunk2 + nRows-nGPU*nStream*chunk;
    int64_t gpuWLStream;
	int64_t d_size_stream;
   
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*gpuWLBlock/tile_y); 
    
	dim3 dimBlock(blockDimx,blockDimy);
//	dim3 dimGrid(gridDimx,gridDimy);
    
    int64_t d_size = sizeof(double)*gpuWLBlock*nCols;
    int64_t d_size_sharedCols = sizeof(double) * gpuWLBlock*gridDimx*2 ;
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
	std::cout<<"GpuKernelPureGpuWithStreams: gpuWLBlock:"<<gpuWLBlock<<std::endl;
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
                gpuWLStream = ((i==(nGPU-1))&&(j==(nStream-1)))? (nRows-ps*chunk) :chunk2;
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

                d_size_stream = sizeof(double)*nCols*gpuWLStream;
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"GpuKernelPureGpuWithStreams multiple streams: stream["<<j<<"]"<<",ps:"<<ps<<", size:"<<d_size_stream<<",d_pos:"<<d_pos/nCols<<",h_pos:"<<h_pos/nCols<<",gpuWLStream:"<<gpuWLStream<<std::endl;
#endif
                err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
                if(err1!=cudaSuccess){
                    std::cout<<"GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                    exit(-1);
                }
#endif
            
                gridDimy_rows = std::ceil(1.0*gpuWLStream/tile_y);
                dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
                gpu_kernel5_stream_cp_rows(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+j*nTile_y*2*nCols, tile_y,gpuWLStream, nCols);
            
#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp rows: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
                blockDimy_cols = (gpuWLStream>NUM_THREADS)?NUM_THREADS:gpuWLStream;
                gridDimy_cols = std::ceil(1.0*gpuWLStream/blockDimy_cols);
        	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
                dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
                int addrCol = j*chunk2*2*gridDimx_cols; 
                gpu_kernel5_stream_cp_cols(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+addrCol, d_sharedRows, tile_x,tile_x,gpuWLStream, nCols);

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp cols: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
	            int gridDimy_stream = std::ceil(1.0*gpuWLStream/tile_y);
	            dim3 dimGrid_stream(gridDimx,gridDimy_stream);
                gpu_kernel5_stream(FRAME(stream[j]) ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+addrCol,d_sharedRows+j*nTile_y*2*nCols,tile_y,gpuWLStream,nCols);
            
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
    uint64_t gpuWL = FRAME(gpuWL);
    int nStream = FRAME(nStream);
    int vnStream = nStream*nGPU;
    
    int tile_y = GRID_TILE_Y;
	int tile_x = NUM_THREADS;
    
    int nTile_y = gpuWL/(tile_y * vnStream);

    int chunk = nTile_y*tile_y;
    int chunk2= chunk+2;
    
    int gpuWLBlock = nStream*chunk2 + gpuWL-nGPU*nStream*chunk;
    int64_t gpuWLStream;
	int64_t d_size_stream;
   
	int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
	int blockDimy = 1;
	int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	int gridDimy = std::ceil(1.0*gpuWLBlock/tile_y); 
    
	dim3 dimBlock(blockDimx,blockDimy);
//	dim3 dimGrid(gridDimx,gridDimy);
    
    int64_t d_size = sizeof(double)*(gpuWLBlock)*nCols;
    int64_t d_size_sharedCols = sizeof(double) * (gpuWLBlock+3)*gridDimx*2 ;
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
	std::cout<<"GpuKernelHybridWithStreams: gpuWLBlock:"<<gpuWLBlock<<std::endl;
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
            gpuWLStream = ((i==(nGPU-1))&&(j==(nStream-1)))? (gpuWL-ps*chunk) :chunk2;
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

            d_size_stream = sizeof(double)*nCols*gpuWLStream;
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams multiple streams: stream["<<j<<"]"<<",ps:"<<ps<<", size:"<<d_size_stream<<",d_pos:"<<d_pos/nCols<<",h_pos:"<<h_pos/nCols<<",gpuWLStream:"<<gpuWLStream<<std::endl;
#endif
            err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
            if(err1!=cudaSuccess){
                std::cout<<"GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                exit(-1);
            }
#endif
        
            gridDimy_rows = std::ceil(1.0*gpuWLStream/tile_y);
            dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
            gpu_kernel5_stream_cp_rows(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+j*nTile_y*2*nCols, tile_y,gpuWLStream, nCols);
        
#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp rows: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        
            blockDimy_cols = (gpuWLStream>NUM_THREADS)?NUM_THREADS:gpuWLStream;
            gridDimy_cols = std::ceil(1.0*gpuWLStream/blockDimy_cols);
    	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
            dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
            int addrCol = j*chunk2*2*gridDimx_cols; 
            gpu_kernel5_stream_cp_cols(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+addrCol, d_sharedRows, tile_x,tile_x,gpuWLStream, nCols);

#ifdef CUDA_ERROR_CHECKING
            err3 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda cp cols: "<<cudaGetErrorString(err3)<<std::endl;
                exit(-1);
            }
#endif
        
	        int gridDimy_stream = std::ceil(1.0*gpuWLStream/tile_y);
	        dim3 dimGrid_stream(gridDimx,gridDimy_stream);
            gpu_kernel5_stream(FRAME(stream[j]) ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+addrCol,d_sharedRows+j*nTile_y*2*nCols,tile_y,gpuWLStream,nCols);
        
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
    uint64_t wlLeft = FRAME(wlLeft);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: GPU finish! and wlLeft: "<<wlLeft<<std::endl;
#endif
    ++FRAME(gpuCnt); 
    if(wlLeft == 0){
        SYNC(Swap);
    }else{

        if(__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,false)&&(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false) )){ 
            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));
            double gpuStepR = FRAME(gpuStepR);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: invoke new GPU"<<std::endl;
#endif
            
            pthread_mutex_lock(&mutex);
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0){
                SYNC(Swap);
            }else{
             
                if(wlLeft <= FRAME(lastCnt)){
                    FRAME(gpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*gpuWL;
                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = gpuWL;
                    }else{
                        t1 = (1+gpuStepR)*gpuWL;
                    }

                    if(wlLeft<t1){
                        FRAME(gpuWL) = wlLeft*rt;
                    }else{
                        FRAME(gpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
                }
                
                FRAME(nGPU)=2;
                FRAME(gpuPos) = nRows-wlLeft ;
                __sync_synchronize();

                SYNC(GpuKernelHybridWithStreams);
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, gpuPos: "<<FRAME(gpuPos)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, gpuWL: "<<FRAME(gpuWL)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams: invoked a new GPU, wlLeft: "<<FRAME(wlLeft)<<std::endl;

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
	uint64_t gpuWL = FRAME(gpuWL);

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

	gpu_mem_valid_t = gpu_mem_avail_t -XMB;
	
	std::cout<<"GpuLoop: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
#endif
	int64_t chunk=gpuWL/nGPU;
	uint64_t nRows_bk = (Id==(nGPU-1))?(gpuWL-Id*chunk):(chunk+2);
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
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
	uint64_t	cpuWL = FRAME(cpuWL);	
	uint64_t	nCPU = FRAME(nCPU);
	uint64_t	chunk = (cpuWL )/nCPU;
	uint64_t	cpuPos = FRAME(cpuPos);
	size_t		pos1 = (cpuPos + chunk*Id)*nCols;
	uint64_t	nRows_bk = (Id==(nCPU-1))? (cpuWL-chunk*Id-1):(chunk +1) ;
	double *src = h_src+pos1;
	double *dst = h_dst+pos1;

	double *d_dst = FRAME(d_dst);
	uint64_t wlLeft = FRAME(wlLeft);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuLoop["<<Id<<"]: cpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: cpuWL:"<<cpuWL<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: chunk:"<<chunk<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nRows_bk:"<<nRows_bk<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: wlLeft:"<<wlLeft<<std::endl;
#endif

	computeInner_stencil2d_v2(dst,src,nRows_bk,nCols);

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
	uint64_t wlLeft = FRAME(wlLeft);
	//uint64_t nRows = FRAME(nRows);
	RESET(CpuSync);

    ++FRAME(cpuCnt);

    if(wlLeft==0.0){
		SYNC(Swap);
	
	}else{
		__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,true);
	
        size_t	gpuPos = FRAME(gpuPos);
		size_t	cpuPos = FRAME(cpuPos);
	
		uint64_t gpuWL = FRAME(gpuWL);
		uint64_t cpuWL = FRAME(cpuWL);
		uint64_t nCols = FRAME(nCols);
		uint64_t nRows = FRAME(nRows);

		double *dst = FRAME(New);
		double *d_dst = FRAME(d_dst);
		double cmCpu = FRAME(cmCpu);
		double cmGpu = FRAME(cmGpu);
        
        double cpuStepR = FRAME(cpuStepR);
        double gpuStepR = FRAME(gpuStepR);
        
        if(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false)) {

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync: gpu kernel finish and invode a new gpu kernel"<<std::endl;
#endif

            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0.0){
                SYNC(Swap);
            }else{
                
                if(wlLeft <= FRAME(lastCnt)){
                    FRAME(gpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{
                    uint64_t t1;
                  
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*gpuWL;
                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = gpuWL;
                    }else{
                        t1 = (1+gpuStepR)*gpuWL;
                    }
                    
                    if(wlLeft<t1){
                        FRAME(gpuWL) = wlLeft*rt;
                    }else{
                        FRAME(gpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
                }
                FRAME(nGPU)=2;
                FRAME(gpuPos) = nRows-wlLeft;
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
			std::cout<<"CpuSync: reqire gpu memory : "<<sizeof(double)*nCols*FRAME(gpuWL)<<std::endl;
#endif		


#ifdef CUDA_DARTS_DEBUG
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: gpuPos: "<<FRAME(gpuPos)<<std::endl;
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: gpuWL: "<<FRAME(gpuWL)<<std::endl;
				std::cout<<"CpuSync new GpuKernelHybridWithStreams: reset rowsLeft: "<<FRAME(wlLeft)<<std::endl;
				std::cout<<"CpuSync: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;

#endif

		}
		if(FRAME(wlLeft)==0){
                SYNC(Swap);
		}else{

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"CpuSync: invoke new CpuLoop!"<<std::endl;
#endif
	
            double rt = FRAME(cmCpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            pthread_mutex_lock(&mutex);
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0){
                SYNC(Swap);
            }else{
                
                if(wlLeft <= FRAME(lastCnt)){
                    FRAME(cpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        t1 = (1-cpuStepR)*cpuWL;
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)){
                        t1 = cpuWL;
                    }else{
                        t1 = (1+cpuStepR)*cpuWL;
                    }

                    if(wlLeft<t1){
                        FRAME(cpuWL) = wlLeft*rt;
                    }else{
                        FRAME(cpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(cpuWL)+2;
                }

				FRAME(cpuPos) = nRows-wlLeft;
                __sync_synchronize();
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"CpuSync: invoked  new CPU, cpuPos: "<<FRAME(cpuPos)<<std::endl;

            std::cout<<"CpuSync: invoked  new CPU, cpuWL: "<<FRAME(cpuWL)<<std::endl;

            std::cout<<"CpuSync: invoked  new CPU, wlLeft: "<<FRAME(wlLeft)<<std::endl;

#endif
			__sync_bool_compare_and_swap(&FRAME(CpuFinish),true,false);
			for(size_t i =0; i<FRAME(nCPU);++i){
				SYNC(CpuLoop[i]);
			}

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync new CpuLoop: cpuPos: "<<FRAME(cpuPos)<<std::endl;
			std::cout<<"cpuSync new CpuLoop: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
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

	uint32_t	nCPU = FRAME(nCPU);
	uint32_t	nGPU = FRAME(nGPU);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Swap: ts: "<<FRAME(ts)<<std::endl;
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
	}else{

		double *d_dst = FRAME(d_dst);
		double *dst = FRAME(New);
		double *src = FRAME(Initial);

		uint64_t	gpuPos = FRAME(gpuPos);
		uint64_t	cpuPos = FRAME(cpuPos);
		
		uint64_t nRows   = FRAME(nRows);
		uint64_t nCols   = FRAME(nCols);
		uint64_t gpuWL = FRAME(gpuWL);
		uint64_t cpuWL = FRAME(cpuWL);
		
		
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
                FRAME(gpuWLInit) = FRAME(gpuWLInit)*(1-FRAME(gpuStepR));
                FRAME(wlLeftInit) = FRAME(nRows) - FRAME(gpuWLInit)-FRAME(cpuWLInit) +4;
                FRAME(cpuPosInit) = FRAME(gpuWLInit)-2;
            } 
            FRAME(gpuPos)	=FRAME(gpuPosInit);
		    FRAME(cpuPos)	=FRAME(cpuPosInit);
		    FRAME(gpuWL)	=FRAME(gpuWLInit);
		    FRAME(cpuWL) =FRAME(cpuWLInit);
		    FRAME(wlLeft)=FRAME(wlLeftInit);
			FRAME(CpuFinish) = false;
			FRAME(GpuFinish) = false;
			FRAME(nCPU) = FRAME(nCPUInit);
			FRAME(nGPU) = FRAME(nGPUInit);
            FRAME(gpuCnt) = 0;
            FRAME(cpuCnt) = 0;
#ifdef CUDA_DARTS_DEBUG
	        std::cout<<"swap: reset gpuPos: "<<FRAME(gpuPos)<<std::endl;
		    std::cout<<"swap: reset cpuPos: "<<FRAME(cpuPos)<<std::endl;
		    std::cout<<"swap: reset gpu rows: "<<FRAME(gpuWL)<<std::endl;
		    std::cout<<"swap: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
		    std::cout<<"swap: reset cpu rowsLeft: "<<FRAME(wlLeft)<<std::endl;

		    std::cout<<"swap: reset CpuFinsh: "<<FRAME(CpuFinish)<<std::endl;

		    std::cout<<"swap: reset GpuFinsh: "<<FRAME(GpuFinish)<<std::endl;

#endif

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



extern "C"
void 
Stencil3D7ptCpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint64_t Id = getID();
	RESET(CpuLoop37[Id]);	
	
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke CPU["<<Id<<"]"<<std::endl;	
#endif
	double	*h_src  = FRAME(Initial);
	double	*h_dst      = FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);
	uint64_t	cpuWL = FRAME(cpuWL);	
	uint64_t	nCPU = FRAME(nCPU);
	uint64_t	chunk = (cpuWL )/nCPU;
	uint64_t	cpuPos = FRAME(cpuPos);
	size_t		pos1 = (cpuPos + chunk*Id)*nCols*nRows;
	uint64_t	nSlicesChunk = (Id==(nCPU-1))? (cpuWL-chunk*Id-1):(chunk +1) ;
	double *src = h_src+pos1;
	double *dst = h_dst+pos1;

	double *d_dst = FRAME(d_dst);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: chunk:"<<chunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nSlicesChunk:"<<nSlicesChunk<<std::endl;
#endif

	computeInner_stencil37(dst,src,nRows,nCols,nSlicesChunk);

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
		SYNC(Swap37);
	}else{
		SYNC(CpuSync37);
	}

//	std::cout<<"cpu["<<Id<<"]: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
	EXIT_TP();
}


extern "C"
void 
Stencil3D7ptSwapCD::fire(void) 
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke swap37"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(Swap37);

	double GpuRatio = FRAME(GpuRatio);

	uint32_t	nCPU = FRAME(nCPU);
	uint32_t	nGPU = FRAME(nGPU);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Swap: ts: "<<FRAME(ts)<<std::endl;
#endif	


	size_t ts = --FRAME(ts);
	if(GpuRatio == 0.0){

		SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
		if(ts!=0){

			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop37[i]);
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
		uint64_t nSlices = FRAME(nSlices);
        uint64_t gpuWL = FRAME(gpuWL);
		uint64_t cpuWL = FRAME(cpuWL);
    
        if(ts!=0){
                    
            if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                FRAME(gpuWLInit) = FRAME(gpuWLInit)*(1-FRAME(gpuStepR));
                FRAME(wlLeftInit) = FRAME(tWL) - FRAME(gpuWLInit)-FRAME(cpuWLInit) +4;
                FRAME(cpuPosInit) = FRAME(gpuWLInit)-2;
            } 
            FRAME(gpuPos)	=FRAME(gpuPosInit);
		    FRAME(cpuPos)	=FRAME(cpuPosInit);
		    FRAME(gpuWL)	=FRAME(gpuWLInit);
		    FRAME(cpuWL) =FRAME(cpuWLInit);
		    FRAME(wlLeft)=FRAME(wlLeftInit);
			FRAME(CpuFinish) = false;
			FRAME(GpuFinish) = false;
			FRAME(nCPU) = FRAME(nCPUInit);
			FRAME(nGPU) = FRAME(nGPUInit);
            FRAME(gpuCnt) = 0;
            FRAME(cpuCnt) = 0;
#ifdef CUDA_DARTS_DEBUG
	        std::cout<<"swap: reset gpuPos: "<<FRAME(gpuPos)<<std::endl;
		    std::cout<<"swap: reset cpuPos: "<<FRAME(cpuPos)<<std::endl;
		    std::cout<<"swap: reset gpu rows: "<<FRAME(gpuWL)<<std::endl;
		    std::cout<<"swap: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
		    std::cout<<"swap: reset cpu rowsLeft: "<<FRAME(wlLeft)<<std::endl;

		    std::cout<<"swap: reset CpuFinsh: "<<FRAME(CpuFinish)<<std::endl;

		    std::cout<<"swap: reset GpuFinsh: "<<FRAME(GpuFinish)<<std::endl;

#endif

		    SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );

            SYNC(GpuKernelHybridWithStreams37);
		    for (size_t i =0;i<nCPU;++i){
		    	SYNC(CpuLoop37[i]);
		    }

		}else{
            SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
            SYNC(sync);
		}
    
    
    }


	EXIT_TP();
}


extern "C"
void
Stencil3D7ptGpuKernelWithAllTimeStepsCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelWithAllTimeSteps37!"<<std::endl;	
#endif
	LOAD_FRAME(StencilTP);
	RESET(GpuKernelWithAllTimeSteps37);

	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double *d_sharedSlices;

	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);
	uint64_t gpuWL = FRAME(gpuWL);

	double d_size = FRAME(d_size);

    int64_t d_size_sharedCols = FRAME(d_size_sharedCols);
	int64_t d_size_sharedRows = FRAME(d_size_sharedRows);
	int64_t d_size_sharedSlices = FRAME(d_size_sharedSlices);
	
    cudaError err1,err2,err3,err4;
	err1 = cudaMalloc( (void **) &d_dst, d_size);
	err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
	err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
	err4 = cudaMalloc( (void **) &d_sharedSlices, d_size_sharedSlices);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}
	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc3: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda malloc4: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

	FRAME(d_dst) = d_dst;
		
	uint64_t gpuPos = FRAME(gpuPos);
	size_t	pos1 = gpuPos*nRows*nCols;
	err1 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif

	int tile_x = FRAME(tile_x); //16
	int tile_y = FRAME(tile_y); //16
    int tile_z = FRAME(tile_z); //100

    int blockDimx = FRAME(blockDimx);   // tile_x
    int blockDimy = FRAME(blockDimy);   // tile_y
    int blockDimz = FRAME(blockDimz);   // 1

	int gridDimx = FRAME(gridDimx);     //  x/tile_x 
	int gridDimy = FRAME(gridDimy);     //  y/tile_y
    int gridDimz = FRAME(gridDimz);     //  z/tile_z
    
    dim3 dimGrid(gridDimx,gridDimy,gridDimz);
	dim3 dimBlock(blockDimx,blockDimy,blockDimz);


#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernelWithAllTimeSteps: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nRows:"<<nRows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nCols:"<<nCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nSlices:"<<nSlices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_x:"<<tile_x<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_y:"<<tile_y<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_z:"<<tile_z<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimx:"<<blockDimx<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimy:"<<blockDimy<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimz:"<<blockDimz<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimx="<<gridDimx<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimy="<<gridDimy<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimz="<<gridDimz<<std::endl;

#endif
	

    
    int numThreads= tile_x*tile_y;

    int blockDimx_slices = numThreads; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/numThreads);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz;

    int blockDimx_rows = numThreads;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/numThreads);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = numThreads;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/numThreads);
    int gridDimz_cols = gridDimz;

#ifdef CUDA_DARTS_DEBUG

	std::cout<<"GpuKernelWithAllTimeSteps: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimy_rows:"<<gridDimy_rows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimz_rows:"<<gridDimz_rows<<std::endl;
	
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimy_cols:"<<gridDimy_cols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimz_cols:"<<gridDimz_cols<<std::endl;

#endif
	
	dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
	dim3 dimBlock_slices(blockDimx_slices,blockDimy_slices,blockDimz_slices);
	
    dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows,blockDimz_rows);

	dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);
	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols,blockDimz_cols);


	size_t ts = FRAME(ts);
	while(ts-- >0){
        gpu_kernel37_cp_slices(dimGrid_slices,dimBlock_slices,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
        gpu_kernel37_cp_rows(dimGrid_rows,dimBlock_rows,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
        gpu_kernel37_cp_cols(dimGrid_cols,dimBlock_cols,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
        gpu_kernel37(dimGrid,dimBlock,d_dst,d_sharedRows,d_sharedCols,d_sharedSlices,nRows,nCols,nSlices,tile_x,tile_y,tile_z);
	}
	
	err1 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda deviceSynchronize: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif
        
#ifdef VERIFICATION
        if(FRAME(ts)%2==0){
    	    SWAP_PTR(&h_dst ,&h_src);
        }
#endif
	//copy from GPU  to CPU
        err1=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif

#ifdef CUDA_DARTS_DEBUG


#endif

	err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedCols);
	err3 = cudaFree(d_sharedRows);
    err4 = cudaFree(d_sharedSlices);
	
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}

	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

	SYNC(sync);

	EXIT_TP();

}


extern "C"
void
Stencil3D7ptGpuKernelPureGpuWithStreamsCD::fire(void)
{

	LOAD_FRAME(StencilTP);
    RESET(GpuKernelPureGpuWithStreams37);
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelPureGpuWithStreams37"<<std::endl;
#endif

	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double *d_sharedSlices;

	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);
	uint64_t gpuWL = FRAME(gpuWL);

	int tile_x = FRAME(tile_x); //16
	int tile_y = FRAME(tile_y); //16
    int tile_z = FRAME(tile_z); //100
    
    uint32_t nGPU = FRAME(nGPU);
    uint64_t gpuPos = FRAME(gpuPos); 
    int nStream = FRAME(nStream);
    int vnStream = nStream*nGPU;

    int nTile_z = std::ceil(1.0*gpuWL/(tile_z*vnStream));
    int chunk = nTile_z*tile_z;
    int nSlicesChunk = nTile_z*tile_z + 2;
    int nSlicesChunkInit = nSlicesChunk;
    
    int blockDimx = FRAME(blockDimx);   // tile_x
    int blockDimy = FRAME(blockDimy);   // tile_y
    int blockDimz = FRAME(blockDimz);   // 1

	int gridDimx = FRAME(gridDimx);     //  x/tile_x 
	int gridDimy = FRAME(gridDimy);     //  y/tile_y
    int gridDimz = nTile_z;             //  chunk/tile_z

    int gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);

	double d_size = sizeof(double)*nRows*nCols*nSlicesChunk;
    
	double d_sizeInit = d_size;
	int64_t d_size_sharedRows   = sizeof(double)*nCols*gridDimy*2*nSlicesChunk;
    int64_t d_size_sharedCols   = sizeof(double)*nRows*gridDimx*2*nSlicesChunk;
	int64_t d_size_sharedSlices = sizeof(double)*nRows*nCols*2*gridDimz2;

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37 : require memory size:"<<(d_size + d_size_sharedCols + d_size_sharedRows)*4<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: Gpupos:"<<FRAME(gpuPos)<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: nSlicesChunk:"<<nSlicesChunk<<std::endl;
#endif

    cudaError err1,err2,err3,err4;
	err1 = cudaMalloc( (void **) &d_dst, d_size*nStream);
	err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols*nStream);
	err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows*nStream);
	err4 = cudaMalloc( (void **) &d_sharedSlices, d_size_sharedSlices*nStream);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}
	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda malloc3: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda malloc4: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

	FRAME(d_dst) = d_dst;
    
    //dim3 dimGrid(gridDimx,gridDimy,gridDimz);
	dim3 dimBlock(blockDimx,blockDimy,blockDimz);

    
    int numThreads= tile_x*tile_y;

    int blockDimx_slices = numThreads; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/numThreads);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz2;

    int blockDimx_rows = numThreads;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/numThreads);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz2;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = numThreads;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/numThreads);
    int gridDimz_cols = gridDimz2;


#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"GpuKernelPureGpuWithStreams3737: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: nRows:"<<nRows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: nCols:"<<nCols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: nSlices:"<<nSlices<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: tile_x:"<<tile_x<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: tile_y:"<<tile_y<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: tile_z:"<<tile_z<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx:"<<blockDimx<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy:"<<blockDimy<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimz:"<<blockDimz<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: gridDimx="<<gridDimx<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: gridDimy="<<gridDimy<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: gridDimz="<<gridDimz<<std::endl;

	std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_slices:"<<blockDimx_slices<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_slices:"<<blockDimy_slices<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_slices:"<<gridDimx_slices<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_slices:"<<gridDimy_slices<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_slices:"<<gridDimz_rows<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_rows:"<<gridDimy_rows<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_rows:"<<gridDimz_rows<<std::endl;
    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_cols:"<<gridDimy_cols<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_cols:"<<gridDimz_cols<<std::endl;

#endif
	
	//dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
	dim3 dimBlock_slices(blockDimx_slices,blockDimy_slices,blockDimz_slices);
	
    //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows,blockDimz_rows);

	//dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);
	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols,blockDimz_cols);


	size_t ts = FRAME(ts);
    uint64_t pos0 = gpuPos*nRows*nCols;
	while(ts-- >0){
	    nSlicesChunk =nSlicesChunkInit ;

        d_size = d_sizeInit;
        gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
        gridDimz_slices = gridDimz2;
        gridDimz_rows = gridDimz2;
        gridDimz_cols = gridDimz2;
        gridDimz = gridDimz2;

        for(size_t i=0;i<nGPU;++i){
            for(size_t j=0;j<nStream;++j){

#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelPureGpuWithStreams37 multiple streams: stream #"<<j<<"-1, error: "<<cudaGetErrorString(err1)<<std::endl;
                    exit(-1);
                }
#endif  
                size_t ps = i*nStream+j;
                
                int64_t h_pos = pos0 + ps*chunk*nRows*nCols;
                int64_t d_pos = j*nSlicesChunk*nRows*nCols; 
                int64_t s_srows_pos=j*nCols*gridDimy*2*nSlicesChunk;
                int64_t s_scols_pos=j*nRows*gridDimx*2*nSlicesChunk;
                int64_t s_sslices_pos=j*nRows*nCols*2*nTile_z;
                int64_t nSlicesLeft = gpuWL-ps*chunk;
                
             

                if(nSlicesLeft>1){
                
                    if(nSlicesLeft<=chunk+2){
                        d_size = sizeof(double)*nSlicesLeft*nRows*nCols;
                        nSlicesChunk = nSlicesLeft;
                        gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
                        gridDimz_slices = gridDimz2;
                        gridDimz_rows = gridDimz2;
                        gridDimz_cols = gridDimz2;
                        gridDimz = gridDimz2;
                    }
	                dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
                    dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
                    dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);

                    dim3 dimGrid (gridDimx,gridDimy,gridDimz);
                

#ifdef CUDA_DARTS_DEBUG
	                std::cout<<"GpuKernelPureGpuWithStreams3737: gpuPos:"<<gpuPos<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: d_size:"<<d_size<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_slices:"<<blockDimx_slices<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_slices:"<<blockDimy_slices<<std::endl;
                    std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_slices:"<<gridDimx_slices<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_slices:"<<gridDimy_slices<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_slices:"<<gridDimz_rows<<std::endl;
                    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_rows:"<<gridDimy_rows<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_rows:"<<gridDimz_rows<<std::endl;
                    std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy_cols:"<<gridDimy_cols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz_cols:"<<gridDimz_cols<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimx:"<<gridDimx<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimy:"<<gridDimy<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: grimDimz:"<<gridDimz2<<std::endl;
                    std::cout<<"GpuKernelPureGpuWithStreams3737: nGpu: "<<i<<",stream:# "<<j<<", h_pos:"<<h_pos<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: nGpu: "<<i<<",stream:# "<<j<<", d_pos:"<<d_pos<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: nGpu: "<<i<<",stream:# "<<j<<", nSlicesChunk:"<<nSlicesChunk<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: nGpu: "<<i<<",stream:# "<<j<<", nSlicesLeft:"<<nSlicesLeft<<std::endl;
	                std::cout<<"GpuKernelPureGpuWithStreams3737: nGpu: "<<i<<",stream:# "<<j<<", chunk:"<<chunk<<std::endl;
#endif
                
                    err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
                    if(err1!=cudaSuccess){
                        std::cout<<"GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                        exit(-1);
                    }
#endif
                    gpu_kernel37_cp_slices_stream(FRAME(stream[j]),dimGrid_slices,dimBlock_slices,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_cp_rows_stream(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_cp_cols_stream(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_stream(FRAME(stream[j]),dimGrid,dimBlock,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);

#ifdef CUDA_ERROR_CHECKING
                    err3 = cudaGetLastError();
                    if(cudaSuccess != err3){
                        std::cout<<"GpuKernelWithStream multiple streams: kernel5 cuda computation: "<<cudaGetErrorString(err3)<<std::endl;
                        exit(-1);
                    }
#endif
            
	                err3=cudaMemcpyAsync(h_dst+h_pos+nRows*nCols, d_dst+d_pos+nRows*nCols,d_size-sizeof(double)*nRows*nCols*2, cudaMemcpyDeviceToHost,FRAME(stream[j]));

#ifdef CUDA_ERROR_CHECKING
                    err3 = cudaGetLastError();
                    if(cudaSuccess != err3){
                        std::cout<<"GpuKernelWithStream multiple streams: kernel5 Asyn Memory copy from device to host: "<<cudaGetErrorString(err3)<<std::endl;
                        exit(-1);
                    }
#endif

                }
            }
        }

	    err1 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
		    std::cout<<"GpuKernelPureGpuWithStreams37: cuda deviceSynchronize: "<<cudaGetErrorString(err1)<<std::endl;
		    exit(-1);
	    }
#endif
    	SWAP_PTR(&h_dst ,&h_src);
    }
	
	err1 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda deviceSynchronize: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif
        
#ifdef VERIFICATION
        if(FRAME(ts)%2==0){
    	    SWAP_PTR(&h_dst ,&h_src);
        }
#endif

#ifdef CUDA_DARTS_DEBUG


#endif

	err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedCols);
	err3 = cudaFree(d_sharedRows);
    err4 = cudaFree(d_sharedSlices);
	
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda memcpy free d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}

	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelPureGpuWithStreams37: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

	SYNC(sync);

	EXIT_TP();
}


extern "C"
void
Stencil3D7ptGpuKernelHybridWithStreamsCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"Invoke GpuKernelHybridWithStreams"<<std::endl;
#endif
    LOAD_FRAME(StencilTP);
    RESET(GpuKernelHybridWithStreams37);

	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double *d_sharedSlices;

	double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);
	uint64_t gpuWL = FRAME(gpuWL);

	int tile_x = FRAME(tile_x); //16
	int tile_y = FRAME(tile_y); //16
    int tile_z = FRAME(tile_z); //100
    
    uint32_t nGPU = FRAME(nGPU);
    uint64_t gpuPos = FRAME(gpuPos); 
    int nStream = FRAME(nStream);
    int vnStream = nStream*nGPU;

    int nTile_z = std::ceil(1.0*gpuWL/(tile_z*vnStream));
    int chunk = nTile_z*tile_z;
    int nSlicesChunk = nTile_z*tile_z + 2;
   
    int blockDimx = FRAME(blockDimx);   // tile_x
    int blockDimy = FRAME(blockDimy);   // tile_y
    int blockDimz = FRAME(blockDimz);   // 1

	int gridDimx = FRAME(gridDimx);     //  x/tile_x 
	int gridDimy = FRAME(gridDimy);     //  y/tile_y
    int gridDimz = nTile_z;             //  chunk/tile_z

    int gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);

	double d_size = sizeof(double)*nRows*nCols*nSlicesChunk;
    
	int64_t d_size_sharedRows   = sizeof(double)*nCols*gridDimy*2*nSlicesChunk;
    int64_t d_size_sharedCols   = sizeof(double)*nRows*gridDimx*2*nSlicesChunk;
	int64_t d_size_sharedSlices = sizeof(double)*nRows*nCols*2*gridDimz2;

#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37 : require memory size:"<<(d_size + d_size_sharedCols + d_size_sharedRows)*4<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: nGPU:"<<nGPU<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: Gpupos:"<<FRAME(gpuPos)<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: nRows:"<<nRows<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: nCols:"<<nCols<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: nSlices:"<<nSlices<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: gpuWL:"<<gpuWL<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: nTile_z:"<<nTile_z<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: nSlicesChunk:"<<nSlicesChunk<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: d_size:"<<d_size<<std::endl;
#endif

    cudaError err1,err2,err3,err4;
	err1 = cudaMalloc( (void **) &d_dst, d_size*nStream);
	err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols*nStream);
	err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows*nStream);
	err4 = cudaMalloc( (void **) &d_sharedSlices, d_size_sharedSlices*nStream);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37:  cuda malloc1: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda malloc2: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}
	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda malloc3: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda malloc4: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

	FRAME(d_dst) = d_dst;
    
    //dim3 dimGrid(gridDimx,gridDimy,gridDimz);
	dim3 dimBlock(blockDimx,blockDimy,blockDimz);

    
    int numThreads= tile_x*tile_y;

    int blockDimx_slices = numThreads; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/numThreads);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz2;

    int blockDimx_rows = numThreads;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/numThreads);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz2;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = numThreads;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/numThreads);
    int gridDimz_cols = gridDimz2;


#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"GpuKernelHybridWithStreams37: GpuPos:"<<pos1/nCols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: nRows:"<<nRows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: nCols:"<<nCols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: nSlices:"<<nSlices<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: tile_x:"<<tile_x<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: tile_y:"<<tile_y<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: tile_z:"<<tile_z<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: blockDimx:"<<blockDimx<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: blockDimy:"<<blockDimy<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: blockDimz:"<<blockDimz<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: gridDimx="<<gridDimx<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: gridDimy="<<gridDimy<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: gridDimz="<<gridDimz<<std::endl;

	std::cout<<"GpuKernelHybridWithStreams37: blockDimx_slices:"<<blockDimx_slices<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_slices:"<<blockDimy_slices<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: grimDimx_slices:"<<gridDimx_slices<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_slices:"<<gridDimy_slices<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_slices:"<<gridDimz_rows<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_rows:"<<gridDimy_rows<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_rows:"<<gridDimz_rows<<std::endl;
    std::cout<<"GpuKernelHybridWithStreams37: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_cols:"<<gridDimy_cols<<std::endl;
	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_cols:"<<gridDimz_cols<<std::endl;
#endif
	
	//dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
	dim3 dimBlock_slices(blockDimx_slices,blockDimy_slices,blockDimz_slices);
	
    //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
	dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows,blockDimz_rows);

	//dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);
	dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols,blockDimz_cols);
    size_t ts = FRAME(ts);
    uint64_t pos0 = gpuPos*nRows*nCols;
    for(size_t i=0;i<nGPU;++i){
        for(size_t j=0;j<nStream;++j){

#ifdef CUDA_ERROR_CHECKING
            err1 = cudaGetLastError();
            if(cudaSuccess != err3){
                std::cout<<"GpuKernelHybridWithStreams37 multiple streams: stream #"<<j<<"-1, error: "<<cudaGetErrorString(err1)<<std::endl;
                exit(-1);
            }
#endif  
            size_t ps = i*nStream+j;
            int64_t h_pos = pos0 + ps*chunk*nRows*nCols;
            int64_t d_pos = j*nSlicesChunk*nRows*nCols; 
            int64_t s_srows_pos=j*nCols*gridDimy*2*nSlicesChunk;
            int64_t s_scols_pos=j*nRows*gridDimx*2*nSlicesChunk;
            int64_t s_sslices_pos=j*nRows*nCols*2*nTile_z;
            int64_t nSlicesLeft = gpuWL-ps*chunk;
           
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"GpuKernelHybridWithStreams37: nGpu: "<<i<<",stream:# "<<j<<", h_pos:"<<h_pos/(nRows*nCols)<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: nGpu: "<<i<<",stream:# "<<j<<", d_pos:"<<d_pos/(nRows*nCols)<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: nGpu: "<<i<<",stream:# "<<j<<", nSlicesChunk:"<<nSlicesChunk<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: nGpu: "<<i<<",stream:# "<<j<<", nSlicesLeft:"<<nSlicesLeft<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: nGpu: "<<i<<",stream:# "<<j<<", chunk:"<<chunk<<std::endl;
#endif
            if(nSlicesLeft>1){
                if(nSlicesLeft<=chunk+2){
                    d_size = sizeof(double)*nSlicesLeft*nRows*nCols;
                    nSlicesChunk = nSlicesLeft;
                    gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
                    gridDimz_slices = gridDimz2;
                    gridDimz_rows = gridDimz2;
                    gridDimz_cols = gridDimz2;
                    gridDimz = gridDimz2;
                }
                dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
                dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
                dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);
                
                dim3 dimGrid (gridDimx,gridDimy,gridDimz);


#ifdef CUDA_DARTS_DEBUG
	            std::cout<<"GpuKernelHybridWithStreams37: blockDimx_slices:"<<blockDimx_slices<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: blockDimy_slices:"<<blockDimy_slices<<std::endl;
                std::cout<<"GpuKernelHybridWithStreams37: grimDimx_slices:"<<gridDimx_slices<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimy_slices:"<<gridDimy_slices<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimz_slices:"<<gridDimz_rows<<std::endl;
                std::cout<<"GpuKernelHybridWithStreams37: blockDimx_rows:"<<blockDimx_rows<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: blockDimy_rows:"<<blockDimy_rows<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimx_rows:"<<gridDimx_rows<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimy_rows:"<<gridDimy_rows<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimz_rows:"<<gridDimz_rows<<std::endl;
                std::cout<<"GpuKernelHybridWithStreams37: blockDimx_cols:"<<blockDimx_cols<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: blockDimy_cols:"<<blockDimy_cols<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimx_cols:"<<gridDimx_cols<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimy_cols:"<<gridDimy_cols<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimz_cols:"<<gridDimz_cols<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimx:"<<gridDimx<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimy:"<<gridDimy<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: grimDimz:"<<gridDimz2<<std::endl;
#endif
                
                err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size, cudaMemcpyHostToDevice,FRAME(stream[j]));
#ifdef CUDA_ERROR_CHECKING
                if(err1!=cudaSuccess){
                    std::cout<<"GpuKernelHybridWithStreams37 multiple streams: cuda MemcpyAsync from host to device: "<<cudaGetErrorString(err1)<<std::endl;
                    exit(-1);
                }
#endif
                gpu_kernel37_cp_slices_stream(FRAME(stream[j]),dimGrid_slices,dimBlock_slices,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                gpu_kernel37_cp_rows_stream(FRAME(stream[j]),dimGrid_rows,dimBlock_rows,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                gpu_kernel37_cp_cols_stream(FRAME(stream[j]),dimGrid_cols,dimBlock_cols,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                gpu_kernel37_stream(FRAME(stream[j]),dimGrid,dimBlock,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelHybridWithStreams37 multiple streams: kernel5 cuda computation: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif
            
                err3=cudaMemcpyAsync(h_dst+h_pos+nRows*nCols, d_dst+d_pos+nRows*nCols,d_size-sizeof(double)*nRows*nCols*2, cudaMemcpyDeviceToHost,FRAME(stream[j]));

#ifdef CUDA_ERROR_CHECKING
                err3 = cudaGetLastError();
                if(cudaSuccess != err3){
                    std::cout<<"GpuKernelHybridWithStreams37 multiple streams: kernel5 Asyn Memory copy from device to host: "<<cudaGetErrorString(err3)<<std::endl;
                    exit(-1);
                }
#endif

            }
        }
    }
	
	err1 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda deviceSynchronize: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif
        
	err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedCols);
	err3 = cudaFree(d_sharedRows);
    err4 = cudaFree(d_sharedSlices);
	
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda memcpy free d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}

	if(err2!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda memcpy free d_sharedCols: "<<cudaGetErrorString(err2)<<std::endl;
		exit(-1);
	}

	if(err3!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err3)<<std::endl;
		exit(-1);
	}

	if(err4!=cudaSuccess){
		std::cout<<"GpuKernelHybridWithStreams37: cuda memcpy free d_sharedRows: "<<cudaGetErrorString(err4)<<std::endl;
		exit(-1);
	}
#endif

    __sync_bool_compare_and_swap(&FRAME(GpuFinish),false,true);
    uint64_t wlLeft = FRAME(wlLeft);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams37: GPU finish! and wlLeft: "<<wlLeft<<std::endl;
#endif
    ++FRAME(gpuCnt); 
    
    if(wlLeft == 0){
        SYNC(Swap37);
    }else{

        if(__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,false)&&(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false) )){ 
            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));
            double gpuStepR = FRAME(gpuStepR);
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams37: invoke new GPU"<<std::endl;
#endif
            
            pthread_mutex_lock(&mutex);
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0){
                SYNC(Swap37);
            }else{
             
                if(wlLeft < FRAME(lastCnt)+FRAME(cpuWLMin)){
                    FRAME(gpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*gpuWL;
                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = gpuWL;
                    }else{
                        t1 = (1+gpuStepR)*gpuWL;
                    }

                    if(wlLeft<t1){
                        FRAME(gpuWL) = wlLeft*rt;
                    }else{
                        FRAME(gpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
                }
                
                FRAME(nGPU)=nGPU;
                FRAME(gpuPos) = nRows-wlLeft ;
                __sync_synchronize();

                SYNC(GpuKernelHybridWithStreams37);
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelHybridWithStreams37: invoked a new GPU, gpuPos: "<<FRAME(gpuPos)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams37: invoked a new GPU, gpuWL: "<<FRAME(gpuWL)<<std::endl;

            std::cout<<"GpuKernelHybridWithStreams37: invoked a new GPU, wlLeft: "<<FRAME(wlLeft)<<std::endl;

#endif
        }
    }

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernelHybridWithStreams37: Swap dependence is "<<FRAME(Swap37).getCounter() <<std::endl;
#endif

//	SYNC(sync);
    
    EXIT_TP();
}


extern "C"
void Stencil3D7ptCpuSyncCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync37 invoke!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
    RESET(CpuSync37);

    ++FRAME(cpuCnt);

	uint64_t wlLeft = FRAME(wlLeft);

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"wlLeft:"<<wlLeft<<std::endl;
#endif
    if(wlLeft==0.0){
		SYNC(Swap37);
	
	}else{
		__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,true);
	
        size_t	gpuPos = FRAME(gpuPos);
		size_t	cpuPos = FRAME(cpuPos);
	
		uint64_t gpuWL = FRAME(gpuWL);
		uint64_t cpuWL = FRAME(cpuWL);
		uint64_t nCols = FRAME(nCols);
		uint64_t nRows = FRAME(nRows);
        uint64_t nSlices = FRAME(nSlices);
		double *dst = FRAME(New);
		double *d_dst = FRAME(d_dst);
		double cmCpu = FRAME(cmCpu);
		double cmGpu = FRAME(cmGpu);
        
        double cpuStepR = FRAME(cpuStepR);
        double gpuStepR = FRAME(gpuStepR);
        
        if(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false)) {

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync37: gpu kernel finish and invode a new gpu kernel"<<std::endl;
#endif

            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0.0){
                SYNC(Swap37);
            }else{
                
                if(wlLeft <= FRAME(lastCnt)+FRAME(gpuWLMin)){
                    FRAME(gpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{
                    uint64_t t1;
                  
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        uint64_t t2;
                        t2 = (1+gpuStepR)*gpuWL;
                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = gpuWL;
                    }else{
                        t1 = (1+gpuStepR)*gpuWL;
                    }
                    
                    if(wlLeft<t1){
                        FRAME(gpuWL) = wlLeft*rt;
                    }else{
                        FRAME(gpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
                }
                FRAME(nGPU)=2;
                FRAME(gpuPos) = FRAME(tWL)-wlLeft;
                __sync_synchronize();
                SYNC(GpuKernelHybridWithStreams37);
            }

#ifdef CUDA_DARTS_DEBUG
			size_t gpu_mem_total_t = 0;
			size_t gpu_mem_avail_t = 0;
			size_t gpu_mem_valid_t = 0;
			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
			gpu_mem_valid_t = gpu_mem_avail_t -XMB;

			std::cout<<"CpuSync37: gpu memory total: "<<gpu_mem_total_t<<std::endl;
			std::cout<<"CpuSync37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
			std::cout<<"CpuSync37: reqire gpu memory : "<<sizeof(double)*nCols*FRAME(gpuWL)<<std::endl;
#endif		


#ifdef CUDA_DARTS_DEBUG
				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: gpuPos: "<<FRAME(gpuPos)<<std::endl;
				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: gpuWL: "<<FRAME(gpuWL)<<std::endl;
				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: reset rowsLeft: "<<FRAME(wlLeft)<<std::endl;
				std::cout<<"CpuSync37: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;

#endif

		}
		if(FRAME(wlLeft)==0){
                SYNC(Swap37);
		}else{

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"CpuSync37: invoke new CpuLoop!"<<std::endl;
#endif
	
            double rt = FRAME(cmCpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
            pthread_mutex_lock(&mutex);
            wlLeft = FRAME(wlLeft);
            if(wlLeft==0){
                SYNC(Swap37);
            }else{
                 
                if(wlLeft <= FRAME(lastCnt)+FRAME(gpuWLMin)){
                    FRAME(cpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;
                }else{

                    uint64_t t1;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        t1 = (1-cpuStepR)*cpuWL;
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)){
                        t1 = cpuWL;
                    }else{
                        t1 = (1+cpuStepR)*cpuWL;
                    }

                    if(wlLeft<t1){
                        FRAME(cpuWL) = wlLeft*rt;
                    }else{
                        FRAME(cpuWL) = t1;
                    }
               
                    FRAME(wlLeft) = wlLeft - FRAME(cpuWL)+2;
                }

				FRAME(cpuPos) = nRows-wlLeft;
                __sync_synchronize();
            }
            pthread_mutex_unlock(&mutex);
            
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"CpuSync37: invoked  new CPU, cpuPos: "<<FRAME(cpuPos)/(FRAME(nRows)*FRAME(nCols))<<std::endl;

            std::cout<<"CpuSync37: invoked  new CPU, cpuWL: "<<FRAME(cpuWL)<<std::endl;

            std::cout<<"CpuSync37: invoked  new CPU, wlLeft: "<<FRAME(wlLeft)<<std::endl;

#endif
			__sync_bool_compare_and_swap(&FRAME(CpuFinish),true,false);
			for(size_t i =0; i<FRAME(nCPU);++i){
				SYNC(CpuLoop37[i]);
			}

#ifdef CUDA_DARTS_DEBUG
			std::cout<<"cpuSync37 new CpuLoop: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
#endif
		
		}

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync37: Swap dependence is "<<FRAME(Swap37).getCounter() <<std::endl;
#endif
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





