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


//	pthread_mutex_lock(&mutex);
//	pthread_mutex_unlock(&mutex);


extern "C"
void 
Stencil3D7ptCpuLoopCD::fire(void)
{
	LOAD_FRAME(StencilTP);
	uint64_t Id = getID();
	RESET(CpuLoop37[Id]);	
	
#ifdef CUDA_DARTS_DEBUG

    pthread_mutex_lock(&mutex2);
	std::cout<<"Invoke CpuLoop37["<<Id<<"]"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif
	double	*h_src		= FRAME(Initial);
	double	*h_dst		= FRAME(New);
	int32_t cpuSlices	= FRAME(cpuSlices);
	int32_t cpuRows		= FRAME(cpuRows);
	int32_t cpuCols		= FRAME(cpuCols);

	int cpuTile_x		= FRAME(cpuTile_x	);
	int cpuTile_y		= FRAME(cpuTile_y	);
    int cpuTile_z		= FRAME(cpuTile_z	);
    int cpuBlockDimx	= FRAME(cpuBlockDimx);
    int cpuBlockDimy	= FRAME(cpuBlockDimy);
    int cpuBlockDimz	= FRAME(cpuBlockDimz);
    int cpuGridDimx		= FRAME(cpuGridDimx	);
    int cpuGridDimy		= FRAME(cpuGridDimy	);
    int cpuGridDimz		= FRAME(cpuGridDimz	);

	int cpuPos			= FRAME(cpuPos);
	int nbxy = cpuGridDimx*cpuGridDimy;
	int posz = (Id/nbxy);
	int posy = (Id%nbxy)/cpuGridDimx;
	int posx = Id - posz*nbxy - posy*cpuGridDimx;

	uint64_t pos = posz*cpuBlockDimz*cpuRows*cpuCols+posy*cpuBlockDimy*cpuCols + posx*cpuBlockDimx;
	double *src = h_src+pos;
	double *dst = h_dst+pos;
	
	int nRowsChunk		= (((posy+1)*cpuBlockDimy)>=cpuRows)?(cpuRows-posy*cpuBlockDimy-1) : (cpuBlockDimy+1) ;
	int nColsChunk		= (((posx+1)*cpuBlockDimx)>=cpuCols)?(cpuCols-posx*cpuBlockDimx-1):(cpuBlockDimx+1) ;
	int nSlicesChunk	= (((posz+1)*cpuBlockDimz)>=cpuSlices)?(cpuSlices-posz*cpuBlockDimz-1):(cpuBlockDimz+1) ;

//	if((nSlicesChunk+cpuPos+chunk*Id)>=nSlices){
//			nSlicesChunk = cpuWL-chunk*Id-2 ;
//	}
//
//
#ifdef CUDA_DARTS_DEBUG

	pthread_mutex_lock(&mutex2);
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPos:"<<cpuPos<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posz:"<<posz<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posy:"<<posy<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posx:"<<posx<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nSlicesChunk:"<<nSlicesChunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nRowsChunk:"<<nRowsChunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nColsChunk:"<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif

//	//computeInner_stencil37(dst,src,nRows,nCols,nSlicesChunk);
	computeBlock_stencil37(dst,src,cpuRows,cpuCols,cpuSlices,nRowsChunk,nColsChunk,nSlicesChunk);


	if(FRAME(GpuRatio)==0.0){
		SYNC(Swap37);
	}else{
		SYNC(CpuSync37);
	}

#ifdef CUDA_DARTS_DEBUG
//	pthread_mutex_lock(&mutex2);
//	std::cout<<"CpuLoop37["<<Id<<"]: finish computing!"<<std::endl;
//	std::cout<<"CpuLoop37["<<Id<<"]: CpuSync37 dependence is "<<FRAME(CpuSync37).getCounter() <<std::endl;
//	pthread_mutex_unlock(&mutex2);
#endif
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
	std::cout<<"Swap37: ts: "<<FRAME(ts)<<std::endl;
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
//
//		double *d_dst = FRAME(d_dst);
//		double *dst = FRAME(New);
//		double *src = FRAME(Initial);
//
//		uint64_t	gpuPos = FRAME(gpuPos);
//		uint64_t	cpuPos = FRAME(cpuPos);
//		
//		uint64_t nRows   = FRAME(nRows);
//		uint64_t nCols   = FRAME(nCols);
//		uint64_t nSlices = FRAME(nSlices);
//        uint64_t gpuWL = FRAME(gpuWL);
//		uint64_t cpuWL = FRAME(cpuWL);
//    
//        if(ts!=0){
//                    
//            if(FRAME(gpuCnt)>FRAME(cpuCnt)){
//
//                FRAME(gpuWLInit) = FRAME(gpuWLInit)*(1+FRAME(gpuStepR));
//                FRAME(wlLeftInit) = FRAME(tWL) - FRAME(gpuWLInit)-FRAME(cpuWLInit) +4;
//                FRAME(cpuPosInit) = FRAME(gpuWLInit)-2;
//            } 
//            FRAME(gpuPos)	=FRAME(gpuPosInit);
//		    FRAME(cpuPos)	=FRAME(cpuPosInit);
//		    FRAME(gpuWL)	=FRAME(gpuWLInit);
//		    FRAME(cpuWL) =FRAME(cpuWLInit);
//		    FRAME(wlLeft)=FRAME(wlLeftInit);
//			FRAME(CpuFinish) = false;
//			FRAME(GpuFinish) = false;
//			FRAME(nCPU) = FRAME(nCPUInit);
//			FRAME(nGPU) = FRAME(nGPUInit);
//            FRAME(gpuCnt) = 0;
//            FRAME(cpuCnt) = 0;
//#ifdef CUDA_DARTS_DEBUG
//	        std::cout<<"swap37: reset gpuPos: "<<FRAME(gpuPos)<<std::endl;
//		    std::cout<<"swap37: reset cpuPos: "<<FRAME(cpuPos)<<std::endl;
//		    std::cout<<"swap37: reset gpu rows: "<<FRAME(gpuWL)<<std::endl;
//		    std::cout<<"swap37: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
//		    std::cout<<"swap37: reset cpu rowsLeft: "<<FRAME(wlLeft)<<std::endl;
//
//		    std::cout<<"swap37: reset CpuFinsh: "<<FRAME(CpuFinish)<<std::endl;
//
//		    std::cout<<"swap37: reset GpuFinsh: "<<FRAME(GpuFinish)<<std::endl;
//
//#endif
//
//		    SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
//
//            SYNC(GpuKernelHybridWithStreams37);
//		    for (size_t i =0;i<nCPU;++i){
//		    	SYNC(CpuLoop37[i]);
//		    }
//
//		}else{
//            SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
//            SYNC(sync);
//		}
//    
//    
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
	const uint64_t gpuRows   = FRAME(gpuRows);
	const uint64_t gpuCols   = FRAME(gpuCols);
    const uint64_t gpuSlices = FRAME(gpuSlices);
	
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
	size_t	pos1 = gpuPos*gpuRows*gpuCols;
	err1 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif

	int gpuTile_x = FRAME(gpuTile_x); //16
	int gpuTile_y = FRAME(gpuTile_y); //16
    int gpuTile_z = FRAME(gpuTile_z); //50

    int gpuBlockDimx = FRAME(gpuBlockDimx);   // tile_x
    int gpuBlockDimy = FRAME(gpuBlockDimy);   // tile_y
    int gpuBlockDimz = FRAME(gpuBlockDimz);   // 1

	int gpuGridDimx = FRAME(gpuGridDimx);     //  x/tile_x 
	int gpuGridDimy = FRAME(gpuGridDimy);     //  y/tile_y
    int gpuGridDimz = FRAME(gpuGridDimz);     //  z/tile_z
    
    dim3 dimGrid (gpuGridDimx, gpuGridDimy, gpuGridDimz);
	dim3 dimBlock(gpuBlockDimx,gpuBlockDimy,gpuBlockDimz);


#ifdef CUDA_DARTS_DEBUG
	std::cout<<"GpuKernelWithAllTimeSteps: GpuPos:"<<pos1/gpuCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nRows:"<<gpuRows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nCols:"<<gpuCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: nSlices:"<<gpuSlices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size:"<<d_size<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_x:"<<gpuTile_x<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_y:"<<gpuTile_y<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: tile_z:"<<gpuTile_z<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimx:"<<gpuBlockDimx<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimy:"<<gpuBlockDimy<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: blockDimz:"<<gpuBlockDimz<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimx="<<gpuGridDimx<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimy="<<gpuGridDimy<<std::endl;
    std::cout<<"GpuKernelWithAllTimeSteps: gridDimz="<<gpuGridDimz<<std::endl;

#endif
	
    int numThreads= gpuTile_x*gpuTile_y;

    int blockDimx_slices = (gpuCols>numThreads)?numThreads:gpuCols; 
	int blockDimy_slices = 1;
	int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*gpuCols/blockDimx_slices);
	int gridDimy_slices = gpuGridDimy;
	int gridDimz_slices = gpuGridDimz;

	int blockDimx_rows = (gpuCols>numThreads)?numThreads:gpuCols;
	int blockDimy_rows = 1;
	int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*gpuCols/blockDimx_rows);
	int gridDimy_rows = gpuGridDimy;
	int gridDimz_rows = gpuGridDimz;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = (gpuRows>numThreads)?numThreads:gpuRows;
	int blockDimz_cols = 1 ;
	int gridDimx_cols = gpuGridDimx;
	int gridDimy_cols = std::ceil(1.0*gpuRows/blockDimy_cols);
	int gridDimz_cols = gpuGridDimz;

#ifdef CUDA_DARTS_DEBUG

    std::cout<<"GpuKernelWithAllTimeSteps: blockDimx_slices:"<<blockDimx_slices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: blockDimy_slices:"<<blockDimy_slices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimx_slices:"<<gridDimx_slices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimy_slices:"<<gridDimy_slices<<std::endl;
	std::cout<<"GpuKernelWithAllTimeSteps: grimDimz_slices:"<<gridDimz_slices<<std::endl;
	
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

    
    int nStream = FRAME(nStream);
    cudaEvent_t *cuEvent = new cudaEvent_t[nStream];
    cudaError *err = new cudaError[nStream];
    for (int i=0;i<nStream;++i){
        err[i] = cudaEventCreate(&cuEvent[i]);
#ifdef CUDA_ERROR_CHECKING
	    if(err[i]!=cudaSuccess){
		    std::cout<<"GpuKernelWithAllTimeSteps: cuda create event ["<<i<<"]:"<<cudaGetErrorString(err[i])<<std::endl;
		    exit(-1);
	    }
#endif
    }

	size_t ts = FRAME(ts);
	while(ts-- >0){
        gpu_kernel37_cp_slices_stream(FRAME(stream[0]),dimGrid_slices,dimBlock_slices,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err1=cudaEventRecord(cuEvent[0],FRAME(stream[0]));
        gpu_kernel37_cp_rows_stream(FRAME(stream[1]),dimGrid_rows,dimBlock_rows,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err2=cudaEventRecord(cuEvent[1],FRAME(stream[1]));
        gpu_kernel37_cp_cols_stream(FRAME(stream[2]),dimGrid_cols,dimBlock_cols,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err3=cudaEventRecord(cuEvent[2],FRAME(stream[2]));
        
#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaEventRecord stream0 errors: "<<cudaGetErrorString(err1)<<std::endl;
	    	exit(-1);
	    }
	    if(err2!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaEventRecord stream1 errors: "<<cudaGetErrorString(err2)<<std::endl;
	    	exit(-1);
	    }
	    if(err3!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaEventRecord stream2 errors: "<<cudaGetErrorString(err3)<<std::endl;
	    	exit(-1);
	    }
#endif
        err1=cudaStreamWaitEvent(FRAME(stream[3]),cuEvent[0],0);
        err2=cudaStreamWaitEvent(FRAME(stream[3]),cuEvent[1],0);
        err3=cudaStreamWaitEvent(FRAME(stream[3]),cuEvent[2],0);
        
#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream0 errors: "<<cudaGetErrorString(err1)<<std::endl;
	    	exit(-1);
	    }
	    if(err2!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream1 errors: "<<cudaGetErrorString(err2)<<std::endl;
	    	exit(-1);
	    }
	    if(err3!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream2 errors: "<<cudaGetErrorString(err3)<<std::endl;
	    	exit(-1);
	    }
#endif
        gpu_kernel37_stream(FRAME(stream[3]),dimGrid,dimBlock,d_dst,d_sharedRows,d_sharedCols,d_sharedSlices,gpuRows,gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);

        err4=cudaEventRecord(cuEvent[3],FRAME(stream[3]));

#ifdef CUDA_ERROR_CHECKING
	    if(err4!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaEventRecord stream3 errors: "<<cudaGetErrorString(err4)<<std::endl;
	    	exit(-1);
	    }
#endif
        err1=cudaStreamWaitEvent(FRAME(stream[0]),cuEvent[3],0);
        err2=cudaStreamWaitEvent(FRAME(stream[1]),cuEvent[3],0);
        err3=cudaStreamWaitEvent(FRAME(stream[2]),cuEvent[3],0);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream0 errors: "<<cudaGetErrorString(err1)<<std::endl;
	    	exit(-1);
	    }
	    if(err2!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream1 errors: "<<cudaGetErrorString(err2)<<std::endl;
	    	exit(-1);
	    }
	    if(err3!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaStreamWaitEvent stream2 errors: "<<cudaGetErrorString(err3)<<std::endl;
	    	exit(-1);
	    }
#endif
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
    //destroy cudaEvent
    for(int i=0;i<nStream;++i){
        err[i]=cudaEventDestroy(cuEvent[i]);
        
#ifdef CUDA_ERROR_CHECKING
	    if(err[i]!=cudaSuccess){
		    std::cout<<"GpuKernelWithAllTimeSteps: cuda destroy event ["<<i<<"]:"<<cudaGetErrorString(err[i])<<std::endl;
		    exit(-1);
	    }
#endif
    }
    delete [] cuEvent;
    delete [] err;
	//copy from GPU  to CPU
        err1=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyDeviceToHost: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif

//#ifdef CUDA_DARTS_DEBUG
//
//
//#endif
//
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

/*
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

    int blockDimx_slices = (nCols>numThreads)?numThreads:nCols; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/blockDimx_slices);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz;

    int blockDimx_rows = (nCols>numThreads)?numThreads:nCols;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = (nRows>numThreads)?numThreads:nRows;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/blockDimy_cols);
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

    int blockDimx_slices = (nCols>numThreads)?numThreads:nCols; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/blockDimx_slices);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz2;

    int blockDimx_rows = (nCols>numThreads)?numThreads:nCols;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/blockDimy_rows);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz2;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = (nRows>numThreads)?numThreads:nRows;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/blockDimy_cols);
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
//	std::cout<<"GpuKernelHybridWithStreams: gpu memory total: "<<gpu_mem_total_t<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37 : require memory size:"<<(d_size + d_size_sharedCols + d_size_sharedRows)*4<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: nGPU:"<<nGPU<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: Gpupos:"<<FRAME(gpuPos)<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: nRows:"<<nRows<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: nCols:"<<nCols<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: nSlices:"<<nSlices<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: gpuWL:"<<gpuWL<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: nTile_z:"<<nTile_z<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: nSlicesChunk:"<<nSlicesChunk<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: d_size:"<<d_size<<std::endl;
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

    int blockDimx_slices = (nCols>numThreads)?numThreads:nCols; 
    int blockDimy_slices = 1;
    int blockDimz_slices = 1;
	int gridDimx_slices = std::ceil(1.0*nCols/blockDimx_slices);
	int gridDimy_slices = gridDimy;
    int gridDimz_slices = gridDimz2;

    int blockDimx_rows = (nCols>numThreads)?numThreads:nCols;
	int blockDimy_rows = 1;
    int blockDimz_rows = 1;
	int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	int gridDimy_rows = gridDimy;
    int gridDimz_rows = gridDimz2;

	int blockDimx_cols = 1 ;
	int blockDimy_cols = (nRows>numThreads)?numThreads:nRows;
	int blockDimz_cols = 1 ;
    int gridDimx_cols = gridDimx;
	int gridDimy_cols = std::ceil(1.0*nRows/blockDimy_cols);
    int gridDimz_cols = gridDimz2;


#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"GpuKernelHybridWithStreams37: nRows:"<<nRows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: nCols:"<<nCols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: nSlices:"<<nSlices<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: d_size:"<<d_size<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: tile_x:"<<tile_x<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: tile_y:"<<tile_y<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: tile_z:"<<tile_z<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: blockDimx:"<<blockDimx<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: blockDimy:"<<blockDimy<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: blockDimz:"<<blockDimz<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: gridDimx="<<gridDimx<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: gridDimy="<<gridDimy<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: gridDimz="<<gridDimz<<std::endl;
//
//	std::cout<<"GpuKernelHybridWithStreams37: blockDimx_slices:"<<blockDimx_slices<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_slices:"<<blockDimy_slices<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: grimDimx_slices:"<<gridDimx_slices<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_slices:"<<gridDimy_slices<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_slices:"<<gridDimz_rows<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: blockDimx_rows:"<<blockDimx_rows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_rows:"<<blockDimy_rows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimx_rows:"<<gridDimx_rows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_rows:"<<gridDimy_rows<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_rows:"<<gridDimz_rows<<std::endl;
//    std::cout<<"GpuKernelHybridWithStreams37: blockDimx_cols:"<<blockDimx_cols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: blockDimy_cols:"<<blockDimy_cols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimx_cols:"<<gridDimx_cols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimy_cols:"<<gridDimy_cols<<std::endl;
//	std::cout<<"GpuKernelHybridWithStreams37: grimDimz_cols:"<<gridDimz_cols<<std::endl;
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
	            std::cout<<"GpuKernelHybridWithStreams37: blockDimz_slices:"<<blockDimz_slices<<std::endl;
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
	            std::cout<<"GpuKernelHybridWithStreams37: h_pos:"<<h_pos<<std::endl;
	            std::cout<<"GpuKernelHybridWithStreams37: d_size:"<<nSlicesLeft<<std::endl;
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
             
                if(wlLeft < FRAME(lastCnt)+FRAME(gpuWLMin)){
                    FRAME(gpuWL) = wlLeft;
                    FRAME(wlLeft) = 0;

                    FRAME(nGPU)=1;
                }else{

                    uint64_t t1,t2;
                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
                        t2 = (1+gpuStepR)*gpuWL;
                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
                        t1 = gpuWL;
                    }else{
                        t2 = (1-gpuStepR)*gpuWL;
                        t1 = t2<FRAME(gpuWLMin)?FRAME(gpuWLMin):t2; 
                    }

                    if(wlLeft<t1){
                        FRAME(gpuWL) = wlLeft;
                        FRAME(wlLeft) = 0;
                    }else{
                        FRAME(gpuWL) = t1;
                        FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
                    }
               
                    FRAME(nGPU)=nGPU;
                }
                
                FRAME(gpuPos) = FRAME(tWL)-wlLeft ;
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
//	std::cout<<"GpuKernelHybridWithStreams37: Swap dependence is "<<FRAME(Swap37).getCounter() <<std::endl;
#endif

//	SYNC(sync);
    
    EXIT_TP();
}

*/
extern "C"
void Stencil3D7ptCpuSyncCD::fire(void)
{

#ifdef CUDA_DARTS_DEBUG
	std::cout<<"CpuSync37 invoke!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
    RESET(CpuSync37);

//    ++FRAME(cpuCnt);
//
//	uint64_t wlLeft = FRAME(wlLeft);
//
//#ifdef CUDA_DARTS_DEBUG
//	std::cout<<"wlLeft:"<<wlLeft<<std::endl;
//#endif
//    if(wlLeft==0.0){
//		SYNC(Swap37);
//	
//	}else{
//		__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,true);
//	
//        size_t	gpuPos = FRAME(gpuPos);
//		size_t	cpuPos = FRAME(cpuPos);
//	
//		uint64_t gpuWL = FRAME(gpuWL);
//		uint64_t cpuWL = FRAME(cpuWL);
//		uint64_t nCols = FRAME(nCols);
//		uint64_t nRows = FRAME(nRows);
//        uint64_t nSlices = FRAME(nSlices);
//		double *dst = FRAME(New);
//		double *d_dst = FRAME(d_dst);
//		double cmCpu = FRAME(cmCpu);
//		double cmGpu = FRAME(cmGpu);
//        
//        double cpuStepR = FRAME(cpuStepR);
//        double gpuStepR = FRAME(gpuStepR);
//        
//        if(__sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false)) {
//
//#ifdef CUDA_DARTS_DEBUG
//			std::cout<<"cpuSync37: gpu kernel finish and invode a new gpu kernel"<<std::endl;
//#endif
//
//            double rt = FRAME(cmGpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
//            wlLeft = FRAME(wlLeft);
//            if(wlLeft==0.0){
//                SYNC(Swap37);
//            }else{
//                
//                if(wlLeft <= FRAME(lastCnt)+FRAME(gpuWLMin)){
//                    FRAME(gpuWL) = wlLeft;
//                    FRAME(wlLeft) = 0;
//                    FRAME(nGPU) = 1;
//                }else{
//                    uint64_t t1;
//                    uint64_t t2;
//                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
//                        t2 = (1+gpuStepR)*gpuWL;
//                        t1 = (t2<FRAME(gpuWLMax))?t2:FRAME(gpuWLMax);
//                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)) {
//                        t1 = gpuWL;
//                    }else{
//                        t2 = (1-gpuStepR)*gpuWL;
//                        t1 = (t2<FRAME(gpuWLMin))?FRAME(gpuWLMin):t2;
//                    }
//                    
//                    if(wlLeft<t1){
//                        FRAME(gpuWL) = wlLeft;
//                        FRAME(wlLeft)=0;
//                    }else{
//                        FRAME(gpuWL) = t1;
//                
//                        FRAME(wlLeft) = wlLeft - FRAME(gpuWL)+2;
//                    }
//                }
//                FRAME(gpuPos) = FRAME(tWL)-wlLeft;
//                __sync_synchronize();
//                SYNC(GpuKernelHybridWithStreams37);
//            }
//
//#ifdef CUDA_DARTS_DEBUG
//			size_t gpu_mem_total_t = 0;
//			size_t gpu_mem_avail_t = 0;
//			size_t gpu_mem_valid_t = 0;
//			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
//			gpu_mem_valid_t = gpu_mem_avail_t -XMB;
//
//			std::cout<<"CpuSync37: gpu memory total: "<<gpu_mem_total_t<<std::endl;
//			std::cout<<"CpuSync37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
//#endif		
//
//
//#ifdef CUDA_DARTS_DEBUG
//				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: gpuPos: "<<FRAME(gpuPos)<<std::endl;
//				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: gpuWL: "<<FRAME(gpuWL)<<std::endl;
//				std::cout<<"CpuSync37 new GpuKernelHybridWithStreams: reset rowsLeft: "<<FRAME(wlLeft)<<std::endl;
////				std::cout<<"CpuSync37: Swap dependence is "<<FRAME(Swap).getCounter() <<std::endl;
//
//#endif
//
//		}
//		if(FRAME(wlLeft)==0){
//                SYNC(Swap37);
//		}else{
//
//#ifdef CUDA_DARTS_DEBUG
//			std::cout<<"CpuSync37: invoke new CpuLoop!"<<std::endl;
//#endif
//	
//            double rt = FRAME(cmCpu)/(FRAME(cmCpu)+FRAME(cmGpu));          
//            pthread_mutex_lock(&mutex);
//            wlLeft = FRAME(wlLeft);
//            if(wlLeft==0){
//                SYNC(Swap37);
//            }else{
//                 
//                if(wlLeft <= FRAME(cpuWLMin)){
//                //if(wlLeft <= FRAME(lastCnt)+FRAME(cpuWLMin)){
//                    FRAME(cpuWL) = wlLeft;
//                    FRAME(wlLeft) = 0;
//					FRAME(nCPU) = wlLeft-1;
//				}else if(wlLeft <= FRAME(lastCnt)+FRAME(cpuWLMin)){
//
//                    FRAME(cpuWL) = wlLeft;
//                    FRAME(wlLeft) = 0;
//				}else{
//
//                    uint64_t t1;
//                    if(FRAME(gpuCnt)>FRAME(cpuCnt)){
//                        t1 = (1-cpuStepR)*cpuWL;
//                    }else if(FRAME(gpuCnt)==FRAME(cpuCnt)){
//                        t1 = cpuWL;
//                    }else{
//                        t1 = (1+cpuStepR)*cpuWL;
//                    }
//
//                    if(wlLeft<t1){
//                        FRAME(cpuWL) = wlLeft;
//                        FRAME(wlLeft) = 0;
//                    }else{
//                        FRAME(cpuWL) = t1;
//                        FRAME(wlLeft) = wlLeft - FRAME(cpuWL)+2;
//                    }
//               
//                }
//
//				FRAME(cpuPos) = FRAME(tWL)-wlLeft;
//                __sync_synchronize();
//            }
//            pthread_mutex_unlock(&mutex);
//            
//#ifdef CUDA_DARTS_DEBUG
//            std::cout<<"CpuSync37: invoked  new CPU, cpuPos: "<<FRAME(cpuPos)<<std::endl;
//
//            std::cout<<"CpuSync37: invoked  new CPU, cpuWL: "<<FRAME(cpuWL)<<std::endl;
//
//            std::cout<<"CpuSync37: invoked  new CPU, wlLeft: "<<FRAME(wlLeft)<<std::endl;
//
//            std::cout<<"CpuSync37: invoked  new CPU, nCPU: "<<FRAME(nCPU)<<std::endl;
//#endif
//			__sync_bool_compare_and_swap(&FRAME(CpuFinish),true,false);
//			for(size_t i =0; i<FRAME(nCPU);++i){
//				SYNC(CpuLoop37[i]);
//			}
//
//#ifdef CUDA_DARTS_DEBUG
//			std::cout<<"cpuSync37 new CpuLoop: reset cpu rows: "<<FRAME(cpuWL)<<std::endl;
//#endif
//		
//		}
//
//#ifdef CUDA_DARTS_DEBUG
////	std::cout<<"CpuSync37: Swap dependence is "<<FRAME(Swap37).getCounter() <<std::endl;
//#endif
//	}

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





