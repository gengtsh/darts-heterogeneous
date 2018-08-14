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
    int     *arrnCpuPos = FRAME(arrnCpuPos);
    int *arrnEdge  = FRAME(arrnEdge);
    int *arrnCpuEdge = FRAME(arrnCpuEdge);
    int *arrnCpuBlock = FRAME(arrnCpuBlock);
    int *arrnCpuGrid  = FRAME(arrnCpuGrid);

	const int nCols		= FRAME(nCols);
	const int nRows		= FRAME(nRows);
    const int nSlices	= FRAME(nSlices);

	int cpuCols		= arrnCpuEdge[0];
	int cpuRows		= arrnCpuEdge[1];
    int cpuSlices	= arrnCpuEdge[2];

    int cpuBlockDimx	= arrnCpuBlock[0];
    int cpuBlockDimy	= arrnCpuBlock[1];
    int cpuBlockDimz	= arrnCpuBlock[2];
    int cpuGridDimx		= arrnCpuGrid[0];
    int cpuGridDimy		= arrnCpuGrid[1];
    int cpuGridDimz		= arrnCpuGrid[2];

	//int cpuPos			= FRAME(cpuPos);
    int cpuPosx         = arrnCpuPos[0]; 
    int cpuPosy         = arrnCpuPos[1];
    int cpuPosz         = arrnCpuPos[2];


    int  arrnCpuEdge2[DIM];

    calcEdge(arrnCpuEdge2,arrnCpuPos,arrnCpuEdge, arrnEdge, 1,0,0,2);
	
    int cpuCols2		= arrnCpuEdge2[0];
	int cpuRows2		= arrnCpuEdge2[1];
    int cpuSlices2  	= arrnCpuEdge2[2];
    
    int cpuPos          = cpuPosx+cpuPosy*nCols+cpuPosz*nCols*nRows;
	
    int nbxy = cpuGridDimx*cpuGridDimy;
	int posz = (Id/nbxy);
	int posy = (Id%nbxy)/cpuGridDimx;
	int posx = Id - posz*nbxy - posy*cpuGridDimx;

	//uint64_t pos = posz*cpuBlockDimz*cpuRows*cpuCols+posy*cpuBlockDimy*cpuCols + posx*cpuBlockDimx;
	uint64_t pos = posz*cpuBlockDimz*nRows*nCols+posy*cpuBlockDimy*nCols + posx*cpuBlockDimx;
	double *src = h_src+cpuPos+pos;
	double *dst = h_dst+cpuPos+pos;
	
	int nRowsChunk		= (((posy+1)*cpuBlockDimy)>=cpuRows2)?(cpuRows2-posy*cpuBlockDimy-1) : (cpuBlockDimy+1) ;
	int nColsChunk		= (((posx+1)*cpuBlockDimx)>=cpuCols2)?(cpuCols2-posx*cpuBlockDimx-1):(cpuBlockDimx+1) ;
	int nSlicesChunk	= (((posz+1)*cpuBlockDimz)>=cpuSlices2)?(cpuSlices2-posz*cpuBlockDimz-1):(cpuBlockDimz+1) ;
#ifdef CUDA_DARTS_DEBUG

	pthread_mutex_lock(&mutex2);
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPos:"<<cpuPos<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPosx:"<<cpuPosx<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPosy:"<<cpuPosy<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: cpuPosz:"<<cpuPosz<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posz:"<<posz<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posy:"<<posy<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: posx:"<<posx<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nSlicesChunk:"<<nSlicesChunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nRowsChunk:"<<nRowsChunk<<std::endl;
	std::cout<<"CpuLoop37["<<Id<<"]: nColsChunk:"<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif

	computeBlock_stencil37(dst,src,nRows,nCols,nSlices,nRowsChunk,nColsChunk,nSlicesChunk);


	if(FRAME(GpuRatio)==0.0){
		SYNC(Swap37);
	}else{
		SYNC(CpuSync37);
	}

#ifdef CUDA_DARTS_DEBUG
	//pthread_mutex_lock(&mutex2);
	//std::cout<<"CpuLoop37["<<Id<<"]: finish computing!"<<std::endl;
	//std::cout<<"CpuLoop37["<<Id<<"]: CpuSync37 dependence is "<<FRAME(CpuSync37).getCounter() <<std::endl;
	//pthread_mutex_unlock(&mutex2);
#endif
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
   
    int nId = FRAME(nId);
    int nIdPlus1 = R(nId+1);
    int nIdMinus1 = R(nId-1);

    int *arrnEdge           = FRAME(arrnEdge);
    int *arrnCpuEdge        = FRAME(arrnCpuEdge);
    int *arrnCpuEdgeMin     = FRAME(arrnCpuEdgeMin);
    int *arrnGpuEdgeMin     = FRAME(arrnGpuEdgeMin);
    //============================= set cpu finish sign ================================//
	__sync_bool_compare_and_swap(&FRAME(CpuFinish),false,true);

    //============================= record cpu execution time ================================//
#ifdef TIMERECORD
    int cnt = RT(FRAME(cpuCnt));
    int *chunk = FRAME(scRec[cnt]).chunk;
    setarrn1Fromarrn2(chunk,arrnCpuEdge);
    FRAME(scRec[cnt]).tEnd = getTime();
    FRAME(scRec[cnt]).tExe=FRAME(scRec[cnt]).tEnd - FRAME(scRec[cnt]).tStart;
#endif

    FRAME(cpuCnt++);
    cnt = RT(FRAME(cpuCnt));
    FRAME(scRec[cnt]).tStart = getTime();

    //============================= calc cpu edge and invoke new cpu/swap codelet ================================//
    
    //int *arrnCpuEdge        = FRAME(arrnCpuEdge);
    int *arrnGpuEdge        = FRAME(arrnGpuEdge);
    int *arrnCpuEdgeLeft    = FRAME(arrnCpuEdgeLeft); 
    int *arrnGpuEdgeLeft    = FRAME(arrnGpuEdgeLeft); 
    int *arrnCpuPos         = FRAME(arrnCpuPos);
    int *arrnGpuPos         = FRAME(arrnGpuPos);
    int *arrnCpuEdgeVar     = FRAME(arrnCpuEdgeVar);
    int *arrnGpuEdgeVar     = FRAME(arrnGpuEdgeVar);
    int cpuCnt              = FRAME(cpuCnt);
    int gpuCnt              = FRAME(gpuCnt);
  

    int *arrnCpuTile = FRAME(arrnCpuTile); 
    int *arrnCpuGridTileBase = FRAME(arrnCpuGridTileBase);
    int *arrnCpuBlock = FRAME(arrnCpuBlock);
    int *arrnCpuGrid  = FRAME(arrnCpuGrid);


#ifdef CUDA_DARTS_DEBUG2
	pthread_mutex_lock(&mutex2);

    std::cout<<"cpu: current cpuCnt = "<<FRAME(cpuCnt)<<std::endl;
    std::cout<<"cpu: current gpuCnt = "<<FRAME(gpuCnt)<<std::endl;
    std::cout<<"cpu: current nCPU = "<<FRAME(nCPU)<<std::endl;
    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current cpuPos["<<i<<"] = "<<arrnCpuPos[i]<<std::endl;
    }
    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current cpuEdge["<<i<<"] = "<<arrnCpuEdge[i]<<std::endl;
    }
    
    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current cpuEdgeLeft["<<i<<"] = "<<arrnCpuEdgeLeft[i]<<std::endl;
    }
    
    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current gpuPos["<<i<<"] = "<<arrnGpuPos[i]<<std::endl;
    }

    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current gpuEdge["<<i<<"] = "<<arrnGpuEdge[i]<<std::endl;
    }
    for (int i=0;i<DIM;++i){
        std::cout<<"cpu: current gpuEdgeLeft["<<i<<"] = "<<arrnGpuEdgeLeft[i]<<std::endl;
    }
    //for (int i=0;i<DIM;++i){
    //    std::cout<<"current cpuBlock["<<i<<"] = "<<arrnCpuBlock[i]<<std::endl;
    //}
    //
    //for (int i=0;i<DIM;++i){
    //    std::cout<<"current cpuGrid["<<i<<"] = "<<arrnCpuGrid[i]<<std::endl;
    //}
    pthread_mutex_unlock(&mutex2);
#endif


    if( checkarrnEqValue(arrnCpuEdgeLeft, 0) && checkarrnEqValue(arrnGpuEdgeLeft,0)){//calculation finish

        //============================= invoke Swap codelet ================================//
		SYNC(Swap37);
    
#ifdef CUDA_DARTS_DEBUG
        std::cout<<"cpu calculation finish!"<<std::endl;
#endif
    
    }else{

        //============================= calc cpu edge  =====================================//
        pthread_mutex_lock(&mutex);
        __sync_synchronize();
        calcNextEP(arrnCpuEdge, arrnGpuEdge,arrnCpuEdgeLeft, arrnGpuEdgeLeft,arrnCpuPos, arrnGpuPos,arrnEdge, arrnCpuEdgeVar, arrnGpuEdgeVar,arrnCpuEdgeMin,arrnGpuEdgeMin,cpuCnt,gpuCnt, nId,"cpu");
        __sync_synchronize();
        pthread_mutex_unlock(&mutex);

        //============================= invoke new cpu codelet ==============================//
        chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
        chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
        calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
		int nCPU = arrnCpuGrid[0]*arrnCpuGrid[1]*arrnCpuGrid[2];
        FRAME(nCPU) = nCPU;
        SETSYNC(CpuSync37,nCPU,nCPU);
   
#ifdef CUDA_DARTS_DEBUG
	    pthread_mutex_lock(&mutex2);
        std::cout<<"next cpuSync37 dep: "<<FRAME(CpuSync37).getCounter()<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif
        
#ifdef CUDA_DARTS_DEBUG2
	    pthread_mutex_lock(&mutex2);

        std::cout<<"cpu: next nCPU = "<<FRAME(nCPU)<<std::endl;
        
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next cpuPos["<<i<<"] = "<<arrnCpuPos[i]<<std::endl;
        }
	    
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next cpuEdge["<<i<<"] = "<<arrnCpuEdge[i]<<std::endl;
        }
        
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next cpuEdgeLeft["<<i<<"] = "<<arrnCpuEdgeLeft[i]<<std::endl;
        }
       
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next gpuPos["<<i<<"] = "<<arrnGpuPos[i]<<std::endl;
        }
	    
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next gpuEdge["<<i<<"] = "<<arrnGpuEdge[i]<<std::endl;
        }
        
        for (int i=0;i<DIM;++i){
            std::cout<<"cpu: next gpuEdgeLeft["<<i<<"] = "<<arrnGpuEdgeLeft[i]<<std::endl;
        }

        //for (int i=0;i<DIM;++i){
        //    std::cout<<"next cpuBlock["<<i<<"] = "<<arrnCpuBlock[i]<<std::endl;
        //}
        //
        //for (int i=0;i<DIM;++i){
        //    std::cout<<"next cpuGrid["<<i<<"] = "<<arrnCpuGrid[i]<<std::endl;
        //}
        pthread_mutex_unlock(&mutex2);
#endif
        for(size_t i=0;i<nCPU; ++i){
	        SYNC(CpuLoop37[i]);
	    }
        
        //============================= reset dep of CpuSync codelet ==============================//
    
    }


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


	SWAP_PTR(&FRAME(New) ,&FRAME(Initial) );
	size_t ts = --FRAME(ts);
	if(GpuRatio == 0.0){

		if(ts!=0){

			for (size_t i =0;i<nCPU;++i){
				SYNC(CpuLoop37[i]);
			}
		}else{
			SYNC(sync);
		}
	}else{

        if(ts!=0){

            int *arrnEdge               = FRAME(arrnEdge);
            
            int *arrnCpuPos             = FRAME(arrnCpuPos);
            int *arrnCpuEdge            = FRAME(arrnCpuEdge);
            int *arrnCpuEdgeMin         = FRAME(arrnCpuEdgeMin);
            int *arrnCpuEdgeLeft        = FRAME(arrnCpuEdgeLeft); 
            int *arrnCpuEdgeVar         = FRAME(arrnCpuEdgeVar);
            int *arrnCpuTile            = FRAME(arrnCpuTile); 
            int *arrnCpuGridTileBase    = FRAME(arrnCpuGridTileBase);
            int *arrnCpuBlock           = FRAME(arrnCpuBlock);
            int *arrnCpuGrid            = FRAME(arrnCpuGrid);
            int cpuCnt                  = FRAME(cpuCnt);
     

            int *arrnGpuPos             = FRAME(arrnGpuPos);
            int *arrnGpuEdge            = FRAME(arrnGpuEdge);
            int *arrnGpuEdgeMin         = FRAME(arrnGpuEdgeMin);
            int *arrnGpuEdgeLeft        = FRAME(arrnGpuEdgeLeft); 
            int *arrnGpuEdgeVar         = FRAME(arrnGpuEdgeVar);
            int *arrnGpuTile            = FRAME(arrnGpuTile); 
            int *arrnGpuGridTileBase    = FRAME(arrnGpuGridTileBase);
            int *arrnGpuBlock           = FRAME(arrnGpuBlock);
            int *arrnGpuGrid            = FRAME(arrnGpuGrid);
            int gpuCnt                  = FRAME(gpuCnt);
           
            int nId   = R(FRAME(nId));
            int nIdP1 = R(nId+1);
            int nIdM1 = R(nId-1);

            //============================= calc cpu gpu EdgeLeft=====================================//
            setarrn1Fromarrn2(arrnCpuEdgeLeft, arrnEdge); 
            setarrn1Fromarrn2(arrnGpuEdgeLeft, arrnEdge); 
            
            //============================= calc cpu pos  =====================================//
            setarrnValue(arrnCpuPos,0); 
            //============================= calc cpu edge  =====================================//
            calcNextEEL(arrnCpuEdge, arrnCpuEdgeLeft,arrnEdge, arrnCpuEdgeVar,arrnCpuEdgeMin,nId, cpuCnt, gpuCnt);
            calcNextEEL(arrnCpuEdge, arrnCpuEdgeLeft,arrnEdge, arrnCpuEdgeVar,arrnCpuEdgeMin,nIdP1, cpuCnt, gpuCnt);
            calcNextEEL(arrnCpuEdge, arrnCpuEdgeLeft,arrnEdge, arrnCpuEdgeVar,arrnCpuEdgeMin,nIdM1, cpuCnt, gpuCnt);
            //============================= calc gpu rnId EdgeLeft = cpu rnId Edgeleft==========//
            arrnGpuEdgeLeft[nId] = arrnCpuEdgeLeft[nId];

            //============================= calc gpu edge  =====================================//
            calcNextEEL(arrnGpuEdge, arrnGpuEdgeLeft,arrnEdge, arrnGpuEdgeVar,arrnGpuEdgeMin,nId, gpuCnt, cpuCnt);
            calcNextEEL(arrnGpuEdge, arrnGpuEdgeLeft,arrnEdge, arrnGpuEdgeVar,arrnGpuEdgeMin,nIdP1, gpuCnt, cpuCnt);
            calcNextEEL(arrnGpuEdge, arrnGpuEdgeLeft,arrnEdge, arrnGpuEdgeVar,arrnGpuEdgeMin,nIdM1, gpuCnt, cpuCnt);
            //============================= calc gpu pos  =====================================//
            calcarrn1Fromarrn2Minusarrn3(arrnGpuPos,arrnEdge,arrnGpuEdge);
            //============================= calc cpu rnId EdgeLeft = cpu gnId Edgeleft==========//
            arrnCpuEdgeLeft[nId] = arrnGpuEdgeLeft[nId];

            //============================= reset cpuCnt and gpuCnt ==============================//
            FRAME(cpuCnt)=0;
            FRAME(gpuCnt)=0;

            //============================= invoke new cpu codelet ==============================//
            chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
            chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
            calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
		    int nCPU = arrnCpuGrid[0]*arrnCpuGrid[1]*arrnCpuGrid[2];
            FRAME(nCPU) = nCPU;
            SETSYNC(CpuSync37,nCPU,nCPU);

            for(size_t i=0;i<nCPU; ++i){
	            SYNC(CpuLoop37[i]);
	        }
            //============================= invoke new gpu codelet ==============================//
            SYNC(GpuKernelHybridWithStreams37);
        
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

	const uint64_t gpuCols   = FRAME(arrnGpuEdge[0]);
    const uint64_t gpuRows   = FRAME(arrnGpuEdge[1]);
    const uint64_t gpuSlices = FRAME(arrnGpuEdge[2]);
	
	const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);


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
    int gpuPosx = FRAME( arrnGpuPos[0]); 
    int gpuPosy = FRAME( arrnGpuPos[1]); 
    int gpuPosz = FRAME( arrnGpuPos[2]); 

	//size_t	pos1 = gpuPos*gpuRows*gpuCols;
    size_t  pos1 = gpuPosx+gpuPosy*nCols+gpuPosz*nRows*nCols;
	err1 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);
#ifdef CUDA_ERROR_CHECKING
	if(err1!=cudaSuccess){
		std::cout<<"GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst: "<<cudaGetErrorString(err1)<<std::endl;
		exit(-1);
	}
#endif

	int gpuTile_x = FRAME(arrnGpuTile[0]); //16
	int gpuTile_y = FRAME(arrnGpuTile[1]); //16
    int gpuTile_z = FRAME(arrnGpuTile[2]); //50

    int gpuBlockDimx = FRAME(arrnGpuBlock[0]);   // tile_x
    int gpuBlockDimy = FRAME(arrnGpuBlock[1]);   // tile_y
    int gpuBlockDimz = FRAME(arrnGpuBlock[2]);   // 1

	int gpuGridDimx = FRAME(arrnGpuGrid[0]);     //  x/tile_x 
	int gpuGridDimy = FRAME(arrnGpuGrid[1]);     //  y/tile_y
    int gpuGridDimz = FRAME(arrnGpuGrid[2]);     //  z/tile_z
    
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

    int tnStream = FRAME(tnStream);
    cudaStream_t *stream ;
	stream = new cudaStream_t[tnStream];
    for(int i=0;i<tnStream;++i){
        cudaStreamCreate(&stream[i]);
    }

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
        gpu_kernel37_cp_slices_stream(stream[0],dimGrid_slices,dimBlock_slices,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err1=cudaEventRecord(cuEvent[0],stream[0]);
        gpu_kernel37_cp_rows_stream(stream[1],dimGrid_rows,dimBlock_rows,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err2=cudaEventRecord(cuEvent[1],stream[1]);
        gpu_kernel37_cp_cols_stream(stream[2],dimGrid_cols,dimBlock_cols,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, gpuRows, gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);
        err3=cudaEventRecord(cuEvent[2],stream[2]);
        
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
        err1=cudaStreamWaitEvent(stream[3],cuEvent[0],0);
        err2=cudaStreamWaitEvent(stream[3],cuEvent[1],0);
        err3=cudaStreamWaitEvent(stream[3],cuEvent[2],0);
        
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
        gpu_kernel37_stream(stream[3],dimGrid,dimBlock,d_dst,d_sharedRows,d_sharedCols,d_sharedSlices,gpuRows,gpuCols,gpuSlices,gpuTile_x,gpuTile_y,gpuTile_z);

        err4=cudaEventRecord(cuEvent[3],stream[3]);

#ifdef CUDA_ERROR_CHECKING
	    if(err4!=cudaSuccess){
	    	std::cout<<"GpuKernelWithAllTimeSteps: cudaEventRecord stream3 errors: "<<cudaGetErrorString(err4)<<std::endl;
	    	exit(-1);
	    }
#endif
        err1=cudaStreamWaitEvent(stream[0],cuEvent[3],0);
        err2=cudaStreamWaitEvent(stream[1],cuEvent[3],0);
        err3=cudaStreamWaitEvent(stream[2],cuEvent[3],0);

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

	for(int i=0;i<tnStream;++i){
		cudaStreamDestroy(stream[i]);
	}
    delete [] stream;
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
	
#ifdef CUDA_DARTS_DEBUG
    std::cout<<"h_src addr: "<<h_src<<std::endl;
    std::cout<<"h_dst addr: "<<h_dst<<std::endl;
#endif

	const uint64_t gpuCols   = FRAME(arrnGpuEdge[0]);
    const uint64_t gpuRows   = FRAME(arrnGpuEdge[1]);
    const uint64_t gpuSlices = FRAME(arrnGpuEdge[2]);
    const uint64_t nRows   = FRAME(nRows);
	const uint64_t nCols   = FRAME(nCols);
    const uint64_t nSlices = FRAME(nSlices);


    const int szdb = sizeof(double);
    
    int nId              = FRAME(nId);
    int *arrnGpuEdgeChunk             = FRAME(arrnGpuEdgeChunk);
    
    cudaError err0,err1,err2,err3,err4;
	int gpuTile_x = FRAME(arrnGpuTile[0]); //16
	int gpuTile_y = FRAME(arrnGpuTile[1]); //16
    int gpuTile_z = FRAME(arrnGpuTile[2]); //16


    int gnStream = FRAME(gnStream);
	
	int gpuTileCols		= arrnGpuEdgeChunk[0];
	int gpuTileRows		= arrnGpuEdgeChunk[1];
	int gpuTileSlices	= arrnGpuEdgeChunk[2];
     
	int gpuTileCols2	= ((gpuTileCols  +2)>gpuCols  )? gpuCols  :(gpuTileCols  +2);
	int gpuTileRows2    = ((gpuTileRows  +2)>gpuRows  )? gpuRows  :(gpuTileRows  +2);   
	int gpuTileSlices2	= ((gpuTileSlices+2)>gpuSlices)? gpuSlices:(gpuTileSlices+2);	
   
    //cudaExtent{width,height,depth} width=gpuTileCols,height=gpuTileRows,depth=gpuTileSlices
    cudaExtent extent = make_cudaExtent(gpuTileCols2*szdb,gpuTileRows2,gpuTileSlices2*gnStream);
    cudaPitchedPtr devPitchedPtr;
    err0=cudaMalloc3D(&devPitchedPtr, extent);
    
#ifdef CUDA_DARTS_DEBUG
    std::cout<<"extend: width: "<<extent.width<<",heigh: "<<extent.height<<", depth: "<<extent.depth<<std::endl;
    std::cout<<"devPitchedPtr: pointer: "<<devPitchedPtr.ptr<<",pitch: "<<devPitchedPtr.pitch<<",xsize: "<<devPitchedPtr.xsize<<", ysize: "<<devPitchedPtr.ysize<<std::endl;
#endif
#ifdef CUDA_ERROR_CHECKING
    gpuErrchk(err0);
#endif

	size_t d_xpitch = devPitchedPtr.pitch/szdb ;
    size_t d_ypitch = extent.height;
    size_t d_zpitch = extent.depth/gnStream;
    int64_t dev_size = d_xpitch * d_ypitch * d_zpitch;

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"dev size unit: "<<dev_size<<std::endl;
    std::cout<<"d_xpitch:"<<d_xpitch<<std::endl;
    std::cout<<"d_ypitch:"<<d_ypitch<<std::endl;
    std::cout<<"d_zpitch:"<<d_zpitch<<std::endl;
#endif

    cudaMemcpy3DParms p ={0};
    p.srcPtr    = make_cudaPitchedPtr((void *)h_src,gpuCols*szdb,gpuCols*szdb,gpuRows);
    p.dstPtr    = devPitchedPtr;
    p.extent.width  = gpuTileCols2*szdb;
    p.extent.height = gpuTileRows2;
    p.extent.depth  = gpuTileSlices2;
    p.kind   = cudaMemcpyHostToDevice;
	
	cudaMemcpy3DParms pdh ={0};
    pdh.srcPtr    = devPitchedPtr;
    pdh.dstPtr    = make_cudaPitchedPtr((void *)h_dst,gpuCols*szdb,gpuCols*szdb,gpuRows);
    pdh.extent.width  = gpuTileCols2*szdb;
    pdh.extent.height = gpuTileRows2;
    pdh.extent.depth  = gpuTileSlices2;
    pdh.kind   = cudaMemcpyDeviceToHost;
   
#ifdef CUDA_DARTS_DEBUG
    std::cout<<"p.srcPtr: pointer: "<<p.srcPtr.ptr<<",pitch: "<<p.srcPtr.pitch<<",xsize: "<<p.srcPtr.xsize<<", ysize: "<<p.srcPtr.ysize<<std::endl;
    std::cout<<"p.dstPtr: pointer: "<<p.dstPtr.ptr<<",pitch: "<<p.dstPtr.pitch<<",xsize: "<<p.dstPtr.xsize<<", ysize: "<<p.dstPtr.ysize<<std::endl;
    std::cout<<"pdh.srcPtr: pointer: "<<pdh.srcPtr.ptr<<",pitch: "<<pdh.srcPtr.pitch<<",xsize: "<<pdh.srcPtr.xsize<<", ysize: "<<pdh.srcPtr.ysize<<std::endl;
    std::cout<<"pdh.dstPtr: pointer: "<<pdh.dstPtr.ptr<<",pitch: "<<pdh.dstPtr.pitch<<",xsize: "<<pdh.dstPtr.xsize<<", ysize: "<<pdh.dstPtr.ysize<<std::endl;
#endif

    
    int gpuBlockDimx = gpuTileCols>gpuTile_x?gpuTile_x:gpuTileCols;
	int gpuBlockDimy = gpuTileRows>gpuTile_y?gpuTile_y:gpuTileRows;
	int gpuBlockDimz = 1;
	int gpuGridDimx2 = std::ceil(1.0*gpuTileCols2/gpuBlockDimx);
	int gpuGridDimy2 = std::ceil(1.0*gpuTileRows2/gpuBlockDimy);
	int gpuGridDimz2 = std::ceil(1.0*gpuTileSlices2/gpuTile_z);


    int64_t srows_size = gpuTileCols2*gpuGridDimy2*gpuTileSlices2*2;
	int64_t scols_size = gpuTileRows2*gpuGridDimx2*gpuTileSlices2*2;
    int64_t sslices_size = gpuTileCols2*gpuTileRows2*gpuGridDimz2*2;


    int64_t t_d_size_sharedCols     = gnStream*szdb*scols_size;    	
    int64_t t_d_size_sharedRows     = gnStream*szdb*srows_size;    	
    int64_t t_d_size_sharedSlices   = gnStream*szdb*sslices_size;	

    //err1 = cudaMalloc( (void **) &d_dst, d_size*gnStream);
	err2 = cudaMalloc(  &d_sharedCols, (t_d_size_sharedCols));
#ifdef CUDA_ERROR_CHECKING
	gpuErrchk(err2);
#endif
    err3 = cudaMalloc( (void **) &d_sharedRows, (t_d_size_sharedRows));
#ifdef CUDA_ERROR_CHECKING
	gpuErrchk(err3);
#endif
    err4 = cudaMalloc( (void **) &d_sharedSlices, (t_d_size_sharedSlices));
#ifdef CUDA_ERROR_CHECKING
	gpuErrchk(err4);
#endif

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"device dst   size :"<< (dev_size*szdb*gnStream)<<std::endl;
    std::cout<<"sharedCols   size :"<< (t_d_size_sharedCols)<<std::endl;
    std::cout<<"sharedRows   size :"<< (t_d_size_sharedRows)<<std::endl;
    std::cout<<"sharedSlices size :"<< (t_d_size_sharedSlices)<<std::endl;
#endif

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"device dst   address:"<<p.dstPtr.ptr   <<" to "<<p.dstPtr.ptr   + (dev_size*szdb*gnStream)   <<std::endl;
    std::cout<<"sharedCols   address:"<<d_sharedCols   <<" to "<<d_sharedCols   + t_d_size_sharedCols/szdb   <<std::endl;
    std::cout<<"sharedRows   address:"<<d_sharedRows   <<" to "<<d_sharedRows   + t_d_size_sharedRows/szdb   <<std::endl;
    std::cout<<"sharedSlices address:"<<d_sharedSlices <<" to "<< d_sharedSlices+ t_d_size_sharedSlices/szdb   <<std::endl;
#endif

	size_t s_xpitch = gpuTileCols2;
    size_t s_ypitch = gpuTileRows2;
    size_t s_zpitch = gpuTileSlices2;

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"s_xpitch:"<<s_xpitch<<std::endl;
    std::cout<<"s_ypitch:"<<s_ypitch<<std::endl;
    std::cout<<"s_zpitch:"<<s_zpitch<<std::endl;
#endif


#ifdef CUDA_DARTS_DEBUG
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
	
	std::cout<<std::setprecision(18)<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: gpu memory total: "<<gpu_mem_total_t<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37: gpu memory available: "<<gpu_mem_avail_t<<std::endl;
	std::cout<<"GpuKernelPureGpuWithStreams37 : require memory size:"<< (dev_size*szdb*gnStream)  + (t_d_size_sharedCols + t_d_size_sharedRows+t_d_size_sharedSlices)<<std::endl;
#endif


    int nGPU     = FRAME(nGPU);


    int numThreads= gpuTile_x*gpuTile_y;
    
    int blockDimx_slices ;
    int blockDimy_slices = 1;
	int blockDimz_slices = 1;
	int gridDimx_slices ;
	int gridDimy_slices ;
    int gridDimz_slices ;

	int blockDimx_rows  ;
	int blockDimy_rows = 1;
	int blockDimz_rows = 1;
	int gridDimx_rows ;
	int gridDimy_rows ;
	int gridDimz_rows ;

	int blockDimx_cols = 1 ;
	int blockDimy_cols; 
    int blockDimz_cols = 1 ;
	int gridDimx_cols ;
	int gridDimy_cols ;
	int gridDimz_cols ;
	
    dim3 dimGrid; 
    dim3 dimBlock;	
	dim3 dimGrid_slices; 
	dim3 dimBlock_slices;
	
    dim3 dimGrid_rows; 
	dim3 dimBlock_rows;

	dim3 dimGrid_cols; 
	dim3 dimBlock_cols;
   


	int nb_x = std::ceil(1.0*gpuCols/gpuTileCols);//x: how many gpu Tile CRS [x] 
	int nb_y = std::ceil(1.0*gpuRows/gpuTileRows);//y: how many gpu Tile CRS [y] 
	int nb_z = std::ceil(1.0*gpuSlices/gpuTileSlices);//z: how many gpu Tile CRS [z] 

	nb_x = (gpuTileCols2  < gpuTileCols +2   )? (1 ):(nb_x );  
	nb_y = (gpuTileRows2  < gpuTileRows + 2  )? (1 ):(nb_y );  
	nb_z = (gpuTileSlices2< gpuTileSlices + 2)? (1 ):(nb_z );


#ifdef CUDA_DARTS_DEBUG

    std::cout<<"gpuCols      : "<< gpuCols   <<std::endl;
    std::cout<<"gpuRows      : "<< gpuRows   <<std::endl;
    std::cout<<"gpuSlices    : "<< gpuSlices <<std::endl;
    
    std::cout<<"gpuTileCols  : "<< gpuTileCols   <<std::endl;
    std::cout<<"gpuTileRows  : "<< gpuTileRows   <<std::endl;
    std::cout<<"gpuTileSlices: "<< gpuTileSlices <<std::endl;
 
    std::cout<<"gpuTileCols2  : "<< gpuTileCols2   <<std::endl;
    std::cout<<"gpuTileRows2  : "<< gpuTileRows2   <<std::endl;
    std::cout<<"gpuTileSlices2: "<< gpuTileSlices2 <<std::endl;

    std::cout<<"gpuBlockDimx : "<< gpuBlockDimx<<std::endl;
    std::cout<<"gpuBlockDimy : "<< gpuBlockDimy<<std::endl;
    std::cout<<"gpuBlockDimz : "<< gpuBlockDimz<<std::endl;

    std::cout<<"gpuGridDimx2 : "<< gpuGridDimx2<<std::endl;
    std::cout<<"gpuGridDimy2 : "<< gpuGridDimy2<<std::endl;
    std::cout<<"gpuGridDimz2 : "<< gpuGridDimz2<<std::endl;

    std::cout<<"nb_x         : "<< nb_x <<std::endl;
    std::cout<<"nb_y         : "<< nb_y <<std::endl;
    std::cout<<"nb_z         : "<< nb_z <<std::endl;
#endif
    
    int nStream  = FRAME(nStream);
    int vnStream = nb_x*nb_y*nb_z;
    int tnStream = vnStream*nStream; 

    cudaStream_t *stream ;
	stream = new cudaStream_t[tnStream];
    for(int i=0;i<tnStream;++i){
        cudaStreamCreate(&stream[i]);
    }


    cudaEvent_t *cuEvent = new cudaEvent_t[tnStream];
    cudaError *err = new cudaError[tnStream];
    for (int i=0;i<tnStream;++i){
        err[i] = cudaEventCreate(&cuEvent[i]);
#ifdef CUDA_ERROR_CHECKING
		gpuErrchk(err[i]);
#endif
    }
    d_dst = devPitchedPtr.ptr;	
	FRAME(d_dst) = d_dst;
    //double pos0 = FRAME(gpuPos); 
    //int vnStream = nGPU*gnStream;

    struct_streams_par *streams_par = new struct_streams_par [tnStream];

    size_t ts = FRAME(ts);
    while(ts-- >0){

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"time step: "<<ts <<std::endl;
        
#endif       
        p.srcPtr.ptr = h_src;
        pdh.dstPtr.ptr = h_dst;
        for(int i=0;i<vnStream;++i){
        //for(int i=0;i<3;i=i+1){
        //    int i=2;{
#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
				int posz = i/(nb_x*nb_y);
				int posy = i/nb_x-posz*nb_y;
                int posx = i%nb_x;
                
                int ps = i%gnStream;

                int nRowsChunk		= (((posy+1)*gpuTileRows)	>=gpuRows)	?(gpuRows-posy*gpuTileRows)	 :(gpuTileRows2) ;
				int nColsChunk		= (((posx+1)*gpuTileCols)	>=gpuCols)	?(gpuCols-posx*gpuTileCols)	 :(gpuTileCols2) ;
				int nSlicesChunk	= (((posz+1)*gpuTileSlices)	>=gpuSlices)?(gpuSlices-posz*gpuTileSlices):(gpuTileSlices2) ;
                
#ifdef CUDA_DARTS_DEBUG

                std::cout<<"stream #: "<<i<<" posx: "<<posx<<std::endl;
                std::cout<<"stream #: "<<i<<" posy: "<<posy<<std::endl;
                std::cout<<"stream #: "<<i<<" posz: "<<posz<<std::endl;
                std::cout<<"stream #: "<<i<<" nColsChunk   = "<< nColsChunk	 <<std::endl;
                std::cout<<"stream #: "<<i<<" nRowsChunk   = "<< nRowsChunk	 <<std::endl;
                std::cout<<"stream #: "<<i<<" nSlicesChunk = "<< nSlicesChunk   <<std::endl;
#endif
                
                blockDimx_slices = (nColsChunk>numThreads)?numThreads:nColsChunk; 
	            blockDimx_rows   = (nColsChunk>numThreads)?numThreads:nColsChunk;
	            blockDimy_cols   = (nRowsChunk>numThreads)?numThreads:nRowsChunk;
               

                //gpuBlockDimx = (nColsChunk-2  )>gpuTile_x?gpuTile_x:(nColsChunk-2);
				//gpuBlockDimy = (nRowsChunk-2  )>gpuTile_y?gpuTile_y:(nRowsChunk-2);
				//gpuBlockDimz = 1;
               
                gpuBlockDimx = gpuTile_x;
				gpuBlockDimy = gpuTile_y;
				gpuBlockDimz = 1;

                gpuGridDimx2 = std::ceil(1.0*nColsChunk/gpuBlockDimx);
                gpuGridDimy2 = std::ceil(1.0*nRowsChunk/gpuBlockDimy);
                gpuGridDimz2 = std::ceil(1.0*nSlicesChunk/gpuTile_z);
	
	            gridDimx_slices = std::ceil(1.0*nColsChunk/blockDimx_slices);
	            gridDimy_slices = gpuGridDimy2;  
                gridDimz_slices = gpuGridDimz2;

	            gridDimx_rows = std::ceil(1.0*nColsChunk/blockDimx_rows);
	            gridDimy_rows = gpuGridDimy2; 
	            gridDimz_rows = gpuGridDimz2; 

	            gridDimx_cols = gpuGridDimx2; 
	            gridDimy_cols = std::ceil(1.0*nRowsChunk/blockDimy_cols);
	            gridDimz_cols = gpuGridDimz2; 
                

#ifdef CUDA_DARTS_DEBUG

                std::cout<<"stream #: "<<i<<" blockDimx_slices:"<<blockDimx_slices<<std::endl;
	            std::cout<<"stream #: "<<i<<" blockDimy_slices:"<<blockDimy_slices<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimx_slices:"<<gridDimx_slices<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimy_slices:"<<gridDimy_slices<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimz_slices:"<<gridDimz_slices<<std::endl;
                std::cout<<"stream #: "<<i<<" blockDimx_rows:"<<blockDimx_rows<<std::endl;
	            std::cout<<"stream #: "<<i<<" blockDimy_rows:"<<blockDimy_rows<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimx_rows:"<<gridDimx_rows<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimy_rows:"<<gridDimy_rows<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimz_rows:"<<gridDimz_rows<<std::endl;
	            
                std::cout<<"stream #: "<<i<<" blockDimx_cols:"<<blockDimx_cols<<std::endl;
	            std::cout<<"stream #: "<<i<<" blockDimy_cols:"<<blockDimy_cols<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimx_cols:"<<gridDimx_cols<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimy_cols:"<<gridDimy_cols<<std::endl;
	            std::cout<<"stream #: "<<i<<" grimDimz_cols:"<<gridDimz_cols<<std::endl;

#endif
                
                int idxSM=i*nStream;
				int nSMM1=nStream-1;
				int idxSMM1 = idxSM+nSMM1;
				if((i/gnStream) !=0){
					int pre = (i-gnStream)*nStream;
					err[pre]=cudaStreamWaitEvent(stream[idxSM],cuEvent[pre+nSMM1],0);
#ifdef CUDA_DARTS_DEBUG
                    std::cout<<"i: "<<i<<" ,event: "<<idxSM<<" waiting for event: "<<pre+nSMM1<<std::endl;
#endif
#ifdef CUDA_ERROR_CHECKING
					gpuErrchk(err[pre]);
#endif
                }           
                int64_t d_posz = ps*gpuTileSlices2;
                
                p.srcPos = make_cudaPos(posx*gpuTileCols*szdb,posy*gpuTileRows,posz*gpuTileSlices);
                p.dstPos = make_cudaPos(0,0,d_posz);

                p.extent.width  = nColsChunk*szdb;
                p.extent.height = nRowsChunk;
                p.extent.depth  = nSlicesChunk;

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"GpuKernelWithStream multiple streams: "<<idxSM <<" memory copy begin!"<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" posx: "<<posx<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" posy: "<<posy<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" posz: "<<posz<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.x: "<<p.srcPos.x<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.y: "<<p.srcPos.y<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.z: "<<p.srcPos.z<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.x: "<<p.dstPos.x<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.y: "<<p.dstPos.y<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.z: "<<p.dstPos.z<<std::endl;

                std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.ptr: "  <<p.srcPtr.ptr  <<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.pitch: "<<p.srcPtr.pitch<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.xsize: "<<p.srcPtr.xsize<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.ysize: "<<p.srcPtr.ysize<<std::endl;

                std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.ptr: "  <<p.dstPtr.ptr  <<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.pitch: "<<p.dstPtr.pitch<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.xsize: "<<p.dstPtr.xsize<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.ysize: "<<p.dstPtr.ysize<<std::endl;

                std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.width: " <<p.extent.width<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.height: "<<p.extent.height<<std::endl;
                std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.depth: " <<p.extent.depth<<std::endl;
#endif
                
                pdh.srcPos = make_cudaPos(1*szdb,1,d_posz+1);
                pdh.dstPos = make_cudaPos((posx*gpuTileCols+1)*szdb,posy*gpuTileRows+1,posz*gpuTileSlices+1);


				int nRowsChunkM2	= (((posy+1)*gpuTileRows)	>=gpuRows)	?(gpuRows-posy*gpuTileRows-2)	 :(gpuTileRows) ;
				int nColsChunkM2	= (((posx+1)*gpuTileCols)	>=gpuCols)	?(gpuCols-posx*gpuTileCols-2)	 :(gpuTileCols) ;
				int nSlicesChunkM2	= (((posz+1)*gpuTileSlices)	>=gpuSlices)?(gpuSlices-posz*gpuTileSlices-2):(gpuTileSlices) ;

                pdh.extent.width  = nColsChunkM2*szdb;
                pdh.extent.height = nRowsChunkM2;
                pdh.extent.depth  = nSlicesChunkM2;

                streams_par[idxSM].htodPtr = &p;
                streams_par[idxSM].dtohPtr = &pdh; 
                streams_par[idxSM].devDstPos    = ps*dev_size;
                streams_par[idxSM].devSRowsPos  = ps*srows_size;
                streams_par[idxSM].devSColsPos  = ps*scols_size;
                streams_par[idxSM].devSSlicesPos= ps*sslices_size; 	 
   
	            streams_par[idxSM].dimBlock_slices = dim3 (blockDimx_slices,blockDimy_slices,blockDimz_slices);
	            streams_par[idxSM].dimBlock_rows   = dim3 (blockDimx_rows,blockDimy_rows,blockDimz_rows);
	            streams_par[idxSM].dimBlock_cols   = dim3 (blockDimx_cols,blockDimy_cols,blockDimz_cols);
	            streams_par[idxSM].dimBlock        = dim3 (gpuBlockDimx,gpuBlockDimy,gpuBlockDimz);


                streams_par[idxSM].dimGrid_slices  = dim3 (gridDimx_slices,gridDimy_slices,gridDimz_slices);
                streams_par[idxSM].dimGrid_rows    = dim3 (gridDimx_rows,gridDimy_rows,gridDimz_rows);
	            streams_par[idxSM].dimGrid_cols    = dim3 (gridDimx_cols,gridDimy_cols,gridDimz_cols);
                streams_par[idxSM].dimGrid         = dim3 (gpuGridDimx2, gpuGridDimy2,gpuGridDimz2);
                streams_par[idxSM].nRowsChunk  = nRowsChunk  ;	
                streams_par[idxSM].nColsChunk  = nColsChunk  ;
                streams_par[idxSM].nSlicesChunk= nSlicesChunk; 

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" device dst   address from: "<<p.dstPtr.ptr  + streams_par[idxSM].devDstPos    <<" to: "<<p.dstPtr.ptr  + streams_par[idxSM].devDstPos    + dev_size      <<std::endl;
                std::cout<<"stream: "<<i<<" sharedSlices address from: "<<d_sharedSlices+ streams_par[idxSM].devSSlicesPos  <<" to: "<<d_sharedSlices+ streams_par[idxSM].devSSlicesPos  + sslices_size    <<std::endl;
                std::cout<<"stream: "<<i<<" sharedRows   address from: "<<d_sharedRows  + streams_par[idxSM].devSRowsPos  <<" to: "<<d_sharedRows  + streams_par[idxSM].devSRowsPos  + srows_size    <<std::endl;
                std::cout<<"stream: "<<i<<" sharedCols   address from: "<<d_sharedCols  + streams_par[idxSM].devSColsPos<<" to: "<<d_sharedCols  + streams_par[idxSM].devSColsPos+ scols_size <<std::endl;
#endif

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream #: "<<i<<" posx: "<<posx<<std::endl;
                std::cout<<"stream #: "<<i<<" posy: "<<posy<<std::endl;
                std::cout<<"stream #: "<<i<<" posz: "<<posz<<std::endl;
                std::cout<<"stream #: "<<i<<" nColsChunkM2   = "<< nColsChunkM2	 <<std::endl;
                std::cout<<"stream #: "<<i<<" nRowsChunkM2   = "<< nRowsChunkM2	 <<std::endl;
                std::cout<<"stream #: "<<i<<" nSlicesChunkM2 = "<< nSlicesChunkM2   <<std::endl;
#endif


#ifdef CUDA_DARTS_DEBUG
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
   
#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
                err0=cudaMemcpy3DAsync(streams_par[idxSM].htodPtr,stream[idxSM]);
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err0);
#endif
 

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif

#ifdef CUDA_DARTS_DEBUG
                {
                    
#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
			        err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                    dim3 dimGrid_t (1,1,1);
	                dim3 dimBlock_t(10,1,1);
                    
                    double *sptr    = p.srcPtr.ptr;
                    int spitch      = p.srcPtr.pitch/szdb;
                    int sposx       = p.srcPos.x/szdb; 
                    int sposy       = p.srcPos.y; 
                    int sposz       = p.srcPos.z;
                    int sxsize      = p.srcPtr.xsize;
                    int sysize      = p.srcPtr.ysize;
                    
                    double *dptr    = p.dstPtr.ptr;
                    int dpitch      = p.dstPtr.pitch/szdb;
                    int dposx       = p.dstPos.x; 
                    int dposy       = p.dstPos.y; 
                    int dposz       = p.dstPos.z;

                    int dxsize      = p.dstPtr.xsize;
                    int dysize      = p.dstPtr.ysize;

                    int stride= 0;
                    std::cout<<"Async cp p: i: "<<i<<std::endl;
                    std::cout<<"Async cp p: p.srcPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,sptr,sposz*spitch*sysize +sposx+1*spitch*sysize+sposy*spitch+2*spitch+stride,1,1);

#ifdef CUDA_ERROR_CHECKING
                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
#endif
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                    std::cout<<"Async cp p: p.dstPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,dptr,dposz*dysize*dpitch+dposx + 1*dpitch*dysize+dposy*dpitch+2*dpitch+stride,1,1);

#ifdef CUDA_ERROR_CHECKING
                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
#endif
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
                }
#endif

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"GpuKernelWithStream multiple streams: "<<idxSM <<" memory copy end!"<<std::endl;
#endif
                
                err[idxSM]=cudaEventRecord(cuEvent[idxSM],stream[idxSM]);
                for(int k=1;k<nSMM1;++k){
                    err[idxSM+k]=cudaStreamWaitEvent(stream[idxSM+k],cuEvent[idxSM],0);
                }
                
#ifdef CUDA_ERROR_CHECKING
                for(int k=1;k<nSMM1;++k){
					gpuErrchk(err[idxSM+k]);
                }
#endif

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedCols   address from: "<<d_sharedCols  + streams_par[idxSM].devSColsPos<<" to: "<<d_sharedCols  + streams_par[idxSM].devSColsPos+ scols_size <<std::endl;
#endif
                gpu_kernel37_cp_cols_stream_p(stream[idxSM+1],streams_par[idxSM].dimGrid_cols,streams_par[idxSM].dimBlock_cols,d_dst+streams_par[idxSM].devDstPos,d_sharedRows+streams_par[idxSM].devSRowsPos, d_sharedCols+ streams_par[idxSM].devSColsPos, d_sharedSlices+ streams_par[idxSM].devSSlicesPos,d_xpitch,d_ypitch, d_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+1]=cudaEventRecord(cuEvent[idxSM+1],stream[idxSM+1]);

#ifdef CUDA_DARTS_DEBUG

                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
                
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedRows   address from: "<<d_sharedRows  + streams_par[idxSM].devSRowsPos  <<" to: "<<d_sharedRows  + streams_par[idxSM].devSRowsPos  + srows_size    <<std::endl;
#endif
                
                gpu_kernel37_cp_rows_stream_p(stream[idxSM+2],streams_par[idxSM].dimGrid_rows,streams_par[idxSM].dimBlock_rows,d_dst+streams_par[idxSM].devDstPos,d_sharedRows+streams_par[idxSM].devSRowsPos, d_sharedCols+ streams_par[idxSM].devSColsPos, d_sharedSlices+ streams_par[idxSM].devSSlicesPos,d_xpitch,d_ypitch, d_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+2]=cudaEventRecord(cuEvent[idxSM+2],stream[idxSM+2]);
                
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
                
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedSlices address from: "<<d_sharedSlices+ streams_par[idxSM].devSSlicesPos  <<" to: "<<d_sharedSlices+ streams_par[idxSM].devSSlicesPos  + sslices_size    <<std::endl;
#endif
                gpu_kernel37_cp_slices_stream_p(stream[idxSM+3],streams_par[idxSM].dimGrid_slices,streams_par[idxSM].dimBlock_slices,d_dst+streams_par[idxSM].devDstPos,d_sharedRows+streams_par[idxSM].devSRowsPos, d_sharedCols+ streams_par[idxSM].devSColsPos, d_sharedSlices+ streams_par[idxSM].devSSlicesPos,d_xpitch,d_ypitch, d_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+3]=cudaEventRecord(cuEvent[idxSM+3],stream[idxSM+3]);
                
#ifdef CUDA_DARTS_DEBUG
                
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
        
#ifdef CUDA_ERROR_CHECKING
                for(int k=1;k<nSMM1;++k){
					gpuErrchk(err[idxSM+k]);
                }
#endif
                for(int k=1;k<nSMM1;++k){
                    err[idxSM+k]=cudaStreamWaitEvent(stream[idxSMM1],cuEvent[idxSM+k],0);
                }
                
#ifdef CUDA_ERROR_CHECKING
                for(int k=1;k<nSMM1;++k){
					gpuErrchk(err[idxSM+k]);
				}
#endif
      
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
  
                gpu_kernel37_stream_p(stream[idxSMM1],streams_par[idxSM].dimGrid,streams_par[idxSM].dimBlock,d_dst+streams_par[idxSM].devDstPos,d_sharedRows+streams_par[idxSM].devSRowsPos, d_sharedCols+ streams_par[idxSM].devSColsPos, d_sharedSlices+ streams_par[idxSM].devSSlicesPos,d_xpitch,d_ypitch, d_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
        
#ifdef CUDA_DARTS_DEBUG
                
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif

#ifdef CUDA_DARTS_DEBUG
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.x: "        <<pdh.srcPos.x<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.y: "        <<pdh.srcPos.y<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.z: "        <<pdh.srcPos.z<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.x: "        <<pdh.dstPos.x<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.y: "        <<pdh.dstPos.y<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.z: "        <<pdh.dstPos.z<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.ptr: "      <<pdh.srcPtr.ptr  <<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.pitch: "    <<pdh.srcPtr.pitch<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.xsize: "    <<pdh.srcPtr.xsize<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.ysize: "    <<pdh.srcPtr.ysize<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.ptr: "      <<pdh.dstPtr.ptr  <<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.pitch: "    <<pdh.dstPtr.pitch<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.xsize: "    <<pdh.dstPtr.xsize<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.ysize: "    <<pdh.dstPtr.ysize<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.width: "    <<pdh.extent.width<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.height: "   <<pdh.extent.height<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.depth: "    <<pdh.extent.depth<<std::endl;

#endif



#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
                err0=cudaMemcpy3DAsync(streams_par[idxSM].dtohPtr,stream[idxSMM1]);
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err0);
#endif

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
				err[idxSMM1]=cudaEventRecord(cuEvent[idxSMM1],stream[idxSMM1]);
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err[idxSMM1]);
#endif

#ifdef CUDA_DARTS_DEBUG
                {

#ifdef CUDA_DARTS_DEBUG
                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
#endif
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                    dim3 dimGrid_t (1,1,1);
	                dim3 dimBlock_t(10,1,1);

                    double *sptr    = pdh.srcPtr.ptr;
                    int spitch      = pdh.srcPtr.pitch/szdb;
                    int sposx       = pdh.srcPos.x/szdb; 
                    int sposy       = pdh.srcPos.y; 
                    int sposz       = pdh.srcPos.z; 
                    int sxsize      = pdh.srcPtr.xsize;
                    int sysize      = pdh.srcPtr.ysize;
                    
                    double *dptr    = pdh.dstPtr.ptr;
                    int dpitch      = pdh.dstPtr.pitch/szdb;
                    int dposx       = pdh.dstPos.x/szdb; 
                    int dposy       = pdh.dstPos.y; 
                    int dposz       = pdh.dstPos.z;
                    int dxsize      = pdh.dstPtr.xsize;
                    int dysize      = pdh.dstPtr.ysize;

                    
                    int stride= 0;
                    std::cout<<"Async cp pdh: i: "<<i<<std::endl;
                    std::cout<<"Async cp pdh: pdh.srcPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,sptr,sposz*spitch*sysize +sposx+0*spitch*sysize+sposy*spitch+1*spitch+stride,1,1);

                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
                    std::cout<<"Async cp pdh: pdh.dstPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,dptr,dposz*dpitch*dysize +dposx+0*dpitch*dysize+dposy*dpitch+1*dpitch+stride,1,1);

                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                }
#endif


#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
        }
		err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
		gpuErrchk(err1);
#endif	
		
		SWAP_PTR(&h_dst ,&h_src);
        //SWAP_PTR(&p.srcPtr.ptr, &pdh.dstPtr.ptr);
		
	}

	err0 = cudaGetLastError();
#ifdef CUDA_ERROR_CHECKING
	gpuErrchk(err0);
#endif	

    //destroy cudaEvent
    for(int i=0;i<tnStream;++i){
        err[i]=cudaEventDestroy(cuEvent[i]);
        
#ifdef CUDA_ERROR_CHECKING
		for(int i=0;i<tnStream;++i){
			gpuErrchk(err[i]);
		}
#endif
    }
    delete [] cuEvent;
    delete [] err;
    delete [] streams_par;

	for(int i=0;i<tnStream;++i){
		cudaStreamDestroy(stream[i]);
	}
    delete [] stream;


#ifdef VERIFICATION
    if(FRAME(ts)%2==0){
		SWAP_PTR(&h_dst ,&h_src);
	}        
#endif

//#ifdef CUDA_DARTS_DEBUG
//
//
//#endif
//
	//err1 = cudaFree(d_dst);
	err1 = cudaFree(d_dst);
	err2 = cudaFree(d_sharedCols);
	err3 = cudaFree(d_sharedRows);
    err4 = cudaFree(d_sharedSlices);
	
#ifdef CUDA_ERROR_CHECKING
	gpuErrchk(err1);
	gpuErrchk(err2);
	gpuErrchk(err3);
	gpuErrchk(err4);
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

    int nId = FRAME(nId);
    int nIdPlus1 = R(nId+1);
    int nIdMinus1 = R(nId-1);

    double	*h_src	= FRAME(Initial);
	double	*h_dst	= FRAME(New);

    const int szdb = FRAME(szdb);
    
    int *arrnEdge = FRAME(arrnEdge);
    int *arrnGpuTile = FRAME(arrnGpuTile);

	const int nCols   = arrnEdge[0];
    const int nRows   = arrnEdge[1];
    const int nSlices = arrnEdge[2];


    double *dptr_dev                = FRAME(dptr_dev              );
	double *dptr_dev_sharedCols     = FRAME(dptr_dev_sharedCols   );
	double *dptr_dev_sharedRows     = FRAME(dptr_dev_sharedRows   );
	double *dptr_dev_sharedSlices   = FRAME(dptr_dev_sharedSlices );

    int64_t srows_size   = FRAME(srows_size   ); 
	int64_t scols_size   = FRAME(scols_size   ); 
    int64_t sslices_size = FRAME(sslices_size ); 

    //int64_t dev_sharedCols_size     = FRAME(dev_sharedCols_size  );  	
    //int64_t dev_sharedRows_size     = FRAME(dev_sharedRows_size  );  	
    //int64_t dev_sharedSlices_size   = FRAME(dev_sharedSlices_size);	


#ifdef CUDA_DARTS_DEBUG
    std::cout<<"h_src addr: "<<h_src<<std::endl;
    std::cout<<"h_dst addr: "<<h_dst<<std::endl;
    
    std::cout<<"device dst   address:"<<dptr_dev             <<std::endl;
    std::cout<<"sharedCols   address:"<<dptr_dev_sharedCols  <<std::endl;
    std::cout<<"sharedRows   address:"<<dptr_dev_sharedRows  <<std::endl;
    std::cout<<"sharedSlices address:"<<dptr_dev_sharedSlices<<std::endl;
#endif

    //========================= cudaMemcpy3DParams ==========================//
   
    //cudaMemcpy3DParms p ={0};
    //p.srcPtr    = make_cudaPitchedPtr((void *)h_src,gpuCols*szdb,gpuCols*szdb,gpuRows);
    //p.dstPtr    = devPitchedPtr;
    //p.extent.width  = gpuTileCols2*szdb;
    //p.extent.height = gpuTileRows2;
    //p.extent.depth  = gpuTileSlices2;
    //p.kind   = cudaMemcpyHostToDevice;
    
    cudaMemcpy3DParms phd ={0};
    phd.srcPtr = make_cudaPitchedPtr((void *)h_src,nCols*szdb,nCols*szdb,nRows);
    phd.dstPtr = FRAME(devPitchedPtr);
    phd.kind   = cudaMemcpyHostToDevice;
	
	cudaMemcpy3DParms pdh ={0};
    pdh.srcPtr  = FRAME(devPitchedPtr);
    pdh.dstPtr  = make_cudaPitchedPtr((void *)h_dst,nCols*szdb,nCols*szdb,nRows);
;
    pdh.kind   = cudaMemcpyDeviceToHost;

#ifdef CUDA_DARTS_DEBUG
    std::cout<<"phd.srcPtr: pointer: "<<phd.srcPtr.ptr<<",pitch: "<<phd.srcPtr.pitch<<",xsize: "<<phd.srcPtr.xsize<<", ysize: "<<phd.srcPtr.ysize<<std::endl;
    std::cout<<"phd.dstPtr: pointer: "<<phd.dstPtr.ptr<<",pitch: "<<phd.dstPtr.pitch<<",xsize: "<<phd.dstPtr.xsize<<", ysize: "<<phd.dstPtr.ysize<<std::endl;
    std::cout<<"pdh.srcPtr: pointer: "<<pdh.srcPtr.ptr<<",pitch: "<<pdh.srcPtr.pitch<<",xsize: "<<pdh.srcPtr.xsize<<", ysize: "<<pdh.srcPtr.ysize<<std::endl;
    std::cout<<"pdh.dstPtr: pointer: "<<pdh.dstPtr.ptr<<",pitch: "<<pdh.dstPtr.pitch<<",xsize: "<<pdh.dstPtr.xsize<<", ysize: "<<pdh.dstPtr.ysize<<std::endl;
#endif

	size_t dev_xpitch = FRAME( dev_xpitch ); 
    size_t dev_ypitch = FRAME( dev_ypitch ); 
    size_t dev_zpitch = FRAME( dev_zpitch ); 
    int64_t dev_size  = FRAME( dev_size);


    //========================= define variables ==========================//

    int *arrnGpuEdgeChunk = FRAME(arrnGpuEdgeChunk);
    int *arrnGpuEdgeChunkAlloc = FRAME(arrnGpuEdgeChunkAlloc);
    int *arrnCpuEdgeMin     = FRAME(arrnCpuEdgeMin);
    int *arrnGpuEdgeMin     = FRAME(arrnGpuEdgeMin);

    int nStream = FRAME(nStream);
    int gnStream = FRAME(gnStream);

	int gpuTile_x = arrnGpuTile[0]; //16
	int gpuTile_y = arrnGpuTile[1]; //16
    int gpuTile_z = arrnGpuTile[2]; //16

    const int gpuColsAlloc     = arrnGpuEdgeChunkAlloc[0];
    const int gpuRowsAlloc     = arrnGpuEdgeChunkAlloc[1];
    const int gpuSlicesAlloc   = arrnGpuEdgeChunkAlloc[2];


    cudaError err0,err1,err2,err3,err4;


    int numThreads= gpuTile_x*gpuTile_y;
    

    int arrnNumThreads[DIM];
    int arrnSSBlock[DIM][DIM];//(cols,rows,slices)(x,y,z)
    int arrnSSGrid[DIM][DIM];
    int arrnGpuBlock[DIM];
    int arrnGpuGrid[DIM];

    dim3 dimGrid; 
    dim3 dimBlock;	
	dim3 dimSSGrid[DIM]; 
	dim3 dimSSBlock[DIM];
    
    setarrnValue(arrnNumThreads,numThreads);

    setarrnValue(arrnSSBlock[0],1); //x: cols 
    setarrnValue(arrnSSBlock[1],1); //y: rows
    setarrnValue(arrnSSBlock[2],1); //z: slices

    //====================================== while loop  begin =================================================//

    int *arrnCpuEdge        = FRAME(arrnCpuEdge);
    int *arrnGpuEdge        = FRAME(arrnGpuEdge);
    int *arrnCpuEdgeLeft    = FRAME(arrnCpuEdgeLeft); 
    int *arrnGpuEdgeLeft    = FRAME(arrnGpuEdgeLeft); 
    int *arrnCpuPos         = FRAME(arrnCpuPos);
    int *arrnGpuPos         = FRAME(arrnGpuPos);
    int *arrnCpuEdgeVar     = FRAME(arrnCpuEdgeVar);
    int *arrnGpuEdgeVar     = FRAME(arrnGpuEdgeVar);
    int cpuCnt              = FRAME(cpuCnt);
    int gpuCnt              = FRAME(gpuCnt);
    
    while(1){
       
        //========================= set gpu begin sign ==========================//
        __sync_bool_compare_and_swap(&FRAME(GpuFinish),true,false);

#ifdef TIMERECORD 
        int cnt = RT(FRAME(gpuCnt));
        FRAME(sgRec[cnt]).tStart = getTime();
#endif


        arrnGpuPos = FRAME(arrnGpuPos);
        arrnGpuEdge= FRAME(arrnGpuEdge);
        int arrnGpuEdge2 [DIM];

        calcEdge(arrnGpuEdge2,arrnGpuPos,arrnGpuEdge, arrnEdge,1,0,0,2);
        
        int gpuPosx = arrnGpuPos[0]; 
        int gpuPosy = arrnGpuPos[1]; 
        int gpuPosz = arrnGpuPos[2]; 
	    


#ifdef CUDA_DARTS_DEBUG
        std::cout<<"gpu: gpuEdge["<<0<<"] ="<<arrnGpuEdge2[0] <<std::endl;
        std::cout<<"gpu: gpuEdge["<<1<<"] ="<<arrnGpuEdge2[1] <<std::endl;
        std::cout<<"gpu: gpuEdge["<<2<<"] ="<<arrnGpuEdge2[2] <<std::endl;
#endif
    
        //========================= calculate gpuChunk  ==========================//
        int arrnGpuEdgeChunk2[DIM];

        chooseSmaller(arrnGpuEdgeChunk2,arrnGpuEdgeChunk,arrnGpuEdge2,2,0,0,2);
        
	    int gpuChunkCols	= arrnGpuEdgeChunk[0];
	    int gpuChunkRows    = arrnGpuEdgeChunk[1];   
	    int gpuChunkSlices	= arrnGpuEdgeChunk[2];	

	    int gpuChunkCols2	= arrnGpuEdgeChunk2[0];
	    int gpuChunkRows2   = arrnGpuEdgeChunk2[1];   
	    int gpuChunkSlices2	= arrnGpuEdgeChunk2[2];	

	    size_t s_xpitch = gpuChunkCols2;
        size_t s_ypitch = gpuChunkRows2;
        size_t s_zpitch = gpuChunkSlices2;
        
        int vnStream;
        
        int arrnNumBk[DIM];

        calcarrnDivCeil(arrnNumBk, arrnGpuEdge2,arrnGpuEdgeChunk);
        int arrnOne[DIM];
        setarrnValue(arrnOne,1);
        cmpsAndassign(arrnNumBk,arrnGpuEdge2, arrnGpuEdgeChunk,arrnOne, arrnNumBk, 0,2,0,0);
	    int numbk_x = arrnNumBk[0]; 
	    int numbk_y = arrnNumBk[1]; 
	    int numbk_z = arrnNumBk[2]; 

        vnStream = numbk_x * numbk_y * numbk_z; 
        int tnStream = vnStream * nStream;

        struct_streams_par *streams_par = new struct_streams_par [tnStream];

#ifdef CUDA_DARTS_DEBUG
        std::cout<<"gpu: nStream: "<<nStream<<std::endl;
        std::cout<<"gpu: vnStream: "<<vnStream<<std::endl;
        std::cout<<"gpu: tnStream: "<<tnStream<<std::endl;
        std::cout<<"gpu: gpuEdgeChunk2["<<0<<"] ="<<arrnGpuEdgeChunk2[0] <<std::endl;
        std::cout<<"gpu: gpuEdgeChunk2["<<1<<"] ="<<arrnGpuEdgeChunk2[1] <<std::endl;
        std::cout<<"gpu: gpuEdgeChunk2["<<2<<"] ="<<arrnGpuEdgeChunk2[2] <<std::endl;
        std::cout<<"gpu: numBlock["<<0<<"] ="<<arrnNumBk[0] <<std::endl;
        std::cout<<"gpu: numBlock["<<1<<"] ="<<arrnNumBk[1] <<std::endl;
        std::cout<<"gpu: numBlock["<<2<<"] ="<<arrnNumBk[2] <<std::endl;
#endif

        //========================= allocate streams,events ==========================//
        cudaStream_t *stream ;
	    stream = new cudaStream_t[tnStream];
        cudaError *err = new cudaError[tnStream];
        cudaEvent_t *cuEvent = new cudaEvent_t[tnStream];
        for(int i=0;i<tnStream;++i){
            err[i]=cudaStreamCreate(&stream[i]);
#ifdef CUDA_ERROR_CHECKING
		    gpuErrchk(err[i]);
#endif
            err[i] = cudaEventCreate(&cuEvent[i]);
#ifdef CUDA_ERROR_CHECKING
		    gpuErrchk(err[i]);
#endif
        }

        //========================= gpu concurrent streams begin ==========================//
     
        for(int i=0;i<vnStream;++i){

#ifdef CUDA_ERROR_CHECKING
            err1 = cudaGetLastError();
		    gpuErrchk(err1);
#endif
            int pos[DIM];

            pos[0] = i%numbk_x;
		    pos[2] = i/(numbk_x*numbk_y);
            pos[1] = i/numbk_x-pos[2]*numbk_y;
            
            int posx = pos[0]; 
		    int posy = pos[1]; 
            int posz = pos[2]; 
            int ps = i%gnStream;
            int arrnChunk[DIM];
            calcEdgeChunk(arrnChunk,pos,arrnGpuEdgeChunk, arrnGpuEdge2, arrnGpuEdgeChunk2,1,0,0,0);
		    int nColsChunk		= arrnChunk[0]; 
            int nRowsChunk		= arrnChunk[1]; 
		    int nSlicesChunk	= arrnChunk[2]; 
                
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"stream #: "<<i<<" posx: "<<posx<<std::endl;
            std::cout<<"stream #: "<<i<<" posy: "<<posy<<std::endl;
            std::cout<<"stream #: "<<i<<" posz: "<<posz<<std::endl;
            std::cout<<"stream #: "<<i<<" nColsChunk   = "<< nColsChunk	 <<std::endl;
            std::cout<<"stream #: "<<i<<" nRowsChunk   = "<< nRowsChunk	 <<std::endl;
            std::cout<<"stream #: "<<i<<" nSlicesChunk = "<< nSlicesChunk   <<std::endl;
#endif
        
            
            //chooseSmaller(arrnGpuBlock,arrnChunk,arrnGpuTile,-2,0,0,-2);
            setarrn1Fromarrn2(arrnGpuBlock,arrnGpuTile); 
            calcarrnDivCeil(arrnGpuGrid,arrnChunk,arrnGpuBlock); 
            arrnGpuBlock[2]=1;
            
            setarrn1Fromarrn2(arrnSSGrid[0],arrnGpuGrid);
            setarrn1Fromarrn2(arrnSSGrid[1],arrnGpuGrid);
            setarrn1Fromarrn2(arrnSSGrid[2],arrnGpuGrid);
        
            //cols, y 
            arrnSSBlock[0][1] = (nRowsChunk>numThreads)?numThreads:nRowsChunk;
            //rows, x
            arrnSSBlock[1][0] = (nColsChunk>numThreads)?numThreads:nColsChunk;
            //slices,x
            arrnSSBlock[2][0] = (nColsChunk>numThreads)?numThreads:nColsChunk; 
            
            arrnSSGrid[0][1] = std::ceil(1.0*nRowsChunk/arrnSSBlock[0][1]);
            arrnSSGrid[1][0] = std::ceil(1.0*nColsChunk/arrnSSBlock[1][0]);
            arrnSSGrid[2][0] = std::ceil(1.0*nColsChunk/arrnSSBlock[2][0]);
                

#ifdef CUDA_DARTS_DEBUG

            std::cout<<"stream #: "<<i<<" arrnSSBlock[0][0]:"<<arrnSSBlock[0][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[0][1]:"<<arrnSSBlock[0][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[0][2]:"<<arrnSSBlock[0][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnSSBlock[1][0]:"<<arrnSSBlock[1][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[1][1]:"<<arrnSSBlock[1][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[1][2]:"<<arrnSSBlock[1][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnSSBlock[2][0]:"<<arrnSSBlock[2][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[2][1]:"<<arrnSSBlock[2][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSBlock[2][2]:"<<arrnSSBlock[2][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnGpuBlock[0]  :"<<arrnGpuBlock[0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnGpuBlock[1]  :"<<arrnGpuBlock[1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnGpuBlock[2]  :"<<arrnGpuBlock[2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnSSGrid[0][0]:"<<arrnSSGrid[0][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[0][1]:"<<arrnSSGrid[0][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[0][2]:"<<arrnSSGrid[0][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnSSGrid[1][0]:"<<arrnSSGrid[1][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[1][1]:"<<arrnSSGrid[1][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[1][2]:"<<arrnSSGrid[1][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnSSGrid[2][0]:"<<arrnSSGrid[2][0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[2][1]:"<<arrnSSGrid[2][1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnSSGrid[2][2]:"<<arrnSSGrid[2][2]<<std::endl;

            std::cout<<"stream #: "<<i<<" arrnGpuGrid[0]  :"<<arrnGpuGrid[0]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnGpuGrid[1]  :"<<arrnGpuGrid[1]<<std::endl;
            std::cout<<"stream #: "<<i<<" arrnGpuGrid[2]  :"<<arrnGpuGrid[2]<<std::endl;

#endif
            int idxSM = i*nStream;

            int nSMM1=nStream-1;
		    int idxSMM1 = idxSM+nSMM1;
		    if((i/gnStream) !=0){
			    int pre = (i-gnStream)*nStream;
			    err[pre]=cudaStreamWaitEvent(stream[idxSM],cuEvent[pre+nSMM1],0);
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"i: "<<i<<" ,event: "<<idxSM<<" waiting for event: "<<pre+nSMM1<<std::endl;
#endif
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err[pre]);
#endif
            }
            int dev_posz = ps*dev_zpitch;

            phd.srcPos = make_cudaPos((gpuPosx+posx*gpuChunkCols)*szdb,gpuPosy+posy*gpuChunkRows,gpuPosz+posz*gpuChunkSlices);
            phd.dstPos = make_cudaPos(0,0,dev_posz);
            
            phd.extent.width  = nColsChunk*szdb;
            phd.extent.height = nRowsChunk;
            phd.extent.depth  = nSlicesChunk;

#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelWithStream multiple streams: "<<idxSM <<" memory copy begin!"<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" posx: "<<posx<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" posy: "<<posy<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" posz: "<<posz<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.x: "<<phd.srcPos.x<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.y: "<<phd.srcPos.y<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPos.z: "<<phd.srcPos.z<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.x: "<<phd.dstPos.x<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.y: "<<phd.dstPos.y<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPos.z: "<<phd.dstPos.z<<std::endl;

            std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.ptr: "  <<phd.srcPtr.ptr  <<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.pitch: "<<phd.srcPtr.pitch<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.xsize: "<<phd.srcPtr.xsize<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" srcPtr.ysize: "<<phd.srcPtr.ysize<<std::endl;

            std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.ptr: "  <<phd.dstPtr.ptr  <<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.pitch: "<<phd.dstPtr.pitch<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.xsize: "<<phd.dstPtr.xsize<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" dstPtr.ysize: "<<phd.dstPtr.ysize<<std::endl;

            std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.width: " <<phd.extent.width<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.height: "<<phd.extent.height<<std::endl;
            std::cout<<"copy3d HToD: stream #"<<i<<" p.extent.depth: " <<phd.extent.depth<<std::endl;
#endif

            pdh.srcPos = make_cudaPos(1*szdb,1,dev_posz+1);
            pdh.dstPos = make_cudaPos((gpuPosx+posx*gpuChunkCols+1)*szdb,gpuPosy+posy*gpuChunkRows+1,gpuPosz+posz*gpuChunkSlices+1);

            int arrnChunk2[DIM];
            calcEdgeChunk(arrnChunk2,pos,arrnGpuEdgeChunk, arrnGpuEdge2, arrnGpuEdgeChunk,1,0,-2,0);
		    int nColsChunk2		= arrnChunk2[0]; 
            int nRowsChunk2		= arrnChunk2[1]; 
		    int nSlicesChunk2	= arrnChunk2[2]; 

            pdh.extent.width  = nColsChunk2*szdb;
            pdh.extent.height = nRowsChunk2;
            pdh.extent.depth  = nSlicesChunk2;

            streams_par[idxSM].htodPtr = &phd;
            streams_par[idxSM].dtohPtr = &pdh; 
            streams_par[idxSM].devDstPos    = ps*dev_size;
            streams_par[idxSM].devSRowsPos  = ps*srows_size;
            streams_par[idxSM].devSColsPos  = ps*scols_size;
            streams_par[idxSM].devSSlicesPos= ps*sslices_size; 	 
   
	        streams_par[idxSM].dimBlock_cols   = dim3 (arrnSSBlock[0][0],arrnSSBlock[0][1],arrnSSBlock[0][2]);
	        streams_par[idxSM].dimBlock_rows   = dim3 (arrnSSBlock[1][0],arrnSSBlock[1][1],arrnSSBlock[1][2]);
            streams_par[idxSM].dimBlock_slices = dim3 	(arrnSSBlock[2][0],arrnSSBlock[2][1],arrnSSBlock[2][2]);
            streams_par[idxSM].dimBlock        = dim3 (arrnGpuBlock[0],arrnGpuBlock[1],arrnGpuBlock[2]);
	        streams_par[idxSM].dimGrid_cols    = dim3 (arrnSSGrid[0][0], arrnSSGrid[0][1],arrnSSGrid[0][2]);
            streams_par[idxSM].dimGrid_rows    = dim3 (arrnSSGrid[1][0], arrnSSGrid[1][1],arrnSSGrid[1][2]);
            streams_par[idxSM].dimGrid_slices  = dim3 
 (arrnSSGrid[2][0], arrnSSGrid[2][1],arrnSSGrid[2][2]);
           streams_par[idxSM].dimGrid         = dim3 (arrnGpuGrid[0], arrnGpuGrid[1],arrnGpuGrid[2]);

            streams_par[idxSM].nRowsChunk  = nRowsChunk  ;	
            streams_par[idxSM].nColsChunk  = nColsChunk  ;
            streams_par[idxSM].nSlicesChunk= nSlicesChunk; 

#ifdef CUDA_DARTS_DEBUG
            std::cout<<"stream: "<<i<<" device dst   address from: "<<phd.dstPtr.ptr  + streams_par[idxSM].devDstPos    <<" to: "<<phd.dstPtr.ptr  + streams_par[idxSM].devDstPos    + dev_size      <<std::endl;
            std::cout<<"stream: "<<i<<" sharedSlices address from: "<<dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos  <<" to: "<<dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos  + sslices_size    <<std::endl;
            std::cout<<"stream: "<<i<<" sharedRows   address from: "<<dptr_dev_sharedRows  + streams_par[idxSM].devSRowsPos  <<" to: "<<dptr_dev_sharedRows  + streams_par[idxSM].devSRowsPos  + srows_size    <<std::endl;
            std::cout<<"stream: "<<i<<" sharedCols   address from: "<<dptr_dev_sharedCols  + streams_par[idxSM].devSColsPos<<" to: "<<dptr_dev_sharedCols  + streams_par[idxSM].devSColsPos+ scols_size <<std::endl;
#endif
        
        
#ifdef CUDA_DARTS_DEBUG
            err1 = cudaDeviceSynchronize();
			gpuErrchk(err1);
#endif        
            err0=cudaMemcpy3DAsync(streams_par[idxSM].htodPtr,stream[idxSM]);
        
#ifdef CUDA_ERROR_CHECKING
            gpuErrchk(err0);
#endif
        
#ifdef CUDA_DARTS_DEBUG
            err1 = cudaDeviceSynchronize();
			gpuErrchk(err1);
#endif

#ifdef CUDA_DARTS_DEBUG
             {
                    
#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
			    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
                dim3 dimGrid_t (1,1,1);
	            dim3 dimBlock_t(10,1,1);
                    
                double *sptr    = phd.srcPtr.ptr;
                int spitch      = phd.srcPtr.pitch/szdb;
                int sposx       = phd.srcPos.x/szdb; 
                int sposy       = phd.srcPos.y; 
                int sposz       = phd.srcPos.z;
                int sxsize      = phd.srcPtr.xsize;
                int sysize      = phd.srcPtr.ysize;
                
                double *dptr    = phd.dstPtr.ptr;
                int dpitch      = phd.dstPtr.pitch/szdb;
                int dposx       = phd.dstPos.x; 
                int dposy       = phd.dstPos.y; 
                int dposz       = phd.dstPos.z;

                int dxsize      = phd.dstPtr.xsize;
                int dysize      = phd.dstPtr.ysize;

                int stride= 0;
                std::cout<<"Async cp phd: i: "<<i<<std::endl;
                std::cout<<"Async cp phd: phd.srcPtr value: "<<std::endl;
                test_copy3d(dimGrid_t,dimBlock_t,sptr,sposz*spitch*sysize +sposx+1*spitch*sysize+sposy*spitch+2*spitch+stride,1,1);

#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
                err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
                std::cout<<"Async cp p: p.dstPtr value: "<<std::endl;
                test_copy3d(dimGrid_t,dimBlock_t,dptr,dposz*dysize*dpitch+dposx + 1*dpitch*dysize+dposy*dpitch+2*dpitch+stride,1,1);

#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
                err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
            }
#endif
        
#ifdef CUDA_DARTS_DEBUG
            std::cout<<"GpuKernelWithStream multiple streams: "<<idxSM <<" memory copy end!"<<std::endl;
#endif
            
            err[idxSM]=cudaEventRecord(cuEvent[idxSM],stream[idxSM]);
            for(int k=1;k<nSMM1;++k){
                err[idxSM+k]=cudaStreamWaitEvent(stream[idxSM+k],cuEvent[idxSM],0);
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err[idxSM+k]);
#endif
            }
            

#ifdef CUDA_DARTS_DEBUG
            err1 = cudaGetLastError();
			gpuErrchk(err1);
#endif
        
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedCols   address from: "<<dptr_dev_sharedCols  + streams_par[idxSM].devSColsPos<<" to: "<<dptr_dev_sharedCols  + streams_par[idxSM].devSColsPos+ scols_size <<std::endl;
#endif
                gpu_kernel37_cp_cols_stream_p(stream[idxSM+1],streams_par[idxSM].dimGrid_cols,streams_par[idxSM].dimBlock_cols,dptr_dev+streams_par[idxSM].devDstPos,dptr_dev_sharedRows+streams_par[idxSM].devSRowsPos, dptr_dev_sharedCols+ streams_par[idxSM].devSColsPos, dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos,dev_xpitch,dev_ypitch, dev_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+1]=cudaEventRecord(cuEvent[idxSM+1],stream[idxSM+1]);

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
                
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedRows   address from: "<<dptr_dev_sharedRows  + streams_par[idxSM].devSRowsPos  <<" to: "<<dptr_dev_sharedRows  + streams_par[idxSM].devSRowsPos  + srows_size    <<std::endl;
#endif
                
                gpu_kernel37_cp_rows_stream_p(stream[idxSM+2],streams_par[idxSM].dimGrid_rows,streams_par[idxSM].dimBlock_rows,dptr_dev+streams_par[idxSM].devDstPos,dptr_dev_sharedRows+streams_par[idxSM].devSRowsPos, dptr_dev_sharedCols+ streams_par[idxSM].devSColsPos, dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos,dev_xpitch,dev_ypitch, dev_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+2]=cudaEventRecord(cuEvent[idxSM+2],stream[idxSM+2]);
                
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
                
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"stream: "<<i<<" sharedSlices address from: "<<dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos  <<" to: "<<dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos  + sslices_size    <<std::endl;
#endif
                gpu_kernel37_cp_slices_stream_p(stream[idxSM+3],streams_par[idxSM].dimGrid_slices,streams_par[idxSM].dimBlock_slices,dptr_dev+streams_par[idxSM].devDstPos,dptr_dev_sharedRows+streams_par[idxSM].devSRowsPos, dptr_dev_sharedCols+ streams_par[idxSM].devSColsPos, dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos,dev_xpitch,dev_ypitch, dev_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
                err[idxSM+3]=cudaEventRecord(cuEvent[idxSM+3],stream[idxSM+3]);
                
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
        
#ifdef CUDA_ERROR_CHECKING
                for(int k=1;k<nSMM1;++k){
					gpuErrchk(err[idxSM+k]);
                }
#endif
                for(int k=1;k<nSMM1;++k){
                    err[idxSM+k]=cudaStreamWaitEvent(stream[idxSMM1],cuEvent[idxSM+k],0);
#ifdef CUDA_ERROR_CHECKING
					gpuErrchk(err[idxSM+k]);
#endif
                }
      
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
  
                gpu_kernel37_stream_p(stream[idxSMM1],streams_par[idxSM].dimGrid,streams_par[idxSM].dimBlock,dptr_dev+streams_par[idxSM].devDstPos,dptr_dev_sharedRows+streams_par[idxSM].devSRowsPos, dptr_dev_sharedCols+ streams_par[idxSM].devSColsPos, dptr_dev_sharedSlices+ streams_par[idxSM].devSSlicesPos,dev_xpitch,dev_ypitch, dev_zpitch,s_xpitch,s_ypitch,s_zpitch,streams_par[idxSM].nRowsChunk, streams_par[idxSM].nColsChunk,streams_par[idxSM].nSlicesChunk,gpuTile_x,gpuTile_y,gpuTile_z);
        
#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
        
#ifdef CUDA_DARTS_DEBUG
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.x: "        <<pdh.srcPos.x<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.y: "        <<pdh.srcPos.y<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPos.z: "        <<pdh.srcPos.z<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.x: "        <<pdh.dstPos.x<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.y: "        <<pdh.dstPos.y<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPos.z: "        <<pdh.dstPos.z<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.ptr: "      <<pdh.srcPtr.ptr  <<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.pitch: "    <<pdh.srcPtr.pitch<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.xsize: "    <<pdh.srcPtr.xsize<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" srcPtr.ysize: "    <<pdh.srcPtr.ysize<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.ptr: "      <<pdh.dstPtr.ptr  <<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.pitch: "    <<pdh.dstPtr.pitch<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.xsize: "    <<pdh.dstPtr.xsize<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dstPtr.ysize: "    <<pdh.dstPtr.ysize<<std::endl;

                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.width: "    <<pdh.extent.width<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.height: "   <<pdh.extent.height<<std::endl;
                std::cout<<"copy3d DToH: stream #"<<i<<" dph.extent.depth: "    <<pdh.extent.depth<<std::endl;

#endif

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
                err0=cudaMemcpy3DAsync(streams_par[idxSM].dtohPtr,stream[idxSMM1]);
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err0);
#endif

#ifdef CUDA_DARTS_DEBUG
                err1 = cudaDeviceSynchronize();
				gpuErrchk(err1);
#endif
				err[idxSMM1]=cudaEventRecord(cuEvent[idxSMM1],stream[idxSMM1]);
#ifdef CUDA_ERROR_CHECKING
                gpuErrchk(err[idxSMM1]);
#endif

#ifdef CUDA_DARTS_DEBUG
                {

#ifdef CUDA_DARTS_DEBUG
                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
#endif
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                    dim3 dimGrid_t (1,1,1);
	                dim3 dimBlock_t(10,1,1);

                    double *sptr    = pdh.srcPtr.ptr;
                    int spitch      = pdh.srcPtr.pitch/szdb;
                    int sposx       = pdh.srcPos.x/szdb; 
                    int sposy       = pdh.srcPos.y; 
                    int sposz       = pdh.srcPos.z; 
                    int sxsize      = pdh.srcPtr.xsize;
                    int sysize      = pdh.srcPtr.ysize;
                    
                    double *dptr    = pdh.dstPtr.ptr;
                    int dpitch      = pdh.dstPtr.pitch/szdb;
                    int dposx       = pdh.dstPos.x/szdb; 
                    int dposy       = pdh.dstPos.y; 
                    int dposz       = pdh.dstPos.z;
                    int dxsize      = pdh.dstPtr.xsize;
                    int dysize      = pdh.dstPtr.ysize;

                    
                    int stride= 0;
                    std::cout<<"Async cp pdh: i: "<<i<<std::endl;
                    std::cout<<"Async cp pdh: pdh.srcPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,sptr,sposz*spitch*sysize +sposx+0*spitch*sysize+sposy*spitch+1*spitch+stride,1,1);

                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				gpuErrchk(err1);
#endif
                    std::cout<<"Async cp pdh: pdh.dstPtr value: "<<std::endl;
                    test_copy3d(dimGrid_t,dimBlock_t,dptr,dposz*dpitch*dysize +dposx+0*dpitch*dysize+dposy*dpitch+1*dpitch+stride,1,1);

                    err1 = cudaGetLastError();
				    gpuErrchk(err1);
                    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
				    gpuErrchk(err1);
#endif
                }
#endif


#ifdef CUDA_ERROR_CHECKING
                err1 = cudaGetLastError();
				gpuErrchk(err1);
#endif
        
        }

		err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
		gpuErrchk(err1);
#endif	
        //========================= gpu concurrent streams finish =====================//
        
        //========================= delete streams,events ==========================//
        delete [] cuEvent;
        delete [] err;
        delete [] streams_par;
	    for(int i=0;i<tnStream;++i){
	    	cudaStreamDestroy(stream[i]);
	    }
        delete [] stream;

        //========================= set gpu finish sign ==========================//
        __sync_bool_compare_and_swap(&FRAME(GpuFinish),false,true);
        
        
        //========================= record gpu execution time ==========================//
#ifdef TIMERECORD 
        int *chunk = FRAME(sgRec[cnt]).chunk;
        setarrn1Fromarrn2(chunk,arrnGpuEdge);
        FRAME(sgRec[cnt]).tEnd=getTime();
        FRAME(sgRec[cnt]).tExe=FRAME(sgRec[cnt]).tEnd - FRAME(sgRec[cnt]).tStart;
#endif
        FRAME(gpuCnt++);

        __sync_synchronize();
        //========================= calc next gpuEdge and EdgeLeft==========================//
        arrnCpuEdge        = FRAME(arrnCpuEdge);
        arrnGpuEdge        = FRAME(arrnGpuEdge);
        arrnCpuEdgeLeft    = FRAME(arrnCpuEdgeLeft); 
        arrnGpuEdgeLeft    = FRAME(arrnGpuEdgeLeft); 
        arrnCpuPos         = FRAME(arrnCpuPos);
        arrnGpuPos         = FRAME(arrnGpuPos);
        arrnCpuEdgeVar     = FRAME(arrnCpuEdgeVar);
        arrnGpuEdgeVar     = FRAME(arrnGpuEdgeVar);
        cpuCnt              = FRAME(cpuCnt);
        gpuCnt              = FRAME(gpuCnt);

#ifdef CUDA_DARTS_DEBUG2
	    pthread_mutex_lock(&mutex2);

        std::cout<<"gpu: current cpuCnt = "<<FRAME(cpuCnt)<<std::endl;
        std::cout<<"gpu: current gpuCnt = "<<FRAME(gpuCnt)<<std::endl;
        
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current cpuPos["<<i<<"] = "<<arrnCpuPos[i]<<std::endl;
        }
	    
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current cpuEdge["<<i<<"] = "<<arrnCpuEdge[i]<<std::endl;
        }
        
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current cpuEdgeLeft["<<i<<"] = "<<arrnCpuEdgeLeft[i]<<std::endl;
        }
       
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current gpuPos["<<i<<"] = "<<arrnGpuPos[i]<<std::endl;
        }
	    
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current gpuEdge["<<i<<"] = "<<arrnGpuEdge[i]<<std::endl;
        }
        
        for (int i=0;i<DIM;++i){
            std::cout<<"gpu: current gpuEdgeLeft["<<i<<"] = "<<arrnGpuEdgeLeft[i]<<std::endl;
        }

        pthread_mutex_unlock(&mutex2);
#endif

        if((checkarrnEqValue(arrnCpuEdgeLeft, 0) && checkarrnEqValue(arrnGpuEdgeLeft,0))){
            break;        
        }else{

            pthread_mutex_lock(&mutex);
            __sync_synchronize();
            calcNextEP(arrnCpuEdge, arrnGpuEdge,arrnCpuEdgeLeft, arrnGpuEdgeLeft,arrnCpuPos, arrnGpuPos,arrnEdge, arrnCpuEdgeVar, arrnGpuEdgeVar,arrnCpuEdgeMin,arrnGpuEdgeMin,cpuCnt,gpuCnt, nId,"gpu");
            __sync_synchronize();
            pthread_mutex_unlock(&mutex);
   
#ifdef CUDA_DARTS_DEBUG2
	        pthread_mutex_lock(&mutex2);

            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next cpuPos["<<i<<"] = "<<arrnCpuPos[i]<<std::endl;
            }
	        
            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next cpuEdge["<<i<<"] = "<<arrnCpuEdge[i]<<std::endl;
            }
            
            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next cpuEdgeLeft["<<i<<"] = "<<arrnCpuEdgeLeft[i]<<std::endl;
            }
       
            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next gpuPos["<<i<<"] = "<<arrnGpuPos[i]<<std::endl;
            }
	        
            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next gpuEdge["<<i<<"] = "<<arrnGpuEdge[i]<<std::endl;
            }
            
            for (int i=0;i<DIM;++i){
                std::cout<<"gpu: next gpuEdgeLeft["<<i<<"] = "<<arrnGpuEdgeLeft[i]<<std::endl;
            }

            pthread_mutex_unlock(&mutex2);
#endif

        }

    }
    //========================================= while loop  end ============================================//


	SYNC(Swap37);

	EXIT_TP();

}


void
SyncCD::fire(void)
{
#ifdef CUDA_DARTS_DEBUG
	std::cout<<"invoke Sync!"<<std::endl;
#endif
	LOAD_FRAME(StencilTP);
	SIGNAL(signalUp);
    EXIT_TP();
}


