#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H
//#define CUDA_APT_PER_THREAD_DEFAULT_STREAM
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "conf.h"
#include "stencil.h"
//#include <math.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define N_CORES TOTAL_NUM_CU
using namespace darts;

#define GPUMETA 0x4

DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_CODELET_ITER(Stencil3D7ptCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(Stencil3D7ptCpuSyncCD,0,SHORTWAIT);
DEF_CODELET(Stencil3D7ptSwapCD,0,SHORTWAIT);


DEF_CODELET(Stencil3D7ptGpuKernelWithAllTimeStepsBigCD,2,LONGWAIT);
DEF_CODELET(Stencil3D7ptGpuKernelWithAllTimeStepsCD,2,LONGWAIT);
DEF_CODELET(Stencil3D7ptGpuKernelPureGpuWithStreamsCD,2,LONGWAIT);
DEF_CODELET(Stencil3D7ptGpuKernelHybridWithStreamsCD,2,LONGWAIT);

DEF_TP(StencilTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
    const uint64_t nSlices;
    double *New;
	int32_t ts;
    int32_t tsInit;	

	bool hard;
    bool IsStatic;
    double GpuRatio;
    bool streamming;
    int32_t nTiles;
	double softGpuRatio=0.0;
	int32_t nCPU = 0;
	int32_t nGPU = 0;
	
	Stencil3D7ptCpuLoopCD *CpuLoop37 = NULL;
    Stencil3D7ptCpuSyncCD	CpuSync37;
	Stencil3D7ptSwapCD Swap37;

	Stencil3D7ptGpuKernelWithAllTimeStepsCD GpuKernelWithAllTimeSteps37;
	Stencil3D7ptGpuKernelWithAllTimeStepsBigCD GpuKernelWithAllTimeStepsBig37;
    Stencil3D7ptGpuKernelPureGpuWithStreamsCD GpuKernelPureGpuWithStreams37;
    Stencil3D7ptGpuKernelHybridWithStreamsCD GpuKernelHybridWithStreams37;
    SyncCD	sync;
    Codelet *signalUp;

    int64_t nId;
    
	int gpuSlicesMin=3;
    int cpuSlicesMin=N_CORES;


    
    double  d_size = 0;
    int64_t d_size_sharedCols = 0;
    int64_t d_size_sharedRows = 0;
    int64_t d_size_sharedSlices = 0;
    double req_size=0;

	double *d_dst = NULL;
	dim3 dimGrid_hack1;
	double cmCpu = 10;	//cpu(total 32 cores) compute ability
	double cmGpu = 10;	//Gpu compute ability
	bool CpuIvGpu = false; // CPU to invoke GPU;
    bool CpuFinish = false; // GPU check CPU computation finish; 
    bool GpuFinish = false; // CPU check GPU computation finish; 
    uint64_t gpuPos=0;
	uint64_t cpuPos=0;

    

    int arrnCpuPos[DIM];
    int arrnGpuPos[DIM];
    
    int arrnEdge[DIM];
    
    int arrnCpuEdgeInit[DIM];
    int arrnGpuEdgeInit[DIM];
    
    int arrnCpuEdge[DIM];
    int arrnGpuEdge[DIM];

    //int arrnCpuEdgeChunk[DIM];
    uint32_t arrnGpuEdgeChunk[DIM];
    int arrnGpuEdgeChunkAlloc[DIM];
   

	int arrnCpuGridTileBase[DIM];
	int arrnGpuGridTileBase[DIM];
	
    int arrnCpuTile[DIM];
	int arrnGpuTile[DIM];
    
    int arrnCpuBlock[DIM];
    int arrnGpuBlock[DIM];

    int arrnCpuGrid[DIM];
    int arrnGpuGrid[DIM];

    //int arrnEdgeLeft[DIM];
    int arrnCpuEdgeLeft[DIM];
    int arrnGpuEdgeLeft[DIM];
   
    int arrnCpuEdgeVar[DIM];
    int arrnGpuEdgeVar[DIM];

    int arrnCpuEdgeMin[DIM];
    int arrnGpuEdgeMin[DIM];
	
    int nCPUInit = 0;
	int nGPUInit = 0;
    int gpuCnt = 0;
    int cpuCnt = 0;
    
	int32_t lastCnt = 0;
    int64_t gpuMemMax = 0;
    int64_t gpuAvMem = GPU_AVMEM_SIZE*GB; 

    int nStream = 0 ;
    int gnStream = 0 ; //group number of stream
    int tnStream = 0 ; //total number of stream 
    int vnStream = 0 ;
    const int szdb = sizeof(double);
    int numPar = NUMPAR;
    bool bIsHybrid = false;
    cudaExtent extent ;
    cudaPitchedPtr devPitchedPtr;
  
	double *dptr_dev ;
	double *dptr_dev_sharedCols ;
	double *dptr_dev_sharedRows ;
	double *dptr_dev_sharedSlices;

	size_t dev_xpitch = 0; 
    size_t dev_ypitch = 0; 
    size_t dev_zpitch = 0; 
    int64_t dev_size  = 0;


    int64_t srows_size   = 0; 
	int64_t scols_size   = 0; 
    int64_t sslices_size = 0; 

    int64_t dev_sharedCols_size     = 0;  	
    int64_t dev_sharedRows_size     = 0;  	
    int64_t dev_sharedSlices_size   = 0;	

#ifdef TIMERECORD
    SRECORD scRec[NUMREC];
    SRECORD sgRec[NUMREC];
#endif

    StencilTP(double * inimatrix,const uint64_t n_rows,const uint64_t n_cols, const uint64_t n_slices,double * newmatrix,uint64_t ts,bool hard,bool IsStatic,double GpuRatio,bool streamming,int32_t nTiles, Codelet *up)
	:Initial(inimatrix)
	,nRows(n_rows)
	,nCols(n_cols)
    ,nSlices(n_slices)
	,New(newmatrix)
	,ts(ts)
	,tsInit(ts)
    ,hard(hard)
    ,IsStatic(IsStatic)
    ,GpuRatio(GpuRatio)
    ,streamming(streamming)
    ,nTiles(nTiles)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{

#ifdef CUDA_DARTS_DEBUG
		std::cout<<"invoke TP!"<<std::endl;	
		std::cout<<std::setprecision(18)<<std::endl;
#endif

        arrnEdge[0] = nCols;
        arrnEdge[1] = nRows;
        arrnEdge[2] = nSlices;

	    arrnCpuGridTileBase[0] = GRID_TILECPU_X;
	    arrnCpuGridTileBase[1] = GRID_TILECPU_Y;
	    arrnCpuGridTileBase[2] = GRID_TILECPU_Z;
	    
        arrnGpuGridTileBase[0] = GRID_TILE37_X;
        arrnGpuGridTileBase[1] = GRID_TILE37_Y;
        arrnGpuGridTileBase[2] = GRID_TILE37_Z;
		
        arrnCpuEdge[0]=arrnEdge[0];
        arrnCpuEdge[1]=arrnEdge[1];
        arrnCpuEdge[2]=arrnEdge[2];
       
        arrnGpuEdge[0]=arrnEdge[0];
        arrnGpuEdge[1]=arrnEdge[1];
        arrnGpuEdge[2]=arrnEdge[2];
        
        int arrnTWO[DIM];
        setarrnValue(arrnTWO,2);
        setarrn1Fromarrn2(arrnCpuEdgeMin,arrnGpuGridTileBase);
        calcarrn1Fromarrn2Marrn3(arrnGpuEdgeMin,arrnGpuGridTileBase,arrnTWO);


        if(GpuRatio == 0.0){

                bIsHybrid = false;
                //nCPU = N_CORES;
	    	    cpuPos = 0;
              
                setarrnValue(arrnCpuPos,0); 

                chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
                chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);

                calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 

                int cpuGridDimx = arrnCpuGrid[0]; 
				int cpuGridDimy = arrnCpuGrid[1]; 
				int cpuGridDimz = arrnCpuGrid[2]; 
                
				nCPU = cpuGridDimx*cpuGridDimy*cpuGridDimz;
				CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                for(size_t i=0;i<nCPU; ++i){
	    	   	    CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    	   	    add(CpuLoop37 + i);
	    	    }

                Swap37 = Stencil3D7ptSwapCD{nCPU,nCPU,this,SHORTWAIT};
#ifdef CUDA_DARTS_DEBUG
		        std::cout<<"cpuBlockDimx = "<<arrnCpuBlock[0]<<std::endl;
		        std::cout<<"cpuBlockDimy = "<<arrnCpuBlock[1]<<std::endl;
		        std::cout<<"cpuBlockDimz = "<<arrnCpuBlock[2]<<std::endl;
		        std::cout<<"cpuGridDimx = "<<arrnCpuGrid[0]<<std::endl;
		        std::cout<<"cpuGridDimy = "<<arrnCpuGrid[1]<<std::endl;
		        std::cout<<"cpuGridDimz = "<<arrnCpuGrid[2]<<std::endl;
#endif
		}else{
           
	    		int deviceCount;
	    		cudaGetDeviceCount(&deviceCount);
	    		cudaDeviceProp props;
	    		cudaGetDeviceProperties(&props,0);
                int ccKernels =props.concurrentKernels; 
#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu device count: "<<deviceCount<<std::endl;
	    		std::cout<<"shared memory per block: "<<props.sharedMemPerBlock/KB<<"KB"<<std::endl;
	    		std::cout<<"registers per Block: "<<props.regsPerBlock<<std::endl;
	    		std::cout<<"Threads per Block:"<<props.maxThreadsPerBlock<<std::endl;
	    		std::cout<<"shared memory per MP: "<<props.sharedMemPerMultiprocessor/KB<<"KB"<<std::endl;
	    		std::cout<<"Threads per MP:"<<props.maxThreadsPerMultiProcessor<<std::endl;
	    		std::cout<<"MB count:"<<props.multiProcessorCount<<std::endl;
	    		std::cout<<"Global memory: "<<props.totalGlobalMem/MB<<"MB"<<std::endl;
                std::cout<<"concurrent kernels ability: "<<ccKernels<<std::endl; 
#endif


	    		size_t gpu_mem_total_t = 0;
	    		size_t gpu_mem_avail_t = 0;
	    		size_t gpu_mem_valid_t = 0;
	    		cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	    		gpu_mem_valid_t = gpu_mem_avail_t - XMB;
                
                gpuMemMax =gpuAvMem> gpu_mem_valid_t?gpu_mem_avail_t: gpuAvMem;

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu avail memory: "<<gpu_mem_avail_t/KB<<"KB"<<std::endl;
	    		std::cout<<"gpu valid memory: "<<gpu_mem_valid_t/KB<<"KB"<<std::endl;
                std::cout<<"gpu memory avail Max: "<<gpuMemMax/KB<<"KB"<<std::endl;
#endif

	    	    cpuPos = 0;
	    	    gpuPos = 0;

                arrnCpuPos[0]=0;
                arrnCpuPos[1]=0;
                arrnCpuPos[2]=0;
                
                arrnGpuPos[0]=0;
                arrnGpuPos[1]=0;
                arrnGpuPos[2]=0;

                int64_t gpuCols     =arrnEdge[0];  
                int64_t gpuRows     =arrnEdge[1];  
                int64_t gpuSlices   =arrnEdge[2];

                chooseSmaller(arrnGpuTile,arrnGpuEdge,arrnGpuGridTileBase,0,0,0,0);
                chooseSmaller(arrnGpuBlock,arrnGpuEdge,arrnGpuTile,-2,0,0,-2);
                calcarrnDivCeil(arrnGpuGrid,arrnGpuEdge,arrnGpuBlock); 
                arrnGpuBlock[2]= 1;
				
                int gpuTile_x = arrnGpuTile[0]; 
				int gpuTile_y = arrnGpuTile[1]; 
				int gpuTile_z = arrnGpuTile[2]; 
                int gpuBlockDimx = arrnGpuBlock[0]; 
                int gpuBlockDimy = arrnGpuBlock[1];
                int gpuBlockDimz = arrnGpuBlock[2];

                int gpuGridDimx = arrnGpuGrid[0];
                int gpuGridDimy = arrnGpuGrid[1];
                int gpuGridDimz = arrnGpuGrid[2];

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpuGridDimx: "<<gpuGridDimx<<std::endl;;
	    		std::cout<<"gpuGridDimy: "<<gpuGridDimy<<std::endl;;
	    		std::cout<<"gpuGridDimz: "<<gpuGridDimz<<std::endl;;
#endif


                d_size = szdb*gpuRows*gpuCols*gpuSlices;  
                d_size_sharedCols   = szdb*gpuRows*gpuSlices*gpuGridDimx*2;
                d_size_sharedRows   = szdb*gpuCols*gpuSlices*gpuGridDimy*2;
                d_size_sharedSlices = szdb*gpuRows*gpuCols  *gpuGridDimz*2;
                req_size = d_size + d_size_sharedCols + d_size_sharedRows + d_size_sharedSlices;
                //gpuWLMax = std::floor(1.0*gpuMemMax/(sizeof(double)*(gpuRows*gpuCols+gpuRows*gpuGridDimx*2+gpuCols*gpuGridDimy*2+1.0*gpuRows*gpuCols*(1/gpuTile_z))));

#ifdef CUDA_DARTS_DEBUG		
	    		
	            std::cout<<"d_size:"<<d_size<<std::endl;
	            std::cout<<"d_size_sharedCols:"<<d_size_sharedCols<<std::endl;
	            std::cout<<"d_size_sharedRows:"<<d_size_sharedRows<<std::endl;
	            std::cout<<"d_size_sharedSlices:"<<d_size_sharedSlices<<std::endl;
                
                std::cout<<"req size: "<<req_size/KB<<"KB,gupMemMax: "<<gpuMemMax/KB<<"KB!"<<std::endl;
#endif
                if (GpuRatio == 1.0){
                    
                    bIsHybrid = false;
                    
                    gpuPos      = 0;
					gpuSlices	= nSlices;
					gpuRows		= nRows;
					gpuCols		= nCols;
                    
                   // if((streamming==false)&&(req_size<gpuMemMax)){
                    if(streamming==false){
                        if(req_size<gpuMemMax){

#ifdef CUDA_DARTS_DEBUG		
	    		            std::cout<<"gpu require size is less than gpuMemMax and not use streamming! "<<std::endl;
#endif
                            nGPU	= 1;
                            gnStream = 1;
						    nStream	= 4;
                            tnStream = nStream;
                            GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                            add(&GpuKernelWithAllTimeSteps37);	
                           }else{
                           //     exit(-1);

#ifdef CUDA_DARTS_DEBUG		
	    		            std::cout<<"gpu require size is larger than gpuMemMax and not use streamming! "<<std::endl;
#endif
                            
                            nGPU = std::ceil(req_size/gpuMemMax);
                            gnStream = 1;
						    nStream	= 4;
                            tnStream = nStream;
                            GpuKernelWithAllTimeStepsBig37 = Stencil3D7ptGpuKernelWithAllTimeStepsBigCD{0,1,this,GPUMETA};
                            add(&GpuKernelWithAllTimeStepsBig37);	
                           }
                    
                    }else{

#ifdef CUDA_DARTS_DEBUG		
	    		        std::cout<<"gpu require size is greater than gpuMemMax or forse to use streamming! "<<std::endl;
#endif
                        nGPU = std::ceil(req_size/gpuMemMax);
                        nStream = NUMASTREAMS;
                        if (nTiles==0){
                            arrnGpuGrid[0] = NUMTILES; 
                            arrnGpuGrid[1] = NUMTILES; 
                            arrnGpuGrid[2] = NUMTILES; 
                            gnStream = NUMPAR ;
                        }else{
                            arrnGpuGrid[0] = nTiles; 
                            arrnGpuGrid[1] = nTiles; 
                            arrnGpuGrid[2] = nTiles; 
                        }

                        arrnGpuEdgeChunk[0]=((arrnGpuGrid[0]*gpuTile_x)>gpuCols  )?gpuCols:  ( arrnGpuGrid[0]*gpuTile_x);
                        arrnGpuEdgeChunk[1]=((arrnGpuGrid[1]*gpuTile_y)>gpuRows  )?gpuRows:  ( arrnGpuGrid[1]*gpuTile_y);
                        arrnGpuEdgeChunk[2]=((arrnGpuGrid[2]*gpuTile_z)>gpuSlices)?gpuSlices:( arrnGpuGrid[2]*gpuTile_z);

                        if(nTiles!=0){
                            double  t_d_size = szdb*arrnGpuEdgeChunk[0]*arrnGpuEdgeChunk[1]*arrnGpuEdgeChunk[2] ; 
                            uint64_t  t_d_size_sharedCols   = szdb*arrnGpuEdgeChunk[1] * arrnGpuEdgeChunk[2] *nTiles*2;
                            uint64_t  t_d_size_sharedRows   = szdb*arrnGpuEdgeChunk[0] * arrnGpuEdgeChunk[2] *nTiles*2;
                            uint64_t  t_d_size_sharedSlices = szdb*arrnGpuEdgeChunk[0] * arrnGpuEdgeChunk[1] *nTiles*2; 
                            uint64_t  t_req_size = t_d_size + t_d_size_sharedCols + t_d_size_sharedRows + t_d_size_sharedSlices;
                            
#ifdef CUDA_DARTS_DEBUG		
                            std::cout<<"temp arrnGpuEdgeChunk[0]: "<<arrnGpuEdgeChunk[0]<<std::endl;
                            std::cout<<"temp arrnGpuEdgeChunk[1]: "<<arrnGpuEdgeChunk[1]<<std::endl;
                            std::cout<<"temp arrnGpuEdgeChunk[2]: "<<arrnGpuEdgeChunk[2]<<std::endl;
                            std::cout<<"temp d_size: "<<t_d_size<<std::endl;
	    		            std::cout<<"temp d_size_sharedCols: "<<t_d_size_sharedCols<<std::endl;
	    		            std::cout<<"temp d_size_sharedRows: "<<t_d_size_sharedRows<<std::endl;
	    		            std::cout<<"temp d_size_sharedSlices: "<<t_d_size_sharedSlices<<std::endl;
	    		            std::cout<<"temp req size: "<<t_req_size/KB<<"KB,gupMemMax: "<<gpuMemMax/KB<<"KB!"<<std::endl;
#endif
                            
                            if(t_req_size>gpuMemMax){

                                arrnGpuGrid[0] = NUMTILESMAX; 
                                arrnGpuGrid[1] = NUMTILESMAX; 
                                arrnGpuGrid[2] = NUMTILESMAX;

                                arrnGpuEdgeChunk[0]=((arrnGpuGrid[0]*gpuTile_x)>gpuCols  )?gpuCols:  ( arrnGpuGrid[0]*gpuTile_x);
                                arrnGpuEdgeChunk[1]=((arrnGpuGrid[1]*gpuTile_y)>gpuRows  )?gpuRows:  ( arrnGpuGrid[1]*gpuTile_y);
                                arrnGpuEdgeChunk[2]=((arrnGpuGrid[2]*gpuTile_z)>gpuSlices)?gpuSlices:( arrnGpuGrid[2]*gpuTile_z);
                                gnStream = NUMPAR ;

                            }else{
                               gnStream = std::ceil(gpuMemMax/t_req_size)-1;
                               gnStream = (gnStream>NUMPARMAX)?NUMPARMAX:gnStream;
                                    
                            }
                            gnStream = (gnStream<1)?1:gnStream;
                        }

                        vnStream = std::ceil(1.0*arrnGpuEdge[0]/arrnGpuEdgeChunk[0])
                                  *std::ceil(1.0*arrnGpuEdge[1]/arrnGpuEdgeChunk[1])
                                  *std::ceil(1.0*arrnGpuEdge[2]/arrnGpuEdgeChunk[2]);

                        gnStream = (vnStream<gnStream)?vnStream:gnStream; 
                        tnStream= vnStream*nStream; 
                        GpuKernelPureGpuWithStreams37 = Stencil3D7ptGpuKernelPureGpuWithStreamsCD{0,1,this,GPUMETA};
                        add(&GpuKernelPureGpuWithStreams37);
                    
#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"gpu count: "<<nGPU<<std::endl;
                        std::cout<<"gpu vnStream: "<<vnStream<<std::endl;
                        std::cout<<"gpu gnStream: "<<gnStream<<std::endl;
                        std::cout<<"gpu total stream: "<<tnStream<<std::endl;
                        std::cout<<"gpu cutting Grid: x: "<< arrnGpuGrid[0]<<",y: "<<arrnGpuGrid[1]<<", z: "<< arrnGpuGrid[2]<<std::endl;
                        std::cout<<"gpu cutting tile: x: "<< arrnGpuEdgeChunk[0]<<",y: "<<arrnGpuEdgeChunk[1]<<", z: "<<arrnGpuEdgeChunk[2]<<std::endl;

#endif
                    }

#ifdef CUDA_DARTS_DEBUG		
	    		    std::cout<<"gpu : gpuSlices "<<gpuSlices<<std::endl;
	    		    std::cout<<"gpu : gpuRows "<<gpuRows<<std::endl;
	    		    std::cout<<"gpu : gpuCols "<<gpuCols<<std::endl;
#endif
                }else{

                    if(req_size< gpuMemMax){
                       
                            bIsHybrid = false;
                           //may change to CPU when req_size less than the L2 cache,still in validation  
                            nGPU	= 1;
                            gnStream = 1;
						    nStream	= 4;
                            tnStream = nStream;
                            GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                            add(&GpuKernelWithAllTimeSteps37);	

                    }else{

#ifdef CUDA_DARTS_DEBUG		
	    		        std::cout<<"CPU&GPU Hybrid! "<<std::endl;
#endif
                       
                        bIsHybrid = true;
                        
                        if(IsStatic == false){

                            chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
                            chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
                            calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
				            nCPU = arrnCpuGrid[0]*arrnCpuGrid[1]*arrnCpuGrid[2];
				            
                            CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                            for(size_t i=0;i<nCPU; ++i){
	    	   	                CpuLoop37[i] = Stencil3D7ptCpuLoopCD{1,1,this,SHORTWAIT,i};
	    	                }

#ifdef CUDA_DARTS_DEBUG		
	    		            std::cout<<"max cpu count: "<<nCPU<<std::endl;
#endif
                        }


                        nGPU = std::ceil(req_size/gpuMemMax);
                        nStream = NUMASTREAMS;
                        gnStream = NUMPAR ;
                        
                        //int arrnCpuEdgeB[DIM]={200, 200, 200};//cols,rows,slices
                        //int arrnGpuEdgeB[DIM]={100, 100, 100};
                        int arrnCpuEdgeB[DIM]={100, 100, 100};//cols,rows,slices
                        int arrnGpuEdgeB[DIM]={100, 100, 100};

                        chooseSmaller(arrnCpuEdgeInit,arrnCpuEdge,arrnCpuEdgeB,0,0,0,0);
                        chooseSmaller(arrnGpuEdgeInit,arrnGpuEdge,arrnGpuEdgeB,0,0,0,0);
                        
                        arrnCpuEdge[0] = arrnCpuEdgeInit[0] ;
                        arrnCpuEdge[1] = arrnCpuEdgeInit[1] ;
                        arrnCpuEdge[2] = arrnCpuEdgeInit[2] ;

                        arrnGpuEdge[0] = arrnGpuEdgeInit[0] ;
                        arrnGpuEdge[1] = arrnGpuEdgeInit[1] ;
                        arrnGpuEdge[2] = arrnGpuEdgeInit[2] ;

                        nId = std::max_element(arrnEdge,arrnEdge+DIM)-arrnEdge;
                         
#ifdef CUDA_DARTS_DEBUG		
	    		        std::cout<<"max edge: "<<nId<<",where 0:Cols, 1:Rows, 2:Slices "<<std::endl;
#endif

                        //=====================calculate cpu/gpu position=======================================// 
                        int rnId = R(nId);
                        int rnIdM1 = R(nId-1);
                        int rnIdP1 = R(nId+1);
                        arrnCpuPos[rnId]  =0;
                        arrnCpuPos[rnIdM1]=0;
                        arrnCpuPos[rnIdP1]=0;
                        
                        arrnGpuPos[rnId]  = arrnEdge[rnId]-arrnGpuEdge[rnId];
                        arrnGpuPos[rnIdM1]= arrnEdge[rnIdM1]-arrnGpuEdge[rnIdM1];
                        arrnGpuPos[rnIdP1]= arrnEdge[rnIdP1]-arrnGpuEdge[rnIdP1];
                        
#ifdef CUDA_DARTS_DEBUG		
	    		        std::cout<<"cpu Position: ["<< arrnCpuPos[0]<<","<<arrnCpuPos[1]<<","<<arrnCpuPos[2]<<"]"<<std::endl;
	    		        std::cout<<"gpu Position: ["<< arrnGpuPos[0]<<","<<arrnGpuPos[1]<<","<<arrnGpuPos[2]<<"]"<<std::endl;
#endif                   
                        //=====================calculate cpu/gpu position=======================================// 

                        //==================== calc edge left for both CPU and GPU =========================// 
                        calcarrnCpuEdgeLeft(arrnCpuEdgeLeft,arrnEdge,arrnCpuPos, arrnCpuEdge);
                        calcarrnGpuEdgeLeft(arrnGpuEdgeLeft,arrnEdge,arrnGpuPos, arrnGpuEdge);
                        int edgeLeft = arrnGpuPos[R(nId)]-arrnCpuPos[R(nId)]-arrnCpuEdge[R(nId)];
                        arrnCpuEdgeLeft[R(nId)] = edgeLeft ;
                        arrnGpuEdgeLeft[R(nId)] = edgeLeft ;
                        
#ifdef CUDA_DARTS_DEBUG
                        for(int i= 0;i<DIM; ++i){
                            std::cout<<"arrnCpuEdgeLeft["<<i<<"] = "<<arrnCpuEdgeLeft[i]<<std::endl;
                        }

                        for(int i= 0;i<DIM; ++i){
                            std::cout<<"arrnGpuEdgeLeft["<<i<<"] = "<<arrnGpuEdgeLeft[i]<<std::endl;
                        }
#endif
                        //==================== calc edge left for both CPU and GPU =========================// 
                        
                        //==================== calc edge change for both CPU and GPU =========================// 
                        int edgeVar[DIM];
                        edgeVar[0] = INCX;
                        edgeVar[1] = INCY;
                        edgeVar[2] = INCZ;

                        setarrn1Fromarrn2(arrnCpuEdgeVar,edgeVar);
                        setarrn1Fromarrn2(arrnGpuEdgeVar,edgeVar);
                        
                        //==================== calc edge change for both CPU and GPU =========================// 
                        
                        
                        //====================allocate GPU memory========================================// 

                        cudaError err0,err1;

                        arrnGpuEdgeChunk[0] = arrnGpuTile[0]*NUMTILES;
                        arrnGpuEdgeChunk[1] = arrnGpuTile[1]*NUMTILES;
                        arrnGpuEdgeChunk[2] = arrnGpuTile[2]*NUMTILES;
                        
#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"arrnGpuTile["<<0<<"]="<<arrnGpuTile[0]<<std::endl;                       
                        std::cout<<"arrnGpuTile["<<1<<"]="<<arrnGpuTile[1]<<std::endl;                       
                        std::cout<<"arrnGpuTile["<<2<<"]="<<arrnGpuTile[2]<<std::endl;                       

                        std::cout<<"arrnGpuEdgeChunk["<<0<<"]="<<arrnGpuEdgeChunk[0]<<std::endl;                       
                        std::cout<<"arrnGpuEdgeChunk["<<1<<"]="<<arrnGpuEdgeChunk[1]<<std::endl;                       
                        std::cout<<"arrnGpuEdgeChunk["<<2<<"]="<<arrnGpuEdgeChunk[2]<<std::endl;                       
#endif 
                        arrnGpuEdgeChunkAlloc[0]=NUMCHUNKALLOCX;
                        arrnGpuEdgeChunkAlloc[1]=NUMCHUNKALLOCY;
                        arrnGpuEdgeChunkAlloc[2]=NUMCHUNKALLOCZ;

                        const int gpuColsAlloc     = NUMCHUNKALLOCX;
                        const int gpuRowsAlloc     = NUMCHUNKALLOCY;
                        const int gpuSlicesAlloc   = NUMCHUNKALLOCZ;

                        //cudaExtent{width,height,depth} width=gpuTileCols,height=gpuTileRows,depth=gpuTileSlices
                        extent = make_cudaExtent(gpuColsAlloc*szdb,gpuRowsAlloc, gpuSlicesAlloc*numPar);
                        err0=cudaMalloc3D(&devPitchedPtr, extent);
                        dptr_dev = devPitchedPtr.ptr; 
#ifdef CUDA_ERROR_CHECKING
                        gpuErrchk(err0);
#endif
	                    dev_xpitch = devPitchedPtr.pitch/szdb ;
                        dev_ypitch = extent.height;
                        dev_zpitch = extent.depth/numPar;
                        dev_size = dev_xpitch * dev_ypitch * dev_zpitch;

#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"dev size unit: "<<dev_size<<std::endl;
                        std::cout<<"dev_xpitch:"<<dev_xpitch<<std::endl;
                        std::cout<<"dev_ypitch:"<<dev_ypitch<<std::endl;
                        std::cout<<"dev_zpitch:"<<dev_zpitch<<std::endl;
#endif

#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"extend: width: "<<extent.width<<",heigh: "<<extent.height<<", depth: "<<extent.depth<<std::endl;
                        std::cout<<"devPitchedPtr: pointer: "<<devPitchedPtr.ptr<<",pitch: "<<devPitchedPtr.pitch<<",xsize: "<<devPitchedPtr.xsize<<", ysize: "<<devPitchedPtr.ysize<<std::endl;
#endif

	                    int gpuGridDimxAlloc = std::ceil(1.0*gpuColsAlloc/gpuTile_x);
	                    int gpuGridDimyAlloc = std::ceil(1.0*gpuRowsAlloc/gpuTile_y);
	                    int gpuGridDimzAlloc = std::ceil(1.0*gpuSlicesAlloc/gpuTile_z);

                        srows_size   = gpuColsAlloc*gpuGridDimyAlloc*gpuSlicesAlloc*2;
	                    scols_size   = gpuRowsAlloc*gpuGridDimxAlloc*gpuSlicesAlloc*2;
                        sslices_size = gpuColsAlloc*gpuRowsAlloc*gpuGridDimzAlloc*2;

                        dev_sharedCols_size     = numPar*szdb*scols_size;    	
                        dev_sharedRows_size     = numPar*szdb*srows_size;    	
                        dev_sharedSlices_size   = numPar*szdb*sslices_size;	


	                    err0 = cudaMalloc(  &dptr_dev_sharedCols, (dev_sharedCols_size));
#ifdef CUDA_ERROR_CHECKING
	                    gpuErrchk(err0);
#endif
                        err0 = cudaMalloc( (void **) &dptr_dev_sharedRows, (dev_sharedRows_size));
#ifdef CUDA_ERROR_CHECKING
	                    gpuErrchk(err0);
#endif
                        err0 = cudaMalloc( (void **) &dptr_dev_sharedSlices, (dev_sharedSlices_size));
#ifdef CUDA_ERROR_CHECKING
	                    gpuErrchk(err0);
#endif

#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"device dst   size :"<< (dev_size*szdb*numPar)<<std::endl;
                        std::cout<<"sharedCols   size :"<< (dev_sharedCols_size)<<std::endl;
                        std::cout<<"sharedRows   size :"<< (dev_sharedRows_size)<<std::endl;
                        std::cout<<"sharedSlices size :"<< (dev_sharedSlices_size)<<std::endl;
                        std::cout<<"total size        :"<< ((dev_size*szdb*numPar)+ (dev_sharedCols_size)+(dev_sharedRows_size)+(dev_sharedSlices_size))<<std::endl;
#endif
                        //====================allocate GPU memory========================================// 
                        


                        //=====================add cpu/gpu codelet=======================================// 
                        GpuKernelHybridWithStreams37 = Stencil3D7ptGpuKernelHybridWithStreamsCD{0,1,this,GPUMETA};
                        add(&GpuKernelHybridWithStreams37);
                        
                        chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
                        chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);
                        calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 
				        nCPU = arrnCpuGrid[0]*arrnCpuGrid[1]*arrnCpuGrid[2];
                        if(IsStatic == true){
                            CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                            for(size_t i=0;i<nCPU; ++i){
	    	   	                CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
                                add(&CpuLoop37[i]);
	    	                }
                        }else{
                            for(size_t i=0;i<nCPU; ++i){
	    	   	                CDSYNC(CpuLoop37[i]);
                            }
                        }

                        CpuSync37 = Stencil3D7ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT}; 
                        
                        Swap37 = Stencil3D7ptSwapCD{2,2,this,SHORTWAIT};

#ifdef CUDA_DARTS_DEBUG
                        std::cout<<"initial cpuLoop37[0] counter: "<<CpuLoop37[0].getCounter()<<std::endl;
                        std::cout<<"initial cpuSync37 counter: "<<CpuSync37.getCounter()<<std::endl;
#endif
                        //=====================add cpu/gpu codelet=======================================// 

                        //=====================get cpu/gpu start time====================================// 
#ifdef TIMERECORD 
                        uint64_t INF = std::numeric_limits<uint64_t>::max();
                        for(int k=0;k<NUMREC;++k){
                            scRec[k].tStart = 0;
                            scRec[k].tEnd = INF; 
                            scRec[k].tExe = INF;
                            
                            sgRec[k].tStart = 0;
                            sgRec[k].tEnd = INF; 
                            sgRec[k].tExe = INF;
                        }
                        scRec[0].tStart = getTime();

#endif
                        //=====================get cpu/gpu start time====================================// 

                    }
                }
            }
        
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"nGPU = "<<nGPU<<std::endl;
		std::cout<<"nCPU = "<<nCPU<<std::endl;
		std::cout<<"nRows = "<<nRows<<std::endl;
		std::cout<<"nSlices = "<<nSlices<<std::endl;
       
#endif
	}
	
	virtual ~StencilTP(){
        delete []CpuLoop37;
        
        if(bIsHybrid){

            cudaError err0;
            err0 = cudaFree(dptr_dev);
#ifdef CUDA_ERROR_CHECKING
	        gpuErrchk(err0);
#endif
            err0 = cudaFree(dptr_dev_sharedCols);
#ifdef CUDA_ERROR_CHECKING
	        gpuErrchk(err0);
#endif
            err0 = cudaFree(dptr_dev_sharedRows);
#ifdef CUDA_ERROR_CHECKING
	        gpuErrchk(err0);
#endif
            err0 = cudaFree(dptr_dev_sharedSlices);
#ifdef CUDA_ERROR_CHECKING
	        gpuErrchk(err0);
#endif
        }

#ifdef CUDA_DARTS_DEBUG
		std::cout<<"~StencilTP finish!"<<std::endl;
#endif
	}
};

/*
	StencilTP(double * inimatrix,const uint64_t n_rows,const uint64_t n_cols, const uint64_t n_slices,double * newmatrix,uint64_t ts,bool hard,double GpuRatio, Codelet *up)
	:Initial(inimatrix)
	,nRows(n_rows)
	,nCols(n_cols)
    ,nSlices(n_slices)
	,New(newmatrix)
	,ts(ts)
	,tsInit(ts)
    ,hard(hard)
    ,GpuRatio(GpuRatio)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{

#ifdef CUDA_DARTS_DEBUG
		std::cout<<"invoke TP!"<<std::endl;	
		std::cout<<std::setprecision(18)<<std::endl;
#endif
	    	if(GpuRatio == 0.0){
                nCPU = N_CORES;
	    	    cpuPos = 0;
                cpuWL=tWL;
                CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                Swap37 = Stencil3D7ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
                for(size_t i=0;i<nCPU; ++i){
	    	   	    CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    	   	    add(CpuLoop37 + i);
	    	    }

                Swap37 = Stencil3D7ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 

            }else{
           
	    		int deviceCount;
	    		cudaGetDeviceCount(&deviceCount);
	    		cudaDeviceProp props;
	    		cudaGetDeviceProperties(&props,0);

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu device count: "<<deviceCount<<std::endl;
	    		std::cout<<"shared memory per block: "<<props.sharedMemPerBlock/KB<<"KB"<<std::endl;
	    		std::cout<<"registers per Block: "<<props.regsPerBlock<<std::endl;
	    		std::cout<<"Threads per Block:"<<props.maxThreadsPerBlock<<std::endl;
	    		std::cout<<"shared memory per MP: "<<props.sharedMemPerMultiprocessor/KB<<"KB"<<std::endl;
	    		std::cout<<"Threads per MP:"<<props.maxThreadsPerMultiProcessor<<std::endl;
	    		std::cout<<"MB count:"<<props.multiProcessorCount<<std::endl;
	    		std::cout<<"Global memory: "<<props.totalGlobalMem/MB<<"MB"<<std::endl;
#endif

	    		size_t gpu_mem_total_t = 0;
	    		size_t gpu_mem_avail_t = 0;
	    		size_t gpu_mem_valid_t = 0;
	    		cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	    		gpu_mem_valid_t = gpu_mem_avail_t - XMB;
                
                gpuMemMax =(2*GB)> gpu_mem_valid_t?gpu_mem_avail_t: 2*GB;

	            tile_x = (nCols>GRID_TILE37_X)?GRID_TILE37_X:nCols; //16
	            tile_y = (nRows>GRID_TILE37_Y)?GRID_TILE37_Y:nRows; //16
                tile_z = (nSlices>GRID_TILE37_Z)?GRID_TILE37_Z:nSlices; //100 tile_z +2 < NUM_THREAD

                blockDimx = (nCols-2)> tile_x? tile_x:(nCols-2);
                blockDimy = (nRows-2)> tile_y? tile_y:(nRows-2);
                blockDimz = 1;

	            gridDimx = std::ceil(1.0*(nCols)/blockDimx);
	            gridDimy = std::ceil(1.0*(nRows)/blockDimy);
                gridDimz = std::ceil(1.0*(nSlices)/tile_z);
                d_size = sizeof(double)*nRows*nCols*nSlices;  
                d_size_sharedCols = sizeof(double)*nRows*nSlices*gridDimx*2;
                d_size_sharedRows = sizeof(double)*nCols*nSlices*gridDimy*2;
                d_size_sharedSlices = sizeof(double)*nRows*nCols*gridDimz*2;
                req_size = d_size + d_size_sharedCols + d_size_sharedRows + d_size_sharedSlices;

                gpuWLMax = std::floor(1.0*gpuMemMax/(sizeof(double)*(nRows*nCols+nRows*gridDimx*2+nCols*gridDimy*2+1.0*nRows*nCols*(1/tile_z))));

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpuWLMax: "<<gpuWLMax<<std::endl;
#endif
                if (GpuRatio == 1.0){
                    
                    if(req_size < gpuMemMax){
                        nCPU = 0;
                        nGPU = 1;
                        gpuWL = tWL;
                        cpuWL = 0;
                        gpuPos = 0;

	    				GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                        add(&GpuKernelWithAllTimeSteps37);	
                    }else{
                        
                        nCPU=0;
                        nGPU = std::ceil(req_size/gpuMemMax);
                        gpuWL = tWL;
                        gpuPos=0;
                        invokeStreams = true;
                        stream = new cudaStream_t[nStream];
                        for(int i=0;i<nStream;++i){
                            cudaStreamCreate(&stream[i]);
                        }
                        GpuKernelPureGpuWithStreams37 = Stencil3D7ptGpuKernelPureGpuWithStreamsCD{0,1,this,GPUMETA};
                        add(&GpuKernelPureGpuWithStreams37);
                    }

                }else{

                    if(req_size< gpuMemMax){
                       
                        nCPU = 0;
                        nGPU = 1;
                        gpuWL = tWL;
                        cpuWL = 0;
                        gpuPos = 0;

	    				GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                        add(&GpuKernelWithAllTimeSteps37);	

                    }else{

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"CPU&GPU Hybrid! "<<std::endl;
#endif
                        
                        nCPU = N_CORES-1;
                        nGPU = std::ceil(req_size/gpuMemMax);
                        invokeStreams = true;

                        cpuWLMin=nCPU+1;

                        stream = new cudaStream_t[nStream];
                        for(int i=0;i<nStream;++i){
                            cudaStreamCreate(&stream[i]);
	    		        }
    

	    				uint64_t t1 = gpuWLBase;
                        //uint64_t t1 = gpuWLMax;
                        //uint64_t t1 = gpuWLMax*ggInitR;
                        uint64_t t2= tWL*gpuInitR;
                        //gpuWL = t1;
                        gpuWL = (t1<gpuWLMax)?t1:gpuWLMax;
                        //gpuWL = (t2<t1)?t2:t1;
                        
                        uint64_t t3=cpuWLBase;
                        //uint64_t t3=tWL*cpuInitR;
	    				//uint64_t t3 = gpuWL*cpuInitR;
	    				//uint64_t t3 = gpuWL*cgInitR;
                        uint64_t t4 = tWL-gpuWL+2;
                        cpuWL = ((t4-2)<=t3)?t4:t3 ;
	    				wlLeft = (t4-2<=t3)?0:(tWL-gpuWL-cpuWL+4);
                        gpuPos = 0;
	    				cpuPos = gpuWL-2;
                    
                        GpuKernelHybridWithStreams37  =  Stencil3D7ptGpuKernelHybridWithStreamsCD(0,1,this,GPUMETA);
                        add(&GpuKernelHybridWithStreams37);	
                        
                        CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                        CpuSync37 = Stencil3D7ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT}; 
                        for(size_t i=0;i<nCPU; ++i){
	    	   	            CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    	   	            add(CpuLoop37 + i);
	    	            }

                        Swap37 = Stencil3D7ptSwapCD{2,2,this,SHORTWAIT}; 

                        lastCnt = 0.1 * tWL;             
	    				gpuPosInit = gpuPos;
	    				cpuPosInit = cpuPos;
	    				gpuWLInit = gpuWL;
	    				cpuWLInit = cpuWL;
	    				wlLeftInit = wlLeft;
	    		        nCPUInit = nCPU;
	    		        nGPUInit = nGPU;
      
                    }

                }
                 
            
            }
        
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"nGPU = "<<nGPU<<std::endl;
		std::cout<<"nCPU = "<<nCPU<<std::endl;
		std::cout<<"nRows = "<<nRows<<std::endl;
		std::cout<<"nSlices = "<<nSlices<<std::endl;
		std::cout<<"gpuPos = "<<gpuPos<<std::endl;
		std::cout<<"gpuWL = "<<gpuWL<<std::endl;
		std::cout<<"cpuPos = "<<cpuPos<<std::endl;
		std::cout<<"cpuWL = "<<cpuWL<<std::endl;
		std::cout<<"wlLeft = "<<wlLeft<<std::endl;
       
        std::cout<<"invokeStreams = "<<invokeStreams<<std::endl;
#endif
	}
	
	virtual ~StencilTP(){
		delete []CpuLoop;
        delete []CpuLoop37;

        if(invokeStreams==true){
			for(int i=0;i<nStream;++i){
				cudaStreamDestroy(stream[i]);
			}
            delete [] stream;
        }
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"~StencilTP finish!"<<std::endl;
#endif
	}
};
*/
#endif
