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
    double GpuRatio;
	double softGpuRatio=0.0;
	int32_t nCPU = 0;
	int32_t nGPU = 0;
	
	Stencil3D7ptCpuLoopCD *CpuLoop37 = NULL;
    Stencil3D7ptCpuSyncCD	CpuSync37;
	Stencil3D7ptSwapCD Swap37;

	Stencil3D7ptGpuKernelWithAllTimeStepsCD GpuKernelWithAllTimeSteps37;
    Stencil3D7ptGpuKernelPureGpuWithStreamsCD GpuKernelPureGpuWithStreams37;
//    Stencil3D7ptGpuKernelHybridWithStreamsCD GpuKernelHybridWithStreams37;
    SyncCD	sync;
    Codelet *signalUp;

    int32_t gpuIdx_xyz;
	int32_t gpuSlices=0;
	int32_t gpuRows=0;
	int32_t gpuCols=0;
    int32_t gpuTileCRS[3];//Columns, Rows,Slices

	int32_t cpuSlices=0;
	int32_t cpuRows=0;
	int32_t cpuCols=0;
    
	int gpuSlicesMin=3;
    int cpuSlicesMin=N_CORES;

	int cpuTile_x = 0;
	int cpuTile_y = 0;
    int cpuTile_z = 0;

    int cpuBlockDimx = 0;
    int cpuBlockDimy = 0;
    int cpuBlockDimz = 0;

    int cpuGridDimx = 0;
    int cpuGridDimy = 0;
    int cpuGridDimz = 0;
	
	int gpuTile_x = 0;
	int gpuTile_y = 0;
    int gpuTile_z = 0;
    int gpuTile_xyz[3];

    int gpuBlockDimx = 0;
    int gpuBlockDimy = 0;
    int gpuBlockDimz = 0;

    int gpuGridDimx = 0;
    int gpuGridDimy = 0;
    int gpuGridDimz = 0;
    int gpuTileGridDim_xyz[3]; 

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

	int32_t nSlicesLeft=0;
	int32_t nRowsLeft=0;
	int32_t nColsLeft=0;

	int32_t gpuPosInit=0;
	int32_t cpuPosInit=0;
	

    int32_t gpuSlicesInit=0;
    int32_t gpuRowsInit=0;
    int32_t gpuColsInit=0;

    int32_t cpuSlicesInit=0;
    int32_t cpuRowsInit=0;
    int32_t cpuColsInit=0;


	int32_t nCPUInit = 0;
	int32_t nGPUInit = 0;
    int32_t gpuCnt = 0;
    int32_t cpuCnt = 0;


//    double 	gpuInitR = 0.5; 
//	double 	cpuInitR = 0.5;
//	double 	ggInitR  = 0.5;
//    double 	cgInitR  = 1.0;
//	double 	gpuStepR = 0.2;
//	double 	cpuStepR = 0.2;
//
//	int32_t	gpuSlicesBase = 100;
//	int32_t	gpuRowsBase   = 100;
//	int32_t	gpuColsBase   = 100;
//
//	int32_t	cpuSlicesBase = 100;	
//	int32_t	cpuRowsBase   = 100;
//	int32_t	cpuColsBase   = 100;
    
	int32_t lastCnt = 0;
    int64_t gpuMemMax = 0;
    int64_t gpuAvMem = GPU_AVMEM_SIZE*GB; 

    int nStream = 0 ;
    int gnStream = 0 ; //group number of stream
    int tnStream = 0 ; //total number of stream 
    cudaStream_t *stream ;


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
                //nCPU = N_CORES;
	    	    cpuPos = 0;
                cpuSlices = nSlices;
				cpuRows = nRows;
				cpuCols = nCols;
				cpuTile_x = (cpuCols>GRID_TILECPU_X)?GRID_TILECPU_X:cpuCols;
				cpuTile_y = (cpuRows>GRID_TILECPU_Y)?GRID_TILECPU_Y:cpuRows;
				cpuTile_z = (cpuSlices>GRID_TILECPU_Z)?GRID_TILECPU_Z:cpuSlices;
				
				cpuBlockDimx = (cpuCols-2)>cpuTile_x?cpuTile_x:(cpuCols-2);
				cpuBlockDimy = (cpuRows-2)>cpuTile_y?cpuTile_y:(cpuRows-2);
				cpuBlockDimz = (cpuSlices-2)>cpuTile_z?cpuTile_z:(cpuSlices-2);
				
				cpuGridDimx = std::ceil(1.0*cpuCols/cpuBlockDimx);
				cpuGridDimy = std::ceil(1.0*cpuRows/cpuBlockDimy);
				cpuGridDimz = std::ceil(1.0*cpuSlices/cpuBlockDimz);

				
				nCPU = cpuGridDimx*cpuGridDimy*cpuGridDimz;
				CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
                Swap37 = Stencil3D7ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
                for(size_t i=0;i<nCPU; ++i){
	    	   	    CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    	   	    add(CpuLoop37 + i);
	    	    }

                Swap37 = Stencil3D7ptSwapCD{nCPU,nCPU,this,SHORTWAIT};
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"cpuBlockDimx = "<<cpuBlockDimx<<std::endl;
		std::cout<<"cpuBlockDimy = "<<cpuBlockDimy<<std::endl;
		std::cout<<"cpuBlockDimz = "<<cpuBlockDimz<<std::endl;
		std::cout<<"cpuGridDimz = "<<cpuGridDimz<<std::endl;
		std::cout<<"cpuGridDimy = "<<cpuGridDimy<<std::endl;
		std::cout<<"cpuGridDimz = "<<cpuGridDimz<<std::endl;
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

	    	    gpuPos = 0;
                gpuSlices = nSlices;
				gpuRows = nRows;
				gpuCols = nCols;
                gpuTileCRS[0]=gpuCols;                
                gpuTileCRS[1]=gpuRows;                
                gpuTileCRS[2]=gpuSlices;
			   

                gpuTile_x = (gpuCols    >GRID_TILE37_X)?GRID_TILE37_X:gpuCols;
				gpuTile_y = (gpuRows    >GRID_TILE37_Y)?GRID_TILE37_Y:gpuRows;
				gpuTile_z = (gpuSlices  >GRID_TILE37_Z)?GRID_TILE37_Z:gpuSlices;
		        gpuTile_xyz[0]=gpuTile_x;
		        gpuTile_xyz[1]=gpuTile_y;
		        gpuTile_xyz[2]=gpuTile_z;

                
                gpuBlockDimx = (gpuCols-2  )>gpuTile_x?gpuTile_x:(gpuCols-2);
				gpuBlockDimy = (gpuRows-2  )>gpuTile_y?gpuTile_y:(gpuRows-2);
				gpuBlockDimz = 1;

				
				gpuGridDimx = std::ceil(1.0*gpuCols/gpuBlockDimx);
				gpuGridDimy = std::ceil(1.0*gpuRows/gpuBlockDimy);
				gpuGridDimz = std::ceil(1.0*gpuSlices/gpuTile_z);
                
                gpuTileGridDim_xyz[0]=gpuGridDimx;
                gpuTileGridDim_xyz[1]=gpuGridDimy;
                gpuTileGridDim_xyz[2]=gpuGridDimz;

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpuGridDimx: "<<gpuGridDimx<<std::endl;;
	    		std::cout<<"gpuGridDimy: "<<gpuGridDimy<<std::endl;;
	    		std::cout<<"gpuGridDimz: "<<gpuGridDimz<<std::endl;;
#endif


                d_size = sizeof(double)*gpuRows*gpuCols*gpuSlices;  
                d_size_sharedCols   = sizeof(double)*gpuRows*gpuSlices*gpuGridDimx*2;
                d_size_sharedRows   = sizeof(double)*gpuCols*gpuSlices*gpuGridDimy*2;
                d_size_sharedSlices = sizeof(double)*gpuRows*gpuCols  *gpuGridDimz*2;
                req_size = d_size + d_size_sharedCols + d_size_sharedRows + d_size_sharedSlices;
                //gpuWLMax = std::floor(1.0*gpuMemMax/(sizeof(double)*(gpuRows*gpuCols+gpuRows*gpuGridDimx*2+gpuCols*gpuGridDimy*2+1.0*gpuRows*gpuCols*(1/gpuTile_z))));

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"req size: "<<req_size/KB<<"KB,gupMemMax: "<<gpuMemMax/KB<<"KB!"<<std::endl;
#endif
                if (GpuRatio == 1.0){
                    gpuPos      = 0;
					gpuSlices	= nSlices;
					gpuRows		= nRows;
					gpuCols		= nCols;
                    if(req_size < gpuMemMax){
#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu require size is less than gpuMemMax! "<<std::endl;
#endif
                        nGPU	= 1;
                        gnStream = 1;
						nStream	= 4;
                        tnStream = nStream;
						stream = new cudaStream_t[tnStream];
                        for(int i=0;i<tnStream;++i){
                            cudaStreamCreate(&stream[i]);
                        }
                        GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                        add(&GpuKernelWithAllTimeSteps37);	
                    }else{

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu require size is greater than gpuMemMax! "<<std::endl;
#endif
                        nGPU = std::ceil(req_size/gpuMemMax);
                        nStream = 5;
                        gnStream= 4;//may change
                        for (int i=0;i<2;i++){ 
                            int gx = gpuTileGridDim_xyz[0];
                            if(std::all_of(gpuTileGridDim_xyz,gpuTileGridDim_xyz+3,[gx](int x){return x==gx;})){
                            //if((gpuTileGridDim_xyz[0]==gpuTileGridDim_xyz[1] )&&(gpuTileGridDim_xyz[0]==gpuTileGridDim_xyz[2] )){
                                gpuIdx_xyz=2;//cut based on slice
                            }else{
                                gpuIdx_xyz=std::max_element(gpuTileGridDim_xyz,gpuTileGridDim_xyz+3)-gpuTileGridDim_xyz;//find the max among  gridDim of row,cols,slices
                            }
                            
#ifdef CUDA_DARTS_DEBUG
                            std::cout<<"gpu cutting " <<i<<" based on : (0:x(columns); 1:y(rows);2:z(slices)):"<<gpuIdx_xyz<<std::endl;
#endif
                            int ncut = (i==0)?nGPU:gnStream;
                            gpuTileGridDim_xyz[gpuIdx_xyz] =  std::ceil(1.0*gpuTileGridDim_xyz[gpuIdx_xyz]/ncut) ;

#ifdef CUDA_DARTS_DEBUG
                            std::cout<<"gpu cutting "<<i<<" Grid: x: "<< gpuTileGridDim_xyz[0]<<",y: "<<gpuTileGridDim_xyz[1]<<", z: "<< gpuTileGridDim_xyz[2]<<std::endl;
#endif
                        }
                        //for(int i=0;i<3;++i){
                        //    gpuTileCRS[i]=gpuTileGridDim_xyz[i]*gpuTile_xyz[i]; 
                        //}
                        gpuTileCRS[0]=gpuTileGridDim_xyz[0]*gpuTile_x;
                        gpuTileCRS[1]=gpuTileGridDim_xyz[1]*gpuTile_y;
                        gpuTileCRS[2]=gpuTileGridDim_xyz[2]*gpuTile_z;

#ifdef CUDA_DARTS_DEBUG

                        std::cout<<"gpu count: "<<nGPU<<std::endl;
                        std::cout<<"gpu gnStream: "<<gnStream<<std::endl;
                        std::cout<<"gpu cutting Grid: x: "<< gpuTileGridDim_xyz[0]<<",y: "<<gpuTileGridDim_xyz[1]<<", z: "<< gpuTileGridDim_xyz[2]<<std::endl;
                        std::cout<<"gpu cutting tile: x: "<< gpuTileCRS[0]<<",y: "<<gpuTileCRS[1]<<", z: "<<gpuTileCRS[2]<<std::endl;

#endif
                        tnStream= nGPU*gnStream*nStream; 
                        stream = new cudaStream_t[tnStream];
                        for(int i=0;i<tnStream;++i){
                            cudaStreamCreate(&stream[i]);
                        }
                        GpuKernelPureGpuWithStreams37 = Stencil3D7ptGpuKernelPureGpuWithStreamsCD{0,1,this,GPUMETA};
                        add(&GpuKernelPureGpuWithStreams37);
                    
                    }

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu : gpuSlices "<<gpuSlices<<std::endl;
	    		std::cout<<"gpu : gpuRows "<<gpuRows<<std::endl;
	    		std::cout<<"gpu : gpuCols "<<gpuCols<<std::endl;
#endif
                }else{

                    if(req_size< gpuMemMax){
//                       
//                        nCPU = 0;
//                        nGPU = 1;
//                        gpuWL = tWL;
//                        cpuWL = 0;
//                        gpuPos = 0;
//
//	    				GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
//                        add(&GpuKernelWithAllTimeSteps37);	
//
                    }else{
//
//#ifdef CUDA_DARTS_DEBUG		
//	    		std::cout<<"CPU&GPU Hybrid! "<<std::endl;
//#endif
//                        
//                        nCPU = N_CORES-1;
//                        nGPU = std::ceil(req_size/gpuMemMax);
//                        invokeStreams = true;
//
//                        cpuWLMin=nCPU+1;
//
//                        stream = new cudaStream_t[nStream];
//                        for(int i=0;i<nStream;++i){
//                            cudaStreamCreate(&stream[i]);
//	    		        }
//    
//
//	    				uint64_t t1 = gpuWLBase;
//                        //uint64_t t1 = gpuWLMax;
//                        //uint64_t t1 = gpuWLMax*ggInitR;
//                        uint64_t t2= tWL*gpuInitR;
//                        //gpuWL = t1;
//                        gpuWL = (t1<gpuWLMax)?t1:gpuWLMax;
//                        //gpuWL = (t2<t1)?t2:t1;
//                        
//                        uint64_t t3=cpuWLBase;
//                        //uint64_t t3=tWL*cpuInitR;
//	    				//uint64_t t3 = gpuWL*cpuInitR;
//	    				//uint64_t t3 = gpuWL*cgInitR;
//                        uint64_t t4 = tWL-gpuWL+2;
//                        cpuWL = ((t4-2)<=t3)?t4:t3 ;
//	    				wlLeft = (t4-2<=t3)?0:(tWL-gpuWL-cpuWL+4);
//                        gpuPos = 0;
//	    				cpuPos = gpuWL-2;
//                    
//                        GpuKernelHybridWithStreams37  =  Stencil3D7ptGpuKernelHybridWithStreamsCD(0,1,this,GPUMETA);
//                        add(&GpuKernelHybridWithStreams37);	
//                        
//                        CpuLoop37 = new Stencil3D7ptCpuLoopCD[nCPU];
//                        CpuSync37 = Stencil3D7ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT}; 
//                        for(size_t i=0;i<nCPU; ++i){
//	    	   	            CpuLoop37[i] = Stencil3D7ptCpuLoopCD{0,1,this,SHORTWAIT,i};
//	    	   	            add(CpuLoop37 + i);
//	    	            }
//
//                        Swap37 = Stencil3D7ptSwapCD{2,2,this,SHORTWAIT}; 
//
//                        lastCnt = 0.1 * tWL;             
//	    				gpuPosInit = gpuPos;
//	    				cpuPosInit = cpuPos;
//	    				gpuWLInit = gpuWL;
//	    				cpuWLInit = cpuWL;
//	    				wlLeftInit = wlLeft;
//	    		        nCPUInit = nCPU;
//	    		        nGPUInit = nGPU;
//      
                    }
//
                }
                 
            
            }
        
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"nGPU = "<<nGPU<<std::endl;
		std::cout<<"nCPU = "<<nCPU<<std::endl;
		std::cout<<"nRows = "<<nRows<<std::endl;
		std::cout<<"nSlices = "<<nSlices<<std::endl;
		std::cout<<"gpuPos = "<<gpuPos<<std::endl;
		std::cout<<"cpuPos = "<<cpuPos<<std::endl;
       
#endif
	}
	
	virtual ~StencilTP(){
        delete []CpuLoop37;
		for(int i=0;i<tnStream;++i){
			cudaStreamDestroy(stream[i]);
		}
        delete [] stream;
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
