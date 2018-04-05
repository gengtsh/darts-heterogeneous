#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H
//#define CUDA_APT_PER_THREAD_DEFAULT_STREAM
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "conf.h"
#include "stencil.h"
//#include <math.h>
#include <cmath>
#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define N_CORES TOTAL_NUM_CU
using namespace darts;

#define GPUMETA 0x4

DEF_CODELET(Stencil2D4ptGpuInitCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelWithAllTimeStepsCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelPureGpuWithStreamsCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelHybridWithStreamsCD,2,LONGWAIT);

DEF_CODELET_ITER(Stencil2D4ptCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(Stencil2D4ptCpuSyncCD,0,SHORTWAIT);
DEF_CODELET(Stencil2D4ptCpuInvokeGpuCD,0,SHORTWAIT);

DEF_CODELET_ITER(Stencil2D4ptGpuLoopCD,0,SHORTWAIT);

DEF_CODELET(Stencil2D4ptSwapCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);



DEF_CODELET_ITER(Stencil3D7ptCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(Stencil3D7ptCpuSyncCD,0,SHORTWAIT);
DEF_CODELET(Stencil3D7ptSwapCD,0,SHORTWAIT);


DEF_CODELET(Stencil3D7ptGpuKernelWithAllTimeStepsCD,2,LONGWAIT);


DEF_TP(StencilTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
    const uint64_t nSlices;
    double *New;
	uint64_t ts;
    uint64_t tsInit;	

	bool hard;
    double GpuRatio;
	double softGpuRatio=0.0;
	uint32_t nCPU = 0;
	uint32_t nGPU = 0;

	Stencil2D4ptGpuKernelCD GpuKernel;
	Stencil2D4ptGpuKernelWithAllTimeStepsCD GpuKernelWithAllTimeSteps;
    Stencil2D4ptGpuKernelPureGpuWithStreamsCD GpuKernelPureGpuWithStreams;
    Stencil2D4ptGpuKernelHybridWithStreamsCD GpuKernelHybridWithStreams;
	Stencil2D4ptCpuLoopCD *CpuLoop = NULL;
	Stencil2D4ptGpuLoopCD *GpuLoop = NULL;
    Stencil2D4ptCpuSyncCD	CpuSync;

	Stencil2D4ptSwapCD Swap;
	
	Stencil3D7ptCpuLoopCD *CpuLoop37 = NULL;
    //Stencil3D7ptCpuSyncCD	CpuSync37;
	Stencil3D7ptSwapCD Swap37;

	Stencil3D7ptGpuKernelWithAllTimeStepsCD GpuKernelWithAllTimeSteps37;

    SyncCD	sync;
    Codelet *signalUp;
	uint64_t nRowsGpu=0;
	uint64_t nRowsCpu=0;
    uint64_t nSlicesGpu=0;
    uint64_t nSlicesCpu=0;


	int tile_x = 0;
	int tile_y = 0;
    int tile_z = 0;
    
    int blockDimx = 0;
    int blockDimy = 0;
    int blockDimz = 0;

    int gridDimx = 0;
    int gridDimy = 0;
    int gridDimz = 0;

    double d_size = 0;
    int64_t d_size_sharedCols = 0;
    int64_t d_size_sharedRows = 0;
    int64_t d_size_sharedSlices = 0;


	double *d_dst = NULL;
	dim3 dimGrid_hack1;
    double req_size=0;
	double cmCpu = 4;	//cpu(total 32 cores) compute ability
	double cmGpu = 10;	//Gpu compute ability
	bool CpuIvGpu = false; // CPU to invoke GPU;
    bool CpuFinish = false; // GPU check CPU computation finish; 
    bool GpuFinish = false; // CPU check GPU computation finish; 
    uint64_t gpuPos=0;
	uint64_t cpuPos=0;

	uint64_t nRowsLeft=0;
    uint64_t nSlicesLeft=0;

	uint64_t gpuPosInit=0;
	uint64_t cpuPosInit=0;
	
    uint64_t nRowsGpuInit=0;
	uint64_t nRowsCpuInit=0;
	uint64_t nRowsLeftInit=0;
   
	uint64_t nSlicesGpuInit=0;
	uint64_t nSlicesCpuInit=0;
	uint64_t nSlicesLeftInit=0;


	uint32_t nCPUInit = 0;
	uint32_t nGPUInit = 0;
    uint32_t gpuCnt = 0;
    uint32_t cpuCnt = 0;

 	double 	gpuInitR = 0.5; 
	double 	cpuInitR = 0.5;
	double 	ggInitR = 0.5;
    double 	cgInitR = 1.0;
	double 	gpuStepR = 0.2;
	double 	cpuStepR = 0.2;
	uint64_t	nRowsGpuBase =2000;
	uint64_t	nRowsCpuBase = 2000;	
	
    uint64_t lastCnt = 0;
    uint64_t nRowsGpuMax= 0;
    size_t gpuMemMax = 0;
    bool invokeStreams = false;
    int nStream = 4 ;
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
        if(nSlices==1){
            //2D4points
            if(nRows==17000){
	            nRowsGpuBase =2000;
	            nRowsCpuBase =4000;	
            }else if(nRows==19000){
	            nRowsGpuBase =1200;
	            nRowsCpuBase =1200;	
            }else if(nRows==23000){
	            nRowsGpuBase =4000;
	            nRowsCpuBase =5000;	
            }else if(nRows==25000){
	            nRowsGpuBase =1000;
	            nRowsCpuBase =800;	
            }else if(nRows==27000){
	            nRowsGpuBase =1000;
	            nRowsCpuBase =1700;	
            }else if(nRows==33000){
	            nRowsGpuBase =3000;
	            nRowsCpuBase =3000;	
            }else if(nRows==37000){
	            nRowsGpuBase =3000;
	            nRowsCpuBase =2000;	
            }else if(nRows==41000){
	            nRowsGpuBase =2000;
	            nRowsCpuBase =3000;	
            }else if(nRows==45000){
	            nRowsGpuBase =4000;
	            nRowsCpuBase =4000;	
            }else if(nRows==47000){
	            nRowsGpuBase =3000;
	            nRowsCpuBase =3000;	
            }else{
	            nRowsGpuBase =1000;
	            nRowsCpuBase =2000;	
            }

	    	if(GpuRatio ==0.0){
	    	    nCPU = N_CORES;
	    	    cpuPos = 0;
	    	    nRowsCpu = nRows;	
              
	    	   CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
	    	   Swap = Stencil2D4ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
	    	   for(size_t i=0;i<nCPU; ++i){
	    	   	CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    	   	add(CpuLoop + i);
	    	   }

	    	   Swap = Stencil2D4ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 

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
                //gpuMemMax = gpu_mem_valid_t;

                tile_y = GRID_TILE_Y;
                tile_x = NUM_THREADS;
                blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
                gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
                gridDimy = std::ceil(1.0*nRows/tile_y); 

	            d_size_sharedCols = sizeof(double)*nRows*gridDimx*2;
	            d_size_sharedRows = sizeof(double)*nCols*gridDimy*2;
                req_size = sizeof(double)* nRows*nCols + d_size_sharedCols + d_size_sharedRows ;

#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"gpu memory total: "<<gpu_mem_total_t/1024<<"KB"<<std::endl;
	    		std::cout<<"gpu memory available: "<<gpu_mem_avail_t/1024<<"KB"<<std::endl;
	    		std::cout<<"required memory size: "<<req_size/1024<<"KB"<<std::endl;
                std::cout<<"2GB:"<<2*GB<<std::endl;
	    		std::cout<<"gpu memory size limition: "<<gpuMemMax/GB<<"GB"<<std::endl;
                if(req_size > 2*gpuMemMax){
	    			std::cout<<"required memory size is larger than 2*gpu_mem_avail_t!"<<std::endl;
	    		}
#endif		
	    		
                nRowsGpuMax =floor((gpuMemMax-2*nStream*(nCols+NUM_THREADS)*4) /(sizeof(double)*(nCols+gridDimx*2 + nCols*2/tile_y )));
#ifdef CUDA_DARTS_DEBUG		
	    		std::cout<<"nRowsGpuMax: "<<nRowsGpuMax<<std::endl;
#endif
                if (GpuRatio == 1.0){
                    nCPU=0;
                    nGPU=std::ceil(req_size/gpuMemMax);
                    nRowsGpu=nRows;
                    nRowsCpu=0;
                    gpuPos=0;
                    invokeStreams = true;       
                    stream = new cudaStream_t[nStream];
	    		    for(int i=0;i<nStream;++i){
	    			    cudaStreamCreate(&stream[i]);
	    		    }
                    GpuKernelPureGpuWithStreams = Stencil2D4ptGpuKernelPureGpuWithStreamsCD{0,1,this,GPUMETA};
                    add(&GpuKernelPureGpuWithStreams);

	    		}else{
	    			if (hard == true){
	    				
	    			}else{
	    			    if(req_size <  gpuMemMax){
	    					nCPU = 0;
	    				    nGPU = 1;
                            nRowsGpu = nRows;

	    			        nRowsCpu = 0;
	    					gpuPos = 0;
	    				
	    				    GpuKernelWithAllTimeSteps = Stencil2D4ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                            add(&GpuKernelWithAllTimeSteps);	
                        }else{
	    					nCPU = N_CORES-1;
                            nGPU = req_size/gpuMemMax +1;
                            invokeStreams = true;
	    				    stream = new cudaStream_t[nStream];
	    		            for(int i=0;i<nStream;++i){
	    			            cudaStreamCreate(&stream[i]);
	    		            }
                            
                            CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
	    			
	    					uint64_t t1 = nRowsGpuBase;
                            //uint64_t t1 = nRowsGpuMax;
                            //uint64_t t1 = nRowsGpuMax*ggInitR;
                            uint64_t t2= nRows*gpuInitR;
                            //nRowsGpu = t1;
                            nRowsGpu = (t1<gpuMemMax)?t1:gpuMemMax;
                            //nRowsGpu = (t2<t1)?t2:t1;
                            
                            uint64_t t3=nRowsCpuBase;
                            //uint64_t t3=nRows*cpuInitR;
	    					//uint64_t t3 = nRowsGpu*cpuInitR;
	    					//uint64_t t3 = nRowsGpu*cgInitR;
                            uint64_t t4 = nRows-nRowsGpu+2;
                            nRowsCpu = (t4<=t3)?t4:t3 ;
	    					nRowsLeft = (t4<=t3)?0:(nRows-nRowsGpu-nRowsCpu+4);
                            gpuPos = 0;
	    					cpuPos = nRowsGpu-2;
                            
                            GpuKernelHybridWithStreams  =  Stencil2D4ptGpuKernelHybridWithStreamsCD(0,1,this,GPUMETA);
                            add(&GpuKernelHybridWithStreams);	
                            
                            CpuSync = Stencil2D4ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT};

	    			        Swap = Stencil2D4ptSwapCD{2,2,this,SHORTWAIT}; 
                        }
	    	
                        lastCnt = 0.1 * nRows;             
	    				gpuPosInit = gpuPos;
	    				cpuPosInit = cpuPos;
	    				nRowsGpuInit = nRowsGpu;
	    				nRowsCpuInit = nRowsCpu;
	    				nRowsLeftInit = nRowsLeft;
	    		        nCPUInit = nCPU;
	    		        nGPUInit = nGPU;
                        
                        for(size_t i=0;i<nCPU; ++i){
	    					CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
	    					add(CpuLoop + i);
	    				}
	    			}

	    		}	

	    	}
        }else{ //3D7Points

	    	if(GpuRatio == 0.0){
                nCPU = N_CORES;
	    	    cpuPos = 0;
                nSlicesCpu=nSlices;
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

	            gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
	            gridDimy = std::ceil(1.0*(nRows-2)/blockDimy);
                gridDimz = std::ceil(1.0*(nSlices-2)/tile_z);
                d_size = sizeof(double)*nRows*nCols*nSlices;  
                d_size_sharedCols = sizeof(double)*nRows*nSlices*gridDimx*2;
                d_size_sharedRows = sizeof(double)*nCols*nSlices*gridDimy*2;
                d_size_sharedSlices = sizeof(double)*nRows*nCols*gridDimz*2;

                req_size = d_size + d_size_sharedCols + d_size_sharedRows + d_size_sharedSlices;

                if (GpuRatio == 1.0){
                    
                    if(req_size < gpuMemMax){
                        nCPU = 0;
                        nGPU = 1;
                        nSlicesGpu = nSlices;
                        nSlicesCpu = 0;
                        gpuPos = 0;

	    				GpuKernelWithAllTimeSteps37 = Stencil3D7ptGpuKernelWithAllTimeStepsCD{0,1,this,GPUMETA};
                        add(&GpuKernelWithAllTimeSteps37);	
                    }else{


                    }


                }else{


                }
                 
            
            }
        }
        
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"nGPU = "<<nGPU<<std::endl;
		std::cout<<"nCPU = "<<nCPU<<std::endl;
		std::cout<<"nRows = "<<nRows<<std::endl;
		std::cout<<"nSlices = "<<nSlices<<std::endl;
		std::cout<<"gpuPos = "<<gpuPos<<std::endl;
		std::cout<<"nRowsGpu = "<<nRowsGpu<<std::endl;
		std::cout<<"cpuPos = "<<cpuPos<<std::endl;
		std::cout<<"nRowsCpu = "<<nRowsCpu<<std::endl;
		std::cout<<"nRowsLeft = "<<nRowsLeft<<std::endl;
       
		std::cout<<"nSlicesLeft = "<<nSlicesLeft<<std::endl;
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

#endif
