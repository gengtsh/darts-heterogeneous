#ifndef STENCILGPUTP_H
#define STENCILGPUTP_H

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "conf.h"
#include "stencil.h"
#include <math.h>
#include <cstdint>
#include "DARTS.h"

//#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define N_CORES TOTAL_NUM_CU
using namespace darts;

#define GPUMETA 0x4

DEF_CODELET(Stencil2D4ptGpuInitCD,2,LONGWAIT);
DEF_CODELET(Stencil2D4ptGpuKernelCD,2,LONGWAIT);

DEF_CODELET_ITER(Stencil2D4ptCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(Stencil2D4ptCpuSyncCD,0,SHORTWAIT);
DEF_CODELET(Stencil2D4ptCpuInvokeGpuCD,0,SHORTWAIT);

DEF_CODELET_ITER(Stencil2D4ptGpuLoopCD,0,SHORTWAIT);

DEF_CODELET(Stencil2D4ptSwapCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(StencilTP)
{
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
	uint64_t ts;
	bool hard;
	double GpuRatio = 1.0 ;
	double softGpuRatio=0.0;
	uint32_t nCPU = 0;
	uint32_t nGPU = 0;

	Stencil2D4ptGpuKernelCD GpuKernel;
	Stencil2D4ptGpuLoopCD *GpuLoop = NULL;
	Stencil2D4ptCpuLoopCD *CpuLoop = NULL;
	Stencil2D4ptCpuSyncCD	CpuSync;
//	Stencil2D4ptCpuInvokeGpuCD	CpuInvokeGpu;

	Stencil2D4ptSwapCD Swap;
	SyncCD	sync;
    Codelet *signalUp;
	uint64_t nRowsGpu;
	uint64_t nRowsCpu;

	double d_size;
	double *d_dst = NULL;
	double *d_src = NULL;
	dim3 dimGrid_hack1;
	double req_size;
	uint64_t nRowsLeft=0;
	double cmCpu = 3;	//cpu(total 32 cores) compute ability
	double cmGpu = 4;	//Gpu compute ability
	bool CpuIvGpu = false; // CPU to invoke GPU;
	uint64_t gpuPos=0;
	uint64_t cpuPos=0;

	StencilTP(double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t ts,bool hard,Codelet *up)
	:Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,ts(ts)
	,hard(hard)
	,sync(1,1,this,LONGWAIT)
	,signalUp(up)
	{
//		std::cout<<"invoke TP!"<<std::endl;	

		if(GpuRatio ==0.0){
			nCPU = N_CORES;
			CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
			Swap = Stencil2D4ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
			cpuPos = 0;
			nRowsCpu = nRows;	
			for(size_t i=0;i<nCPU; ++i){
				CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
				add(CpuLoop + i);
			}

			Swap = Stencil2D4ptSwapCD{nCPU,nCPU,this,SHORTWAIT}; 
		}else{

			int deviceCount;
			cudaGetDeviceCount(&deviceCount);
//			std::cout<<"gpu device count: "<<deviceCount<<std::endl;
#ifdef CUDA_DARTS_DEBUG		
			cudaDeviceProp props;
			cudaGetDeviceProperties(&props,0);
			int KB=1024;
			int MB=KB*KB;
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
			cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);

#ifdef CUDA_DARTS_DEBUG		
			std::cout<<"gpu memory total: "<<gpu_mem_total_t/1024<<"KB"<<std::endl;
			std::cout<<"gpu memory available: "<<gpu_mem_avail_t/1024<<"KB"<<std::endl;
#endif		
			req_size = sizeof(double)* nRows*nCols*2;
			
			if (GpuRatio == 1.0){
				nCPU = 0;
				nGPU=std::ceil(req_size/gpu_mem_avail_t);
				nRowsGpu = nRows;
				nRowsCpu = 0;
				gpuPos = 0;
				//std::cout<<"gpu loop number:" <<nGPU<<std::endl;
				GpuLoop = new Stencil2D4ptGpuLoopCD[nGPU];
				for(size_t i=0;i<nGPU; ++i){
					GpuLoop[i] = Stencil2D4ptGpuLoopCD{0,1,this,GPUMETA,i};

					add(GpuLoop + i);
				}
			
				Swap = Stencil2D4ptSwapCD{1,1,this,SHORTWAIT}; 
			}else{
				if (hard == true){
					
				}else{
					if(req_size < gpu_mem_avail_t){
						nCPU = 0;
						nRowsGpu = nRows;
						gpuPos = 0;
					}else if (req_size<2 * gpu_mem_avail_t){
						nCPU = N_CORES-1;
						CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
						nRowsGpu = nRows*cmGpu/(cmGpu+cmCpu);
						nRowsCpu = nRows - nRowsGpu ;
						gpuPos = 0;
						cpuPos = nRowsGpu;

						CpuSync = Stencil2D4ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT};
					}else{
						nCPU = N_CORES-1;
						CpuLoop = new Stencil2D4ptCpuLoopCD[nCPU];
						nRowsGpu=floor(gpu_mem_avail_t/(nCols*2*sizeof(double)));
						nRowsCpu = nRowsGpu*(cmCpu/cmGpu);
						nRowsLeft=nRows - nRowsGpu-nRowsCpu;			
//						std::cout<<"nRowsLeft initial:" <<nRowsLeft <<std::endl;
						if (nRowsLeft< (1/cmGpu)*nRowsGpu){
							nRowsCpu = nRows - nRowsGpu;
							nRowsLeft = 0;
						}else{
							CpuIvGpu = true ;
							//CpuInvokeGpu = Stencil2D4ptCpuInvokeGpuCD{1,1,this,GPUMETA};

						}
						gpuPos = 0;
						cpuPos = nRowsGpu;

						CpuSync = Stencil2D4ptCpuSyncCD{nCPU,nCPU,this,SHORTWAIT};

					}
					nGPU = 1;
		
					GpuKernel = Stencil2D4ptGpuKernelCD{0,1,this,GPUMETA};
					add(&GpuKernel);
					for(size_t i=0;i<nCPU; ++i){
						CpuLoop[i] = Stencil2D4ptCpuLoopCD{0,1,this,SHORTWAIT,i};
						add(CpuLoop + i);
					}
				}

				Swap = Stencil2D4ptSwapCD{2,2,this,SHORTWAIT}; 
			}	

		}
#ifdef CUDA_DARTS_DEBUG
		std::cout<<"nGPU = "<<nGPU<<std::endl;
		std::cout<<"nCPU = "<<nCPU<<std::endl;
		std::cout<<"nRows = "<<nRows<<std::endl;
		std::cout<<"gpuPos = "<<gpuPos<<std::endl;
		std::cout<<"nRowsGpu = "<<nRowsGpu<<std::endl;
		std::cout<<"cpuPos = "<<cpuPos<<std::endl;
		std::cout<<"nRowsCpu = "<<nRowsCpu<<std::endl;
		std::cout<<"nRowsLeft = "<<nRowsLeft<<std::endl;
#endif
	}
	
	virtual ~StencilTP(){
		if(GpuRatio == 1.0){
			delete []GpuLoop;
		}else {
			delete []CpuLoop;
		}

	}
};

#endif
