#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "Stencil2DRowDecomposition.h"
#include "Stencil2DKernel.h"
#include <math.h>
#include <inttypes.h>
#include <sched.h>

pthread_mutex_t mutex2;

void 
Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	uint64_t Id = getID();
    RESET(compute[Id]);

#ifdef STENCIL2DDEBUG
	std::cout<<"Invoke CPU["<<Id<<"]"<<std::endl;	
#endif

    double        *h_src  = FRAME(Initial);
	double        *h_dst  = FRAME(New);
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
#ifdef STENCIL2DDEBUG

	pthread_mutex_lock(&mutex2);
	std::cout<<"CpuLoop["<<Id<<"]: cpuPos:"<<cpuPos<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: cpuPosx:"<<cpuPosx<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: cpuPosy:"<<cpuPosy<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: posx:"<<posx<<std::endl;
    std::cout<<"CpuLoop["<<Id<<"]: posy:"<<posy<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nRowsChunk:"<<nRowsChunk<<std::endl;
	std::cout<<"CpuLoop["<<Id<<"]: nColsChunk:"<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif

	computeBlock_stencil25(dst,src,nRows,nCols,nRowsChunk,nColsChunk);

	SYNC(sync);

	EXIT_TP();
}


//void 
//Stencil2DRowLoop::fire(void) 
//{
//	LOAD_FRAME(Stencil2DRowDecomposition);
//	uint64_t nCPU = FRAME(nCPU);	
//	double *Initial = FRAME(Initial); //matrix pointer initial Matrix[M][n]
//	const uint64_t n_rows = FRAME(nRows); // matrix M row
//	const uint64_t n_cols = FRAME(nCols); // Matrix N column
//	double *New = FRAME(New);
//	uint64_t bp = n_cols + 1;//block position begin in the [row1,cols1]
//	uint64_t BlockM = n_rows-2;
//	uint64_t BlockN = n_cols-2;
//	uint64_t Id = getID();	
//
//	const uint64_t rows_ini = BlockM / nCPU; // initially,the total number of rows in every nCPU
//	uint64_t BlockM_final = ((Id==(nCPU-1))? (BlockM%nCPU):0) + rows_ini;//blockM from 1 to (nCPU-1), is the same (nr_initial), but the last one is different with the former, it needs add the remains.
//	uint64_t pos = Id*rows_ini * n_cols;//beginning position of every nCPU 
//
//
//	computeInner_stencil2d(bp+pos,Initial,New, BlockM_final,BlockN,n_cols);
//
////	uint64_t block_sz = BlockM_final * BlockN;//block size
////	/**
////	* if block size < TILE_SZ then compute_inner()
////	* else decompose the block into sub_blocks.
////	*/
////	if(block_sz <= TOTAL_TILE_SZ ){
////		computeInner_stencil2d(bp+pos, BlockM_final,BlockN,Initial,New, n_cols);
////	}else{
////		uint64_t M = BlockM_final/N_ROWS_TILE_SZ;
////		uint64_t N = BlockN/N_COLS_TILE_SZ;
////
////		for (size_t i=0;i<= M;i+=N_ROWS_TILE_SZ)
////			for(size_t j=0;j<= N;j+=N_COLS_TILE_SZ){
////				uint64_t sub_m =(i==M)?(BlockM_final-i*N_ROWS_TILE_SZ):( (BlockM_final<N_ROWS_TILE_SZ)? BlockM_final:N_ROWS_TILE_SZ);
////				uint64_t sub_n = (j==N)?(BlockN-j*N_ROWS_TILE_SZ):((BlockN<N_COLS_TILE_SZ)?BlockN:N_COLS_TILE_SZ);
////				uint64_t sub_pos = i*n_cols*N_COLS_TILE_SZ+j*sub_n;
////				computeInner_stencil2d(bp+pos+sub_pos, sub_m,sub_n,Initial,New, n_cols);
////			}
////		
////	}
////
//
////	printf("Line CPU %d: %" PRIu64 " - %" PRIu64 "\n",(int)getID(),(BlockPosition+pos_final),BlockM_final);
//
//	//computeInner_stencil2d(bp+pos_final, BlockM_final,BlockN,Initial,New, n_cols);
//	
////	std::cout<<"core num: "<<sched_getcpu()<<" compute: "<<Id<<"\n"<<std::endl;		
////	RESET(compute[Id]);    
//	SYNC(sync);
//	EXIT_TP();
//}

void
Stencil2DRowSync::fire(void)
{
    LOAD_FRAME(Stencil2DRowDecomposition);
	
//	std::cout<<"core num: "<<sched_getcpu()<<" Tp: "<<FRAME(nTp)<<" RowSync! "<<" \n"  <<std::endl;
//	if(*(FRAME(timeStep))!=1)
//		SYNC(syncToMaster);
	SIGNAL(signalUp);
//	std::cout<<"signalUp: "<<FRAME(signalUp)->getCounter() <<std::endl;
//	RESET(sync);
    EXIT_TP();
}

//void
//Stencil2DSyncToMaster::fire(void)
//{
//
//    LOAD_FRAME(Stencil2DRowDecomposition);
//
//	std::cout<<"core num: "<<sched_getcpu()<< " SyncToMaster TP "<<FRAME(nTp)<<" begin!\n"  <<std::endl;
//	uint64_t nCPU = FRAME(nCPU);	
//	while(1){
//
//		std::cout<<"core num: "<<sched_getcpu()<< " SyncToMaster TP "<<FRAME(nTp)<<" Reset signal: " <<*FRAME(reset)<<"\n"  <<std::endl;
//
//		if((*FRAME(reset))==true){
//			for(size_t i=0;i<nCPU; ++i){
//				SYNC(compute[i]);
//			}
//		//	(*FRAME(reset))=false;
//			(*FRAME(counter))++;
//			if((*FRAME(counter))==g_nSU){
//				(*FRAME(reset))=false;
//				*(FRAME(counter))=0;
//			}
//
//			RESET(syncToMaster);
//			break;
//		}
//	}
//
//	std::cout<<"core num: "<<sched_getcpu()<< " SyncToMaster TP "<<FRAME(nTp)<<" finish!\n"  <<std::endl;
//	EXIT_TP();
//}
//
