#define __STDC_FORMAT_MACROS
#include <cstdint>
#include "Stencil2DRowDecomposition.h"
#include "Stencil2DKernel.h"
//#include <cmath>
#include <inttypes.h>

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

	SYNC(syncSwap);

	EXIT_TP();
}

void
Stencil2DRowSyncSwap::fire(void)
{
	LOAD_FRAME(Stencil2DRowDecomposition);

    SWAP_PTR(&FRAME(New),&FRAME(Initial));
    if ( FRAME(timeStep)-- > 0 ) {
		RESET(syncSwap);
        uint64_t nCPU = FRAME(nCPU);	
		for(size_t i = 0; i < nCPU; ++i)
			SYNC(compute[i]);
    } else {
		SIGNAL(signalUp);
    }
	
	EXIT_TP();
}
