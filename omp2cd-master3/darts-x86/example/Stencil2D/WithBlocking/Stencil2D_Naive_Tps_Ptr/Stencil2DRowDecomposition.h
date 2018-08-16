#ifndef DARTS_SPENDIAL2D_ROW_SPLIT_H
#define DARTS_SPENDIAL2D_ROW_SPLIT_H

#include <stdint.h>
#include "DARTS.h"
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

using namespace darts;

DEF_CODELET_ITER(Stencil2DRowLoop,0,SHORTWAIT);
DEF_CODELET(Stencil2DRowSync,2,LONGWAIT);
//DEF_CODELET(Stencil2DSyncToMaster,1,LONGWAIT);

DEF_TP(Stencil2DRowDecomposition)
{
	uint32_t nTp;
	double *Initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	double *New;
    uint64_t *timeStep;
	Stencil2DRowLoop *compute;
    Stencil2DRowSync  sync;
	//Stencil2DSyncToMaster syncToMaster;
	//bool *reset;	
	//uint64_t *counter;
	uint32_t nCPU;
	Codelet*   signalUp;

    int arrnCpuPos[DIM];
    int arrnCpuEdge[DIM];
    int arrnCpuTile[DIM];
    int arrnCpuBlock[DIM];
    int arrnCpuGrid[DIM];
	int arrnCpuGridTileBase[DIM];

    int arrnEdge[DIM];


	Stencil2DRowDecomposition(uint32_t ntp, double *inimatrix,const uint64_t inim,const uint64_t inin,double *newmatrix,uint64_t *ts,Codelet *up, Codelet *cpt=0)
	:nTp(ntp)
	,Initial(inimatrix)
	,nRows(inim)
	,nCols(inin)
	,New(newmatrix)
	,timeStep(ts)
	//,compute(static_cast<Stencil2DRowLoop*>(cpt?cpt:new Stencil2DRowLoop[g_nCU]))
	//,sync(g_nCU,g_nCU,this,LONGWAIT)
	,signalUp(up)
	{

        arrnEdge[0] = nCols;
        arrnEdge[1] = nRows;

	    arrnCpuGridTileBase[0] = GRID_TILECPU_X;
	    arrnCpuGridTileBase[1] = GRID_TILECPU_Y;

        arrnCpuEdge[0]=arrnEdge[0];
        arrnCpuEdge[1]=arrnEdge[1];
        
        setarrnValue(arrnCpuPos,0); 
   

        chooseSmaller(arrnCpuTile,arrnCpuEdge,arrnCpuGridTileBase,0,0,0,0);
        chooseSmaller(arrnCpuBlock,arrnCpuEdge,arrnCpuTile,-2,0,0,-2);

        calcarrnDivCeil(arrnCpuGrid,arrnCpuEdge,arrnCpuBlock); 

        int cpuGridDimx = arrnCpuGrid[0]; 
		int cpuGridDimy = arrnCpuGrid[1]; 
        
		nCPU = cpuGridDimx*cpuGridDimy;

	    compute =static_cast<Stencil2DRowLoop*>(cpt?cpt:new Stencil2DRowLoop[nCPU]);
		sync = Stencil2DRowSync{nCPU,nCPU,this,SHORTWAIT}; 
		for ( size_t i = 0; i < nCPU; ++i ) {
            compute[i] = Stencil2DRowLoop {0,1,this,SHORTWAIT,i};
            add( compute + i );
        }

	}


	~Stencil2DRowDecomposition(){if ((*timeStep)==0) delete []compute;}
};

#endif
