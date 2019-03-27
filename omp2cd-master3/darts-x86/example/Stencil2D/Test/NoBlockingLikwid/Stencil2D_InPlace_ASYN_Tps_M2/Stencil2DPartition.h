#ifndef DARTS_SPENDIAL2D_PARTITION_H
#define DARTS_SPENDIAL2D_PARTITION_H

#include <stdint.h>
#include <stdlib.h>
#include <darts.h>
#include "Stencil2D_main.h"
#include "DARTS.h"
#include "Stencil2DKernel.h"
using namespace darts;


DEF_CODELET_ITER(Stencil2DPartitionCheckTP,1,SHORTWAIT);
DEF_CODELET_ITER(Stencil2DRowCheck,1,SHORTWAIT);

DEF_CODELET(Stencil2DPartitionChunks,0,SHORTWAIT);
DEF_CODELET(Stencil2DPartitionSync,1,LONGWAIT);
DEF_TP(Stencil2DPartition)
{
	double *initial; //matrix pointer initial matrix[M][N]
	const uint64_t nRows; // matrix M row
	const uint64_t nCols; // Matrix N column
	uint64_t timeStep;	
	
	Stencil2DPartitionChunks nChunks;
	Stencil2DPartitionCheckTP *checkTP;
	Stencil2DRowCheck **check;
	Stencil2DPartitionSync sync;
	Codelet*   signalUp;

	Stencil2DPartition(double *inimatrix,const uint64_t n_rows,const uint64_t n_cols,uint64_t ts, Codelet *up)
	:initial(inimatrix)
	,nRows(n_rows)
	,nCols(n_cols)
	,timeStep(ts)
	,nChunks(0,0,this,SHORTWAIT)
	,checkTP(new Stencil2DPartitionCheckTP[g_nSU+1])
	,check(new Stencil2DRowCheck*[g_nSU])
	,sync(g_nSU,g_nSU,this,LONGWAIT)
	,signalUp(up)
	{
		add(&nChunks);
		//memset(check,0,sizeof(Stencil2DRowCheck*)*g_nSU);
		setRef(9999);		
		for(size_t i=0;i<g_nSU;++i){
			checkTP[i] = Stencil2DPartitionCheckTP {2,2,this,SHORTWAIT,i};
			add(checkTP+i);
			check[i] = new Stencil2DRowCheck[2];
		
		}

		checkTP[g_nSU] = Stencil2DPartitionCheckTP {2,2,this,SHORTWAIT,g_nSU};
		add(checkTP+g_nSU);
	}

	virtual	~Stencil2DPartition(){
		
		delete [] checkTP;
		
		for (size_t i=0;i<g_nSU;++i){
			delete [] check[i];
		}
		delete [] check;
	
		//std::cout<<"PartitonFinish!"<<std::endl;	
	}
};

#endif

