#ifndef DARTS_CALCACCELERATIONFORNODESTP_H
#define DARTS_CALCACCELERATIONFORNODESTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

DEF_CODELET_ITER(CalcAccelerationForNodesIterCD,0,SHORTWAIT);

DEF_TP(CalcAccelerationForNodesTP)
{
	Domain *domain;
	Index_t numNode;	
	CalcAccelerationForNodesIterCD *calcAccelerationForNodesIter;
	Codelet *signalUp;
	
	CalcAccelerationForNodesTP(Domain *domain,Index_t numNode,Codelet *up)
		:domain(domain)
		,numNode(numNode)
		,calcAccelerationForNodesIter(new  CalcAccelerationForNodesIterCD[N_CORES])
		,signalUp(up)
		{
	////		for ( size_t i = 0; i < N_CORES; ++i ) {
	////			calcAccelerationForNodesIter[i]=CalcAccelerationForNodesIterCD{0,1,this,SHORTWAIT,i};
	////			add ( calcAccelerationForNodesIter+ i);
	////		}

			//size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < N_TREE; ++i ) {
				calcAccelerationForNodesIter[i]=CalcAccelerationForNodesIterCD{0,1,this,SHORTWAIT,i};
				add ( calcAccelerationForNodesIter+ i);
			}
		}
	virtual ~CalcAccelerationForNodesTP(){
		delete [] calcAccelerationForNodesIter;
	}

};


#endif
