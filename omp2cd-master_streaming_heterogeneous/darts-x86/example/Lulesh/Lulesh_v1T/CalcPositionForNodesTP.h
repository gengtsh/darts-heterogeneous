#ifndef DARTS_CALCPOSITIONFORNODESTP_H
#define DARTS_CALCPOSITIONFORNODESTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcPositionForNodesIterCD,0,SHORTWAIT);

DEF_TP(CalcPositionForNodesTP)
{
	Domain *domain;
	
	CalcPositionForNodesIterCD *calcPositionForNodesIter;
	Codelet *signalUp;
	
	CalcPositionForNodesTP(Domain *domain,Codelet *up)
		:domain(domain)
		,calcPositionForNodesIter(new CalcPositionForNodesIterCD[N_CORES])
		,signalUp(up)
		{
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////			calcPositionForNodesIter[i]= CalcPositionForNodesIterCD{0,1,this,SHORTWAIT,i};
////			add ( calcPositionForNodesIter+ i);
////			}

			//size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < N_TREE; ++i ) {
				calcPositionForNodesIter[i]= CalcPositionForNodesIterCD{0,1,this,SHORTWAIT,i};
				add ( calcPositionForNodesIter+ i);
			}
		
		}
	virtual ~CalcPositionForNodesTP(){
		delete [] calcPositionForNodesIter;
	}

};


#endif
