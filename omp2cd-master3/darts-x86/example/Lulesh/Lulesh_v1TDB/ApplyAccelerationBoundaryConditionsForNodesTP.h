#ifndef DARTS_APPLYACCELERATIONBOUNDARYCONDITIONSFORNODESTP_H
#define DARTS_APPLYACCELERATIONBOUNDARYCONDITIONSFORNODESTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(ApplyAccelerationBoundaryConditionsForNodesIterCD,0,SHORTWAIT);

DEF_TP(ApplyAccelerationBoundaryConditionsForNodesTP)
{
	Domain *domain;
	
	ApplyAccelerationBoundaryConditionsForNodesIterCD *applyAccelerationBoundaryConditionsForNodesIter;
	Codelet *signalUp;
	
	ApplyAccelerationBoundaryConditionsForNodesTP(Domain *domain,Codelet *up)
		:domain(domain)
		,applyAccelerationBoundaryConditionsForNodesIter(new ApplyAccelerationBoundaryConditionsForNodesIterCD[N_CORES])
		,signalUp(up)
		{
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////			applyAccelerationBoundaryConditionsForNodesIter[i]= ApplyAccelerationBoundaryConditionsForNodesIterCD{0,1,this,SHORTWAIT,i};
////			add ( applyAccelerationBoundaryConditionsForNodesIter+ i);
////			}


			//size_t tree = MIN(g_treeBarrier,N_CORES);
			
			for ( size_t i = 0; i < N_TREE; ++i ) {
				applyAccelerationBoundaryConditionsForNodesIter[i]= ApplyAccelerationBoundaryConditionsForNodesIterCD{0,1,this,SHORTWAIT,i};
				add ( applyAccelerationBoundaryConditionsForNodesIter+ i);
			}
		
		}
	virtual ~ApplyAccelerationBoundaryConditionsForNodesTP(){
		delete [] applyAccelerationBoundaryConditionsForNodesIter;
	}

};


#endif
