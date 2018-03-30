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

			calcPositionForNodesIter[0]= CalcPositionForNodesIterCD{0,1,this,SHORTWAIT,0};
			add ( calcPositionForNodesIter+ 0);
			if(N_CORES>1){
				calcPositionForNodesIter[1]= CalcPositionForNodesIterCD{0,1,this,SHORTWAIT,1};
				add ( calcPositionForNodesIter+ 1);
			}
		
		}
	virtual ~CalcPositionForNodesTP(){
		delete [] calcPositionForNodesIter;
	}

};


#endif
