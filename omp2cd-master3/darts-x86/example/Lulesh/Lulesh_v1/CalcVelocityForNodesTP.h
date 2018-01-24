#ifndef DARTS_CALCVELOCITYFORNODESTP_H
#define DARTS_CALCVELOCITYFORNODESTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcVelocityForNodesIterCD,0,SHORTWAIT);

DEF_TP(CalcVelocityForNodesTP)
{
	Domain *domain;
	
	CalcVelocityForNodesIterCD *calcVelocityForNodesIter;
	Codelet *signalUp;
	
	CalcVelocityForNodesTP(Domain *domain,Codelet *up)
		:domain(domain)
		,calcVelocityForNodesIter(new CalcVelocityForNodesIterCD[N_CORES])
		,signalUp(up)
		{
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////				calcVelocityForNodesIter[i]= CalcVelocityForNodesIterCD{0,1,this,SHORTWAIT,i};
////				add ( calcVelocityForNodesIter+ i);
////			}

			calcVelocityForNodesIter[0]= CalcVelocityForNodesIterCD{0,1,this,SHORTWAIT,0};
			add ( calcVelocityForNodesIter+ 0);
			if(N_CORES>1){
				calcVelocityForNodesIter[1]= CalcVelocityForNodesIterCD{0,1,this,SHORTWAIT,1};
				add ( calcVelocityForNodesIter+ 1);
			}
		}
	virtual ~CalcVelocityForNodesTP(){
		delete [] calcVelocityForNodesIter;
	}

};


#endif
