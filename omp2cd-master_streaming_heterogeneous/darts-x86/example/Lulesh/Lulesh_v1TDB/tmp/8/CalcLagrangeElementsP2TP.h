#ifndef DARTS_CALCLAGRANGEELEMENTSP2TP_H
#define DARTS_CALCLAGRANGEELEMENTSP2TP_H
 

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcLagrangeElementsP2IterCD,0,SHORTWAIT);

DEF_TP(CalcLagrangeElementsP2TP)
{
	Domain *domain;
	Real_t *vnew;	
	CalcLagrangeElementsP2IterCD *calcLagrangeElementsP2Iter;
	Codelet *signalUp;
	
	CalcLagrangeElementsP2TP(Domain *domain,Real_t *vnew,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,calcLagrangeElementsP2Iter(new CalcLagrangeElementsP2IterCD[N_CORES])
		,signalUp(up)
		{
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////				calcLagrangeElementsP2Iter[i]= CalcLagrangeElementsP2IterCD{0,1,this,SHORTWAIT,i};
////				add ( calcLagrangeElementsP2Iter+ i);
////			}

			calcLagrangeElementsP2Iter[0]= CalcLagrangeElementsP2IterCD{0,1,this,SHORTWAIT,0};
			add ( calcLagrangeElementsP2Iter+ 0);
			if(N_CORES>1){
				calcLagrangeElementsP2Iter[1]= CalcLagrangeElementsP2IterCD{0,1,this,SHORTWAIT,1};
				add ( calcLagrangeElementsP2Iter+ 1);
			}
		}
	virtual ~CalcLagrangeElementsP2TP(){
		delete [] calcLagrangeElementsP2Iter;
	}

};


#endif
