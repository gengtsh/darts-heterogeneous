#ifndef DARTS_CALCMONOTONICQREGIONFORELEMSTP_H
#define DARTS_CALCMONOTONICQREGIONFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcMonotonicQRegionForElemsIterCD,0,SHORTWAIT);

DEF_CODELET(CalcMonotonicQRegionForElemsSyncCD,0,SHORTWAIT);

DEF_TP(CalcMonotonicQRegionForElemsTP)
{
	Domain *domain;
	Real_t *vnew;
	Index_t numZone;
	const Real_t ptiny = Real_t(1.e-36) ;
	Index_t numReg;
	Index_t numZoneInit;	
	CalcMonotonicQRegionForElemsIterCD *calcMonotonicQRegionForElemsIter;
	CalcMonotonicQRegionForElemsSyncCD calcMonotonicQRegionForElemsSync;
	Codelet *signalUp;
	
	CalcMonotonicQRegionForElemsTP(Domain *domain,Real_t *vnew,Index_t numZone,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,numZone(numZone)
		,numZoneInit(numZone)
		,calcMonotonicQRegionForElemsIter(new CalcMonotonicQRegionForElemsIterCD[N_CORES])
		,calcMonotonicQRegionForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			
			numReg  = domain->numReg();	
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////				calcMonotonicQRegionForElemsIter[i]= CalcMonotonicQRegionForElemsIterCD{0,1,this,SHORTWAIT,i};
////				add ( calcMonotonicQRegionForElemsIter+ i);
////			}
////			add (&calcMonotonicQRegionForElemsSync);
				
			//size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < N_TREE; ++i ) {
				calcMonotonicQRegionForElemsIter[i]= CalcMonotonicQRegionForElemsIterCD{0,1,this,SHORTWAIT,i};
				add ( calcMonotonicQRegionForElemsIter+ i);
			}
		}
	virtual ~CalcMonotonicQRegionForElemsTP(){
		delete [] calcMonotonicQRegionForElemsIter;
	}

};


#endif
