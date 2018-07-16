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
	Index_t numReg;
	Real_t *vnew;
	Real_t ptiny;
	
	CalcMonotonicQRegionForElemsIterCD *calcMonotonicQRegionForElemsIter;
	CalcMonotonicQRegionForElemsSyncCD calcMonotonicQRegionForElemsSync;
	Codelet *signalUp;
	
	CalcMonotonicQRegionForElemsTP(Domain *domain,Index_t numReg,Real_t *vnew,Real_t ptiny,Codelet *up)
		:domain(domain)
		,numReg(numReg)
		,vnew(vnew)
		,ptiny(ptiny)
		,calcMonotonicQRegionForElemsIter(new CalcMonotonicQRegionForElemsIterCD[N_CORES])
		,calcMonotonicQRegionForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			for ( size_t i = 0; i < N_CORES; ++i ) {
			calcMonotonicQRegionForElemsIter[i]= CalcMonotonicQRegionForElemsIterCD{0,1,this,SHORTWAIT,i};
				add ( calcMonotonicQRegionForElemsIter+ i);
			}
			add (&calcMonotonicQRegionForElemsSync);
		}
	virtual ~CalcMonotonicQRegionForElemsTP(){
		delete [] calcMonotonicQRegionForElemsIter;
	}

};


#endif
