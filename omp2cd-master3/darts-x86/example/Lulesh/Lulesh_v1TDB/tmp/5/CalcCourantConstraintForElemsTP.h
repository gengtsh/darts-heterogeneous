#ifndef DARTS_CALCCOURANTCONSTRAINTFORELEMSTP_H
#define DARTS_CALCCOURANTCONSTRAINTFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcCourantConstraintForElemsP1IterCD,0,SHORTWAIT);
DEF_CODELET(CalcCourantConstraintForElemsSyncCD,1,SHORTWAIT);

DEF_TP(CalcCourantConstraintForElemsTP)
{
	Domain *domain;
	Index_t regElemSize;
	Index_t *regElemlist;
	Index_t *courant_elem_per_thread;
	Real_t *dtcourant_per_thread;
	
	CalcCourantConstraintForElemsP1IterCD *calcCourantConstraintForElemsP1Iter;
	CalcCourantConstraintForElemsSyncCD calcCourantConstraintForElemsSync;
	Codelet *signalUp;                                                     
	CalcCourantConstraintForElemsTP(Domain *domain,Index_t regElemSize,Index_t *regElemlist, Index_t *courant_elem_per_thread,Real_t *dtcourant_per_thread , Codelet *up)
		:domain(domain)
		,regElemSize(regElemSize)
		,regElemlist(regElemlist)
		,courant_elem_per_thread(courant_elem_per_thread ) 
		,dtcourant_per_thread(dtcourant_per_thread ) 
		,calcCourantConstraintForElemsP1Iter(new CalcCourantConstraintForElemsP1IterCD[N_CORES])
		,calcCourantConstraintForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)	
		,signalUp(up)
		{
			
			for ( size_t i = 0; i < N_CORES; ++i ) {
			calcCourantConstraintForElemsP1Iter[i]= CalcCourantConstraintForElemsP1IterCD{0,1,this,SHORTWAIT,i};
				add ( calcCourantConstraintForElemsP1Iter+ i);
			}                               
			add (&calcCourantConstraintForElemsSync);
		}                                   
	virtual ~CalcCourantConstraintForElemsTP(){
		delete [] calcCourantConstraintForElemsP1Iter;
	}

};


#endif
