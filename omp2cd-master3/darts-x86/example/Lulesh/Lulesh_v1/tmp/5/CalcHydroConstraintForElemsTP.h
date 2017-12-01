#ifndef DARTS_CALCHYDROCONSTRAINTFORELEMSTP_H
#define DARTS_CALCHYDROCONSTRAINTFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcHydroConstraintForElemsP1IterCD,0,SHORTWAIT);
DEF_CODELET(CalcHydroConstraintForElemsSyncCD,1,SHORTWAIT);

DEF_TP(CalcHydroConstraintForElemsTP)
{
	Domain *domain;
	Index_t regElemSize;
	Index_t *regElemlist;
	Index_t *hydro_elem_per_thread;
	Real_t *dthydro_per_thread;
	CalcHydroConstraintForElemsP1IterCD *calcHydroConstraintForElemsP1Iter;
	CalcHydroConstraintForElemsSyncCD calcHydroConstraintForElemsSync;
	Codelet *signalUp;                                                     
	CalcHydroConstraintForElemsTP(Domain *domain,Index_t regElemSize,Index_t *regElemlist, Index_t *hydro_elem_per_thread,Real_t *dthydro_per_thread, Codelet *up)
		:domain(domain)
		,regElemSize(regElemSize)
		,regElemlist(regElemlist)
		,hydro_elem_per_thread(hydro_elem_per_thread )
		,dthydro_per_thread( dthydro_per_thread)
		,calcHydroConstraintForElemsP1Iter(new CalcHydroConstraintForElemsP1IterCD[N_CORES])
		,calcHydroConstraintForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)	
		,signalUp(up)
		{
			
			for ( size_t i = 0; i < N_CORES; ++i ) {
			calcHydroConstraintForElemsP1Iter[i]= CalcHydroConstraintForElemsP1IterCD{0,1,this,SHORTWAIT,i};
				add ( calcHydroConstraintForElemsP1Iter+ i);
			}                               
			add (&calcHydroConstraintForElemsSync);
		}                                   
	virtual ~CalcHydroConstraintForElemsTP(){
		delete [] calcHydroConstraintForElemsP1Iter;
	}

};


#endif
