#include "CalcTimeConstraintsForElemsTP.h"
#include "CalcCourantConstraintForElemsTP.h"
#include "CalcHydroConstraintForElemsTP.h"

void CalcTimeConstraintsForElemsCD::fire(void)
{
	LOAD_FRAME(CalcTimeConstraintsForElemsTP);
	Domain *domain=FRAME(domain);

	Index_t *courant_elem_per_thread =FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread =FRAME(dtcourant_per_thread);
	Index_t *hydro_elem_per_thread=FRAME( hydro_elem_per_thread);
	Real_t *dthydro_per_thread=FRAME( dthydro_per_thread);


	for (Index_t r=0 ; r < domain->numReg() ; ++r) {
	   /* evaluate time constraint */
		INVOKE(CalcCourantConstraintForElemsTP, domain, domain->regElemSize(r),domain->regElemlist(r), courant_elem_per_thread, dtcourant_per_thread,&FRAME(calcTimeConstraintsForElemsSync )) ;
	
	   /* check hydro constraint */
	   INVOKE(CalcHydroConstraintForElemsTP,domain, domain->regElemSize(r),domain->regElemlist(r),hydro_elem_per_thread, dthydro_per_thread,&FRAME(calcTimeConstraintsForElemsSync)) ;
	}
	
	EXIT_TP();
}
void CalcTimeConstraintsForElemsSyncCD::fire(void)
{
	LOAD_FRAME(CalcTimeConstraintsForElemsTP);
	SIGNAL(signalUp);

	EXIT_TP();
}
