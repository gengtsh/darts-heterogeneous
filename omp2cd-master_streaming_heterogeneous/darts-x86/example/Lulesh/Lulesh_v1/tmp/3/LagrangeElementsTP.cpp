#include "LagrangeElementsTP.h"
#include "CalcLagrangeElementsTP.h"
#include "CalcQForElemsTP.h"

void CalcLagrangeElementsCD::fire(void)
{
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	
	INVOKE(CalcLagrangeElementsTP, domain,vnew,&FRAME(calcQForElems)) ;

	EXIT_TP();
}

void CalcQForElemsCD::fire(void)
{
	
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	
	INVOKE(CalcQForElemsTP, domain,vnew,&FRAME(lagrangeElementsSync )) ;

	EXIT_TP();
}
void LagrangeElementsSyncCD::fire(void)
{
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Index_t numElem=domain->numElem();
	Real_t *vnew = FRAME(vnew) ; 
	
	//CalcLagrangeElements(*domain, vnew) ;
	
	/* Calculate Q.  (Monotonic q option requires communication) */
//	CalcQForElems(*domain, vnew) ;
	
	ApplyMaterialPropertiesForElems(*domain, vnew) ;
	
	UpdateVolumesForElems(*domain, vnew,
	                      domain->v_cut(), numElem) ;
	
	SIGNAL(signalUp);	
	EXIT_TP();
}
