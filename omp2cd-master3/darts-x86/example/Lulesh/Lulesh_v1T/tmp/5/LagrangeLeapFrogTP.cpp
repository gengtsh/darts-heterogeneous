#include "LagrangeLeapFrogTP.h"
#include "LagrangeNodalTP.h"
#include "LagrangeElementsTP.h"
#include "CalcTimeConstraintsForElemsTP.h"

void LagrangeNodalCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	
	Domain *domain=FRAME(domain);
	
	INVOKE(LagrangeNodalTP, domain,&FRAME(lagrangeElements)) ;
	
	EXIT_TP();

}


void LagrangeElementsCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	INVOKE(LagrangeElementsTP,domain,domain->numElem(),&FRAME(lagrangeCalcTimeConstraintsForElems));

	EXIT_TP();
}

void LagrangeCalcTimeConstraintsForElemsCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	INVOKE(CalcTimeConstraintsForElemsTP,domain,&FRAME(sync));
	EXIT_TP();
}
void LeapFrogSyncCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
