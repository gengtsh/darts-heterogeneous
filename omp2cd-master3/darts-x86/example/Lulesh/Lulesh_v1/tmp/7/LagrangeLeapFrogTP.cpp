#include "LagrangeLeapFrogTP.h"
#include "LagrangeNodalTP.h"
#include "LagrangeElementsTP.h"
#include "CalcTimeConstraintsForElemsTP.h"

static int Iteration=0;
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

	std::cout<<"LagrangeNodal finish, Iteration is "<<Iteration<<std::endl;
	EXIT_TP();
}

void LagrangeCalcTimeConstraintsForElemsCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	INVOKE(CalcTimeConstraintsForElemsTP,domain,&FRAME(sync));

	std::cout<<"LagrangeElements finish, Iteration is "<<Iteration<<std::endl;
	EXIT_TP();
}

void LeapFrogSyncCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	SIGNAL(signalUp);
	std::cout<<"LagrangeCalcTimeConstraintsForElems finish, Iteration is "<<++Iteration<<std::endl;
	EXIT_TP();
}
