#include "LagrangeLeapFrogTP.h"
#include "LagrangeNodalTP.h"
#include "LagrangeElementsTP.h"
#include "CalcTimeConstraintsForElemsTP.h"

static int Iteration=0;
void LagrangeNodalCD::fire(void)
{

	std::cout<<"LagrangeNodalTP begin in iteration: "<<Iteration<<std::endl;
	LOAD_FRAME(LagrangeLeapFrogTP);
	
	Domain *domain=FRAME(domain);
	
	INVOKE(LagrangeNodalTP, domain,&FRAME(lagrangeElements)) ;
	
	EXIT_TP();

}


void LagrangeElementsCD::fire(void)
{
	std::cout<<"LagrangeElements begin in iteration: "<<Iteration<<std::endl;
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	INVOKE(LagrangeElementsTP,domain,domain->numElem(),&FRAME(lagrangeCalcTimeConstraintsForElems));

	EXIT_TP();
}

void LagrangeCalcTimeConstraintsForElemsCD::fire(void)
{
	
	std::cout<<"LagrangeCalcTimeConstraintsForElemsTP begin: "<<Iteration<<std::endl;
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	INVOKE(CalcTimeConstraintsForElemsTP,domain,&FRAME(sync));
	EXIT_TP();
}
void LeapFrogSyncCD::fire(void)
{
	std::cout<<"LagrangeLeapFrogTPSync iteration: "<<++Iteration<<std::endl;
	LOAD_FRAME(LagrangeLeapFrogTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
