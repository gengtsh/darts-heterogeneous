#include "lulesh.h"
#include "LagrangeLeapFrogTP.h"

//void NodalCD::fire(void)
//{
//	LOAD_FRAME(LagrangeLeapFrogTP);
//	Domain *domain=FRAME(domain);
//	//LagrangeNodal(*domain);
//	
//	const Real_t delt = domain->deltatime() ;
//	Real_t u_cut = domain->u_cut() ;
//	
//	/* time of boundary condition evaluation is beginning of step for force and
//	* acceleration boundary conditions. */
//	CalcForceForNodes(*domain);
//	CalcAccelerationForNodes(*domain, domain->numNode());
//	ApplyAccelerationBoundaryConditionsForNodes(*domain);
//	CalcVelocityForNodes( *domain, delt, u_cut, domain->numNode()) ;
//	CalcPositionForNodes( *domain, delt, domain->numNode() );
//	
//	SYNC(elements);
//	EXIT_TP();
//
//}

void NodalCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	//LagrangeNodal(*domain);
	
	const Real_t delt = domain->deltatime() ;
	Real_t u_cut = domain->u_cut() ;
	
	/* time of boundary condition evaluation is beginning of step for force and
	* acceleration boundary conditions. */
	CalcForceForNodes(*domain);
	CalcAccelerationForNodes(*domain, domain->numNode());
	ApplyAccelerationBoundaryConditionsForNodes(*domain);
	CalcVelocityForNodes( *domain, delt, u_cut, domain->numNode()) ;
	CalcPositionForNodes( *domain, delt, domain->numNode() );
	
	SYNC(elements);
	EXIT_TP();

}


void ElementsCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	LagrangeElements(*domain,domain->numElem());

	SYNC(timeConstraints);
	EXIT_TP();
}

void TimeConstraintsCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	Domain *domain=FRAME(domain);
	CalcTimeConstraintsForElems(*domain);
	SYNC(sync);
	EXIT_TP();
}
void LeapFrogSyncCD::fire(void)
{
	LOAD_FRAME(LagrangeLeapFrogTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
