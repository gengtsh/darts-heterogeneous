#include "CalcVelocityForNodesTP.h"


void  CalcVelocityForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcVelocityForNodesTP);
	Domain *domain =FRAME(domain);	
	Index_t numNode=domain->numNode();	
	
	const Real_t delt = domain->deltatime() ;
	Real_t u_cut = domain->u_cut() ;
	
	size_t	Chunk = numNode/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;

	CalcVelocityForNodes_darts( *domain, delt, u_cut, lw,hi) ;

	SIGNAL(signalUp);

	EXIT_TP();
}
