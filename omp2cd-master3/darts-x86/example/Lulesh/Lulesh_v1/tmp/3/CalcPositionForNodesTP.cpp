#include "CalcPositionForNodesTP.h"


void  CalcPositionForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcPositionForNodesTP);
	Domain *domain =FRAME(domain);	
	
	Index_t numNode = domain->numNode() ;
	
	size_t	Chunk = numNode/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
	
    const Real_t delt = domain->deltatime() ;
	
	CalcPositionForNodes_darts(*domain, delt, lw,hi);

	SIGNAL(signalUp);

	EXIT_TP();
}
