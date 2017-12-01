#include "CalcAccelerationForNodesTP.h"


void  CalcAccelerationForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcAccelerationForNodesTP);
	Domain *domain =FRAME(domain);	
	Index_t numNode=FRAME(numNode);
	
	size_t	Chunk = numNode/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
	
	for ( Index_t i=lw ; i<hi ; ++i ) {
		domain->xdd(i) = domain->fx(i) / domain->nodalMass(i);
		domain->ydd(i) = domain->fy(i) / domain->nodalMass(i);
		domain->zdd(i) = domain->fz(i) / domain->nodalMass(i);
	}
	
	SIGNAL(signalUp);

	EXIT_TP();
}
