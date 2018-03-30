#include "CalcAccelerationForNodesTP.h"


void  CalcAccelerationForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcAccelerationForNodesTP);
	size_t	Id	= getID();

	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcAccelerationForNodesIter[IdL]) = CalcAccelerationForNodesIterCD{0,1,getTP(),SHORTWAIT,IdL};
		ADD(calcAccelerationForNodesIter+IdL );
	}
	if(IdR<N_CORES){
		FRAME(calcAccelerationForNodesIter[IdR]) = CalcAccelerationForNodesIterCD{0,1,getTP(),SHORTWAIT,IdR};
		ADD(calcAccelerationForNodesIter+IdR );
	}

	Domain *domain =FRAME(domain);	
	Index_t numNode=FRAME(numNode);
	
	size_t	Chunk = numNode/N_CORES;
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
