#include "CalcAccelerationForNodesTP.h"


void  CalcAccelerationForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcAccelerationForNodesTP);
	size_t	Id	= getID();
	std::cout<<"CalcAccelerationForNodesIter["<<Id<<"] is running!"<<std::endl;

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcAccelerationForNodesIter[i]) = CalcAccelerationForNodesIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD(calcAccelerationForNodesIter+i );
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
