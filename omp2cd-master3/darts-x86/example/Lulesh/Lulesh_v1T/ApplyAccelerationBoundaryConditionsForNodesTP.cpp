#include "ApplyAccelerationBoundaryConditionsForNodesTP.h"


void  ApplyAccelerationBoundaryConditionsForNodesIterCD ::fire(void)
{
	LOAD_FRAME( ApplyAccelerationBoundaryConditionsForNodesTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(applyAccelerationBoundaryConditionsForNodesIter[IdL])= ApplyAccelerationBoundaryConditionsForNodesIterCD{0,1,getTP(),SHORTWAIT,IdL};
//		ADD ( applyAccelerationBoundaryConditionsForNodesIter+ IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(applyAccelerationBoundaryConditionsForNodesIter[IdR])= ApplyAccelerationBoundaryConditionsForNodesIterCD{0,1,getTP(),SHORTWAIT,IdR};
//		ADD ( applyAccelerationBoundaryConditionsForNodesIter+ IdR);
//	}
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(applyAccelerationBoundaryConditionsForNodesIter[i])= ApplyAccelerationBoundaryConditionsForNodesIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( applyAccelerationBoundaryConditionsForNodesIter+ i);
	}

	Domain *domain =FRAME(domain);	
	
	Index_t size = domain->sizeX();
	Index_t numNodeBC = (size+1)*(size+1) ;
	
	size_t	Chunk = numNodeBC/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNodeBC%N_CORES):0) + (Id+1)* Chunk;
	ApplyAccelerationBoundaryConditionsForNodes_darts(*domain,lw,hi);
	SIGNAL(signalUp);

	EXIT_TP();
}
