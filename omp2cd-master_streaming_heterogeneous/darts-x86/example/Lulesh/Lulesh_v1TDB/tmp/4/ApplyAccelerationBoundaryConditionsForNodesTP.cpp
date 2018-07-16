#include "ApplyAccelerationBoundaryConditionsForNodesTP.h"


void  ApplyAccelerationBoundaryConditionsForNodesIterCD ::fire(void)
{
	LOAD_FRAME( ApplyAccelerationBoundaryConditionsForNodesTP);
	Domain *domain =FRAME(domain);	
	
	Index_t size = domain->sizeX();
	Index_t numNodeBC = (size+1)*(size+1) ;
	
	size_t	Chunk = numNodeBC/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNodeBC%N_CORES):0) + (Id+1)* Chunk;
	

//    if (!domain->symmXempty() != 0) {
//       for(Index_t i=lw ; i<hi ; ++i)
//          domain->xdd(domain->symmX(i)) = Real_t(0.0) ;
//    }
//
//    if (!domain->symmYempty() != 0) {
//       for(Index_t i=lw ; i<hi ; ++i)
//          domain->ydd(domain->symmY(i)) = Real_t(0.0) ;
//    }
//
//    if (!domain->symmZempty() != 0) {
//       for(Index_t i=lw ; i<hi ; ++i)
//          domain->zdd(domain->symmZ(i)) = Real_t(0.0) ;
//    }
	ApplyAccelerationBoundaryConditionsForNodes_darts(*domain,lw,hi);

	SIGNAL(signalUp);

	EXIT_TP();
}
