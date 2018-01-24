#include "CalcPositionForNodesTP.h"


void  CalcPositionForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcPositionForNodesTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcPositionForNodesIter[IdL])= CalcPositionForNodesIterCD{0,1,getTP(),SHORTWAIT,IdL};
//		ADD ( calcPositionForNodesIter+ IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(calcPositionForNodesIter[IdR])= CalcPositionForNodesIterCD{0,1,getTP(),SHORTWAIT,IdR};
//		ADD ( calcPositionForNodesIter+ IdR);
//	}

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcPositionForNodesIter[i])= CalcPositionForNodesIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( calcPositionForNodesIter+ i);
	}

	Domain *domain =FRAME(domain);	
	
	Index_t numNode = domain->numNode() ;
	
	size_t	Chunk = numNode/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
	
    const Real_t delt = domain->deltatime() ;
	
	CalcPositionForNodes_darts(*domain, delt, lw,hi);

	SIGNAL(signalUp);

	EXIT_TP();
}
