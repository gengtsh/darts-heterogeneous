#include "CalcPositionForNodesTP.h"


void  CalcPositionForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcPositionForNodesTP);
	size_t	Id	= getID();
	std::cout<<" CalcPositionForNodesIter["<<Id<<"] is running!"<<std::endl;
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
