#include "IntegrateStressForElemsP2TP.h"

void IntegrateStressForElemsP2IterCD ::fire(void)
{
	LOAD_FRAME(IntegrateStressForElemsP2TP);
	size_t	Id	= getID();
	std::cout<<"IntegrateStressForElemsP2Iter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(integrateStressForElemsP2Iter[i])=IntegrateStressForElemsP2IterCD{0,1,getTP(),SHORTWAIT,i};	
		ADD ( integrateStressForElemsP2Iter +i);
	}

	Domain *domain=FRAME(domain);
	Real_t *fx_elem=FRAME(fx_elem);
	Real_t *fy_elem=FRAME(fy_elem);
	Real_t *fz_elem=FRAME(fz_elem);
	Index_t numNode=FRAME(numNode);
	
	size_t	Chunk = numNode/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
    
    IntegrateStressForElemsP2(*domain,fx_elem,fy_elem,fz_elem,lw,hi) ;

	SIGNAL(signalUp);
	EXIT_TP();
}


	
