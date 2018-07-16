#include "IntegrateStressForElemsP2TP.h"

void IntegrateStressForElemsP2IterCD ::fire(void)
{
	LOAD_FRAME(IntegrateStressForElemsP2TP);
	
	Domain *domain=FRAME(domain);
	Real_t *fx_elem=FRAME(fx_elem);
	Real_t *fy_elem=FRAME(fy_elem);
	Real_t *fz_elem=FRAME(fz_elem);
	Index_t numNode=FRAME(numNode);
	
	size_t	Chunk = numNode/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
    
    IntegrateStressForElemsP2(*domain,fx_elem,fy_elem,fz_elem,lw,hi) ;

	SIGNAL(signalUp);
	EXIT_TP();
}


	
