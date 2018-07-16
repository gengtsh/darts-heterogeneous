#include "IntegrateStressForElemsP2TP.h"

void IntegrateStressForElemsP2IterCD ::fire(void)
{
	LOAD_FRAME(IntegrateStressForElemsP2TP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(integrateStressForElemsP2Iter[IdL])=IntegrateStressForElemsP2IterCD{0,1,getTP(),SHORTWAIT,IdL};	
		ADD ( integrateStressForElemsP2Iter +IdL);
	}
	if(IdR<N_CORES){
		FRAME(integrateStressForElemsP2Iter[IdR])=IntegrateStressForElemsP2IterCD{0,1,getTP(),SHORTWAIT,IdR};	
		ADD ( integrateStressForElemsP2Iter +IdR);
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


	
