#include "IntegrateStressForElemsTP.h"
#include "IntegrateStressForElemsP2TP.h"

// call elemlib stress integration loop to produce nodal forces from material stresses.
void IntegrateStressForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(IntegrateStressForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *sigxx =FRAME(sigxx); 
	Real_t *sigyy =FRAME(sigyy);
	Real_t *sigzz =FRAME(sigzz);
	Real_t *determ =FRAME(determ);
	Index_t numElem=FRAME(numElem); 		
	Real_t *fx_elem=FRAME(fx_elem);
	Real_t *fy_elem=FRAME(fy_elem);
	Real_t *fz_elem=FRAME(fz_elem);

	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
    
	IntegrateStressForElemsP1(*domain,sigxx, sigyy, sigzz, determ, fx_elem,fy_elem,fz_elem,lw,hi) ;
	
	SYNC(integrateStressForElemsP2);


	EXIT_TP();
}

void  IntegrateStressForElemsP2CD ::fire(void)
{
	
	LOAD_FRAME(IntegrateStressForElemsTP);

	Domain *domain=FRAME(domain);
	Index_t numNode=FRAME(numNode);
	Real_t *fx_elem=FRAME(fx_elem);
	Real_t *fy_elem=FRAME(fy_elem);
	Real_t *fz_elem=FRAME(fz_elem);
	
//    IntegrateStressForElemsP2(*domain,numNode,fx_elem,fy_elem,fz_elem) ;
//	SIGNAL(signalUp);
	INVOKE(IntegrateStressForElemsP2TP, domain,numNode,fx_elem,fy_elem,fz_elem,FRAME(signalUp)) ;
	EXIT_TP();
}

	
