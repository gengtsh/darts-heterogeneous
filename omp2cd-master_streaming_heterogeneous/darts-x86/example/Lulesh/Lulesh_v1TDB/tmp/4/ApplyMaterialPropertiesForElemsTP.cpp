#include "ApplyMaterialPropertiesForElemsTP.h"

void ApplyMaterialPropertiesForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *vnew  =FRAME(vnew);
	Index_t numElem=domain->numElem(); 		

	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
    
	ApplyMaterialPropertiesForElemsP1_darts(*domain,vnew,lw,hi) ;
	
	SYNC(applyMaterialPropertiesForElemsP2);

	EXIT_TP();
}



void  ApplyMaterialPropertiesForElemsP2CD ::fire(void)
{
	
	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);

	Domain *domain	=FRAME(domain);
	Real_t *vnew	=FRAME(vnew);
	
    ApplyMaterialPropertiesForElemsP2_darts(*domain,vnew) ;
	SIGNAL(signalUp);
	EXIT_TP();
}

	
