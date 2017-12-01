#include "CalcLagrangeElementsP2TP.h"


/* Sum contributions to total stress tensor */
void  CalcLagrangeElementsP2IterCD ::fire(void)
{
	LOAD_FRAME( CalcLagrangeElementsP2TP);
	Domain *domain =FRAME(domain);	
	Real_t *vnew =FRAME(vnew);
	Index_t numElem=domain->numElem();
	
	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
	CalcLagrangeElementsP2_darts(*domain,vnew,lw,hi);
	
	SIGNAL(signalUp);

	EXIT_TP();
}
