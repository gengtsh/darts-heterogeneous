#include "UpdateVolumesForElemsTP.h"


/* Sum contributions to total stress tensor */
void  UpdateVolumesForElemsIterCD ::fire(void)
{
	LOAD_FRAME( UpdateVolumesForElemsTP);
	Domain *domain =FRAME(domain);	
	Real_t *vnew   =FRAME(vnew);
	Real_t v_cut	=FRAME(v_cut);	
	Index_t numElem=FRAME(numElem);
	
	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
	
	UpdateVolumesForElems_darts(*domain, vnew,v_cut,lw,hi) ;

	SIGNAL(signalUp);

	EXIT_TP();
}
