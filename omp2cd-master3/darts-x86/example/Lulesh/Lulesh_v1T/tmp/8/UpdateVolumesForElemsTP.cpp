#include "UpdateVolumesForElemsTP.h"


/* Sum contributions to total stress tensor */
void  UpdateVolumesForElemsIterCD ::fire(void)
{
	LOAD_FRAME( UpdateVolumesForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(upadteVolumesForElemsIter[IdL])= UpdateVolumesForElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
		ADD( upadteVolumesForElemsIter+ IdL);
	}
	if(IdR<N_CORES){
		FRAME(upadteVolumesForElemsIter[IdR])= UpdateVolumesForElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
		ADD( upadteVolumesForElemsIter+ IdR);
	}
	Domain *domain =FRAME(domain);	
	Real_t *vnew   =FRAME(vnew);
	Real_t v_cut	=FRAME(v_cut);	
	Index_t numElem=FRAME(numElem);
	
	size_t	Chunk = numElem/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
	
	UpdateVolumesForElems_darts(*domain, vnew,v_cut,lw,hi) ;

	SIGNAL(signalUp);

	EXIT_TP();
}
