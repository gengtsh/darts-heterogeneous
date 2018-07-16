#include "UpdateVolumesForElemsTP.h"


/* Sum contributions to total stress tensor */
void  UpdateVolumesForElemsIterCD ::fire(void)
{
	LOAD_FRAME( UpdateVolumesForElemsTP);
	size_t	Id	= getID();
	std::cout<<"UpdateVolumesForElemsIter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(upadteVolumesForElemsIter[i])= UpdateVolumesForElemsIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD( upadteVolumesForElemsIter+ i);

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
