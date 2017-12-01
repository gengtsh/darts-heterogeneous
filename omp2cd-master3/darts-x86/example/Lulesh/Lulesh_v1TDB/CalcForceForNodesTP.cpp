#include "CalcForceForNodesTP.h"
#include "CalcVolumeForceForElemsTP.h"

void CalcForceForNodesP1IterCD::fire(void)
{
	LOAD_FRAME(CalcForceForNodesTP);

	size_t	Id	= getID();
	std::cout<<"CalcForceForNodelsP1Iter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcForceForNodesP1Iter[i])=CalcForceForNodesP1IterCD{0,1,getTP(),SHORTWAIT,i};
		ADD (calcForceForNodesP1Iter + i);
	}
	
	Domain *domain=FRAME(domain);
	Index_t numNode = domain->numNode() ;
	size_t	Chunk = numNode/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;

	CalcForceForNodesP1(*domain,lw,hi);

	SYNC(calcVolumeForceForElems);
	EXIT_TP();
}

/* Calcforce calls partial, force, hourq */
void CalcVolumeForceForElemsCD::fire(void)
{
	std::cout<<"CalcVolumeForceForElemsTP begin!"<<std::endl;	
	LOAD_FRAME(CalcForceForNodesTP);
	Domain *domain =FRAME(domain);

	/* Calcforce calls partial, force, hourq */
	INVOKE(CalcVolumeForceForElemsTP, domain,FRAME(signalUp)) ;

	EXIT_TP();

}
