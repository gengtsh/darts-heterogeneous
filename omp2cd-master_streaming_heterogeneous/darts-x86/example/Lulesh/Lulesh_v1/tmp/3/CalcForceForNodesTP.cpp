#include "CalcForceForNodesTP.h"
#include "CalcVolumeForceForElemsTP.h"

void CalcForceForNodesP1IterCD::fire(void)
{
	LOAD_FRAME(CalcForceForNodesTP);
	Domain *domain=FRAME(domain);
	Index_t numNode = domain->numNode() ;
	size_t	Chunk = numNode/N_CORES;
	size_t	Id	= getID();
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
	LOAD_FRAME(CalcForceForNodesTP);
	Domain *domain =FRAME(domain);

	/* Calcforce calls partial, force, hourq */
	//CalcVolumeForceForElems(*domain) ;
	//SIGNAL(signalUp);
	INVOKE(CalcVolumeForceForElemsTP, domain,FRAME(signalUp)) ;

	EXIT_TP();

}
