#include "CalcMonotonicQRegionForElemsTP.h"


/* Sum contributions to total stress tensor */
void  CalcMonotonicQRegionForElemsIterCD ::fire(void)
{
	LOAD_FRAME( CalcMonotonicQRegionForElemsTP);
	Domain *domain	=FRAME(domain);	
	Real_t *vnew	=FRAME(vnew);
	Real_t	ptiny	=FRAME(ptiny);
	Real_t	numReg	=FRAME(numReg);
	Index_t numElemSize	= domain->regElemSize(numReg);
	
	size_t	Chunk = numElemSize/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemSize%N_CORES):0) + (Id+1)* Chunk;
	CalcMonotonicQRegionForElems_darts(*domain,numReg,vnew,ptiny,lw,hi);
	SYNC(calcMonotonicQRegionForElemsSync);	

	EXIT_TP();
}

void CalcMonotonicQRegionForElemsSyncCD::fire(void)
{
	LOAD_FRAME( CalcMonotonicQRegionForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
