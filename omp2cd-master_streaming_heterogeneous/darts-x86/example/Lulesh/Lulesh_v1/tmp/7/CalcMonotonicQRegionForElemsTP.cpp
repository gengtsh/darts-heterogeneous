#include "CalcMonotonicQRegionForElemsTP.h"


/* Sum contributions to total stress tensor */
void  CalcMonotonicQRegionForElemsIterCD ::fire(void)
{
	LOAD_FRAME( CalcMonotonicQRegionForElemsTP);
	Domain *domain	=FRAME(domain);	
	Real_t *vnew	=FRAME(vnew);
	Real_t	ptiny	=FRAME(ptiny);
	Real_t	numZone	=FRAME(numZone);
	Index_t numElemSize	= domain->regElemSize(numZone);
	
	size_t	Chunk = numElemSize/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemSize%N_CORES):0) + (Id+1)* Chunk;
	CalcMonotonicQRegionForElems_darts(*domain,numZone,vnew,ptiny,lw,hi);

	RESET(calcMonotonicQRegionForElemsIter[Id]);

	SYNC(calcMonotonicQRegionForElemsSync);	
	EXIT_TP();
}

void CalcMonotonicQRegionForElemsSyncCD::fire(void)
{
	LOAD_FRAME( CalcMonotonicQRegionForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t numZone=FRAME(numZone);

	Index_t idx;
	Index_t numReg= domain->numReg();
	
	if (numZone==(numReg-1)){
		SIGNAL(signalUp);
	}else{
		for(idx=numZone+1;idx<numReg;++idx){
			if(domain->regElemSize(idx)){
				for(size_t i=0;i<N_CORES;++i){
					SYNC( calcMonotonicQRegionForElemsIter[i] );
				}
				RESET(calcMonotonicQRegionForElemsSync);
				FRAME(numZone) = idx;
				break;
			}
		}
		if((FRAME(numZone)!=(numReg-1))&&(idx==numReg)){
			SIGNAL(signalUp);
		}
	}
	
	EXIT_TP();
}
