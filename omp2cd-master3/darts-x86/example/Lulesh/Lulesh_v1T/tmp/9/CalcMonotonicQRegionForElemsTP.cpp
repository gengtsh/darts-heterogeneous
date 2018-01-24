#include "CalcMonotonicQRegionForElemsTP.h"


/* Sum contributions to total stress tensor */
void  CalcMonotonicQRegionForElemsIterCD ::fire(void)
{
	LOAD_FRAME( CalcMonotonicQRegionForElemsTP);
	
	Index_t numZone=FRAME(numZone);
	
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(numZone == FRAME(numZoneInit)){
//		if(IdL<N_CORES){
//			FRAME( calcMonotonicQRegionForElemsIter[IdL]) = CalcMonotonicQRegionForElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
//			ADD( calcMonotonicQRegionForElemsIter+IdL);
//		}
//		if(IdR<N_CORES){
//			FRAME( calcMonotonicQRegionForElemsIter[IdR]) = CalcMonotonicQRegionForElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
//			ADD( calcMonotonicQRegionForElemsIter+IdR);
//		}
//	}else{
//			if(IdL<N_CORES){
//				SYNC( calcMonotonicQRegionForElemsIter[IdL]);
//			}
//			if(IdR<N_CORES){
//				SYNC( calcMonotonicQRegionForElemsIter[IdR]);
//			}
//	}


	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	if(numZone == FRAME(numZoneInit)){
		for (size_t i=IdC0;i<tree;++i){
			FRAME( calcMonotonicQRegionForElemsIter[i]) = CalcMonotonicQRegionForElemsIterCD{0,1,getTP(),SHORTWAIT,i};
			ADD( calcMonotonicQRegionForElemsIter+i);
		}
	}else{
		for (size_t i=IdC0;i<tree;++i){
			SYNC( calcMonotonicQRegionForElemsIter[i]);
		}
	}

	Domain *domain	=FRAME(domain);	
	Real_t *vnew	=FRAME(vnew);
	Real_t	ptiny	=FRAME(ptiny);
	Index_t numElemSize	= domain->regElemSize(numZone);
	
	size_t	Chunk = numElemSize/N_CORES;
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

	RESET(calcMonotonicQRegionForElemsSync);
	
	Index_t idx;
	Index_t numReg= domain->numReg();
	
	if (numZone==(numReg-1)){
		SIGNAL(signalUp);
	}else{
		for(idx=numZone+1;idx<numReg;++idx){
			if(domain->regElemSize(idx)>0){
////				for(size_t i=0;i<N_CORES;++i){
////					SYNC( calcMonotonicQRegionForElemsIter[i] );
////				}
				
				//size_t tree = MIN(g_treeBarrier,N_CORES);
				for ( size_t i = 0; i < N_TREE; ++i ) {
					SYNC(calcMonotonicQRegionForElemsIter[i] );
				}
				
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
