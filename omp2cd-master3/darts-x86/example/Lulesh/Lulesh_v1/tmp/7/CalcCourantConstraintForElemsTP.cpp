#include "CalcCourantConstraintForElemsTP.h"


void  CalcCourantConstraintForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	Domain *domain =FRAME(domain);

	//Index_t regElemSize	=FRAME(regElemSize);
	//Index_t *regElemlist=FRAME(regElemlist);

	Index_t numZone = FRAME(numZone);
	Index_t regElemSize = domain->regElemSize(numZone);
	Index_t *regElemlist = domain->regElemlist(numZone);
	Index_t *courant_elem_per_thread	=FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread	= FRAME(dtcourant_per_thread );

	size_t	Chunk = regElemSize/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (regElemSize%N_CORES):0) + (Id+1)* Chunk;
	CalcCourantConstraintForElemsP1_darts(*domain,regElemlist,domain->qqc(),domain->dtcourant(),courant_elem_per_thread, dtcourant_per_thread,Id,lw,hi );

	SYNC(calcCourantConstraintForElemsSync );
	
	Index_t numRegMinus1  = FRAME(numRegMinus1);
	if (numZone < numRegMinus1){
		RESET(calcCourantConstraintForElemsP1Iter[Id] );
	}
	EXIT_TP();
}

void CalcCourantConstraintForElemsSyncCD::fire(void)
{

	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t *courant_elem_per_thread=FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread= FRAME(dtcourant_per_thread );
	
	std::cout<<"calcCourantConstraintForElemsTP Iter.counter init:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	CalcCourantConstraintForElemsP2_darts(domain->dtcourant(),courant_elem_per_thread,dtcourant_per_thread,N_CORES);
	std::cout<<"calcCourantConstraintForElemsTP Iter.counter midd:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	Index_t numZone = FRAME(numZone);
	Index_t numRegMinus1  = FRAME(numRegMinus1);
	if (numZone < numRegMinus1){
		numZone= ++FRAME(numZone);
//		FRAME(regElemSize) = domain->regElemSize(numZone);
//		FRAME(regElemlist) = domain->regElemlist(numZone);
		RESET( calcCourantConstraintForElemsSync );
		for(size_t i=0;i<N_CORES;++i){
			SYNC(calcCourantConstraintForElemsP1Iter[i] );
		}
	}else{
		SIGNAL(signalUp);

	}

	EXIT_TP();
}


