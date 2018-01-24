#include "CalcCourantConstraintForElemsTP.h"


void  CalcCourantConstraintForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	Index_t numZone = FRAME(numZone);
	
	if(numZone == FRAME(numZoneInit)){
		if(IdL<N_CORES){
			FRAME(calcCourantConstraintForElemsP1Iter[IdL])= CalcCourantConstraintForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdL};
			ADD( calcCourantConstraintForElemsP1Iter+ IdL);
		}

		if(IdR<N_CORES){
			FRAME(calcCourantConstraintForElemsP1Iter[IdR])= CalcCourantConstraintForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdR};
			ADD( calcCourantConstraintForElemsP1Iter+ IdR);
		}
	}else{
		if(IdL<N_CORES){
			SYNC( calcCourantConstraintForElemsP1Iter[IdL]);
		}
		if(IdR<N_CORES){
			SYNC( calcCourantConstraintForElemsP1Iter[IdR]);
		}
	}
	
	
	Domain *domain =FRAME(domain);

	//Index_t regElemSize	=FRAME(regElemSize);
	//Index_t *regElemlist=FRAME(regElemlist);
	
	Index_t regElemSize = domain->regElemSize(numZone);
	Index_t *regElemlist = domain->regElemlist(numZone);
	Index_t *courant_elem_per_thread	=FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread	= FRAME(dtcourant_per_thread );

	size_t	Chunk = regElemSize/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (regElemSize%N_CORES):0) + (Id+1)* Chunk;
	CalcCourantConstraintForElemsP1_darts(*domain,regElemlist,domain->qqc(),domain->dtcourant(),courant_elem_per_thread, dtcourant_per_thread,Id,lw,hi );

	SYNC(calcCourantConstraintForElemsSync );
	
	RESET(calcCourantConstraintForElemsP1Iter[Id] );

	EXIT_TP();
}

void CalcCourantConstraintForElemsSyncCD::fire(void)
{

	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t *courant_elem_per_thread=FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread= FRAME(dtcourant_per_thread );

	RESET( calcCourantConstraintForElemsSync );
	
	//std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;

	CalcCourantConstraintForElemsP2_darts(domain->dtcourant(),courant_elem_per_thread,dtcourant_per_thread,N_CORES);

	//std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter again:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	Index_t numZone = FRAME(numZone);
	Index_t numRegMinus1  = FRAME(numRegMinus1);
	if (numZone < numRegMinus1){
		++FRAME(numZone);
	//	FRAME(regElemSize) = domain->regElemSize(numZone);
	//	FRAME(regElemlist) = domain->regElemlist(numZone);
////		for(size_t i=0;i<N_CORES;++i){
////			SYNC(calcCourantConstraintForElemsP1Iter[i] );
////		}

		SYNC(calcCourantConstraintForElemsP1Iter[0] );
		if(N_CORES>1){
			SYNC(calcCourantConstraintForElemsP1Iter[1] );
		}
	
	}else{
		SIGNAL(signalUp);

	}
	//std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<< ",numZone :" <<numZone<<std::endl;
	EXIT_TP();
}


