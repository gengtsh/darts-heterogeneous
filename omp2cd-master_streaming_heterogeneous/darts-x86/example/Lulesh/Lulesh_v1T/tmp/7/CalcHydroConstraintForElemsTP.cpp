#include "CalcHydroConstraintForElemsTP.h"

void CalcHydroConstraintForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcHydroConstraintForElemsTP);
	Domain *domain =FRAME(domain);

	//Index_t regElemSize	=FRAME(regElemSize);
	//Index_t *regElemlist=FRAME(regElemlist);
	
	Index_t numZone = FRAME(numZone);
	Index_t	regElemSize = domain->regElemSize(numZone);
	Index_t *regElemlist = domain->regElemlist(numZone);
	Index_t *hydro_elem_per_thread=FRAME(hydro_elem_per_thread);
	Real_t *dthydro_per_thread	= FRAME(dthydro_per_thread);

	size_t	Chunk = regElemSize/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (regElemSize%N_CORES):0) + (Id+1)* Chunk;
	CalcHydroConstraintForElemsP1_darts(*domain,regElemlist,domain->dvovmax(),domain-> dthydro(),hydro_elem_per_thread, dthydro_per_thread,Id,lw,hi );

	SYNC(calcHydroConstraintForElemsSync );

	Index_t numRegMinus1  = FRAME(numRegMinus1);
	if (numZone < numRegMinus1){
		RESET(calcHydroConstraintForElemsP1Iter[Id] );
	}
	
	EXIT_TP();
}

void  CalcHydroConstraintForElemsSyncCD::fire(void)
{

	LOAD_FRAME( CalcHydroConstraintForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t *hydro_elem_per_thread=FRAME(hydro_elem_per_thread);
	Real_t *dthydro_per_thread= FRAME(dthydro_per_thread );
	
	std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter: "<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	CalcHydroConstraintForElemsP2_darts(domain->dthydro(),hydro_elem_per_thread,dthydro_per_thread,N_CORES);

	std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter again: "<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	Index_t numZone = FRAME(numZone);
	Index_t numRegMinus1  = FRAME(numRegMinus1);
	if (numZone < numRegMinus1){
		numZone= ++FRAME(numZone);
	//	FRAME(regElemSize) = domain->regElemSize(numZone);
	//	FRAME(regElemlist) = domain->regElemlist(numZone);
		RESET(calcHydroConstraintForElemsSync );
		for (size_t i=0; i<N_CORES;++i){
			SYNC(calcHydroConstraintForElemsP1Iter[i]);
		}
	}else{
		SIGNAL(signalUp);
	}

	std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter: "<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<",numZone :" <<numZone<<std::endl;
	EXIT_TP();
}


