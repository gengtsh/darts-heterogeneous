#include "CalcHydroConstraintForElemsTP.h"

void CalcHydroConstraintForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcHydroConstraintForElemsTP);
	
	size_t	Id	= getID();
	Index_t numZone = FRAME(numZone);
	std::cout<<"CalcHydroConstraintForElemsP1Iter["<<Id<<"] in Zone:"<<numZone<<" , is running!"<<std::endl;

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	if(numZone == FRAME(numZoneInit)){
		for (size_t i=IdC0;i<tree;++i){
			FRAME(calcHydroConstraintForElemsP1Iter[i])= CalcHydroConstraintForElemsP1IterCD{0,1,getTP(),SHORTWAIT,i};
			ADD ( calcHydroConstraintForElemsP1Iter+ i);
		}
	}else{
		for (size_t i=IdC0;i<tree;++i){
			SYNC( calcHydroConstraintForElemsP1Iter[i]);
		}
	}

	Domain *domain =FRAME(domain);

	std::cout <<"CalcHydroConstraintForElemsTP IterCD["<<Id<<"].counter : "<<FRAME(calcHydroConstraintForElemsP1Iter[Id]).getCounter()<<",numZone :" <<numZone<<std::endl;
	//Index_t regElemSize	=FRAME(regElemSize);
	//Index_t *regElemlist=FRAME(regElemlist);
	
	Index_t	regElemSize = domain->regElemSize(numZone);
	Index_t *regElemlist = domain->regElemlist(numZone);
	Index_t *hydro_elem_per_thread=FRAME(hydro_elem_per_thread);
	Real_t *dthydro_per_thread	= FRAME(dthydro_per_thread);

	size_t	Chunk = regElemSize/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (regElemSize%N_CORES):0) + (Id+1)* Chunk;

	RESET(calcHydroConstraintForElemsP1Iter[Id] );
	CalcHydroConstraintForElemsP1_darts(*domain,regElemlist,domain->dvovmax(),domain-> dthydro(),hydro_elem_per_thread, dthydro_per_thread,Id,lw,hi );

	SYNC(calcHydroConstraintForElemsSync );

	
	
	EXIT_TP();
}

void  CalcHydroConstraintForElemsSyncCD::fire(void)
{

	std::cout <<"CalcHydroConstraintForElemsSync is running!"<<std::endl;
	LOAD_FRAME( CalcHydroConstraintForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t *hydro_elem_per_thread=FRAME(hydro_elem_per_thread);
	Real_t *dthydro_per_thread= FRAME(dthydro_per_thread );
	
	RESET(calcHydroConstraintForElemsSync );
//	std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter: "<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	CalcHydroConstraintForElemsP2_darts(domain->dthydro(),hydro_elem_per_thread,dthydro_per_thread,N_CORES);

//	std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter again: "<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	Index_t numZone = FRAME(numZone);
	Index_t numRegMinus1  = FRAME(numRegMinus1);
//	std::cout<<"CalcHydroConstraintForElemsTP, numZone:"<<numZone<<",numRegMinus1:"<<numRegMinus1<<std::endl;
	if (numZone < numRegMinus1){
// 		++FRAME(numZone);
	//	FRAME(regElemSize) = domain->regElemSize(numZone);
	//	FRAME(regElemlist) = domain->regElemlist(numZone);
////		for (size_t i=0; i<N_CORES;++i){
////			SYNC(calcHydroConstraintForElemsP1Iter[i]);
////		}
	
//		std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter again2:"<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<",numZone:"<<numZone<<std::endl;

		//size_t tree = MIN(g_treeBarrier,N_CORES);
		for ( size_t i = 0; i < N_TREE; ++i ) {
			SYNC(calcHydroConstraintForElemsP1Iter[i]);
		}
		
//		std::cout <<"CalcHydroConstraintForElemsTP IterCD.counter again3:"<<FRAME(calcHydroConstraintForElemsP1Iter[0]).getCounter()<<",numZone:"<<numZone<<std::endl;
	
	}else{
		SIGNAL(signalUp);

		std::cout <<"CalcHydroConstraintForElemsTP finished! "<<std::endl;
	}

	++FRAME(numZone);
	
	//for (size_t i=0;i<N_TREE;++i){
	//	std::cout <<"CalcHydroConstraintForElemsTP IterCD["<<i<<"].counter final: "<<FRAME(calcHydroConstraintForElemsP1Iter[i]).getCounter()<<",numZone :" <<numZone<<std::endl;
	//}
	EXIT_TP();
}


