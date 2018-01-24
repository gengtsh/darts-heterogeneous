#include "CalcCourantConstraintForElemsTP.h"


void  CalcCourantConstraintForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	
	size_t	Id	= getID();
	Index_t numZone = FRAME(numZone);
	std::cout<<"CalcCourantConstraintForElemsP1Iter["<<Id<<"] in Zone:"<< numZone<<" , is running!"<<std::endl;
	
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	if(numZone == FRAME(numZoneInit)){
		for (size_t i=IdC0;i<tree;++i){
			FRAME(calcCourantConstraintForElemsP1Iter[i])= CalcCourantConstraintForElemsP1IterCD{0,1,getTP(),SHORTWAIT,i};
			ADD( calcCourantConstraintForElemsP1Iter+ i);
		}
	}else{
		for (size_t i=IdC0;i<tree;++i){
			SYNC( calcCourantConstraintForElemsP1Iter[i]);
		}

	}

	Domain *domain =FRAME(domain);

	std::cout <<"CalcCourantConstraintForElemsTP IterCD["<<Id<<"].counter :"<<FRAME(calcCourantConstraintForElemsP1Iter[Id]).getCounter()<< ",numZone :" <<numZone<<std::endl;
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

	RESET(calcCourantConstraintForElemsP1Iter[Id] );
	
	CalcCourantConstraintForElemsP1_darts(*domain,regElemlist,domain->qqc(),domain->dtcourant(),courant_elem_per_thread, dtcourant_per_thread,Id,lw,hi );

	SYNC(calcCourantConstraintForElemsSync );
	

	EXIT_TP();
}

void CalcCourantConstraintForElemsSyncCD::fire(void)
{

	std::cout <<"CalcCourantConstraintForElemsSync is running!" <<std::endl;
	LOAD_FRAME( CalcCourantConstraintForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t *courant_elem_per_thread=FRAME(courant_elem_per_thread);
	Real_t *dtcourant_per_thread= FRAME(dtcourant_per_thread );

	RESET( calcCourantConstraintForElemsSync );
//	std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;

	CalcCourantConstraintForElemsP2_darts(domain->dtcourant(),courant_elem_per_thread,dtcourant_per_thread,N_CORES);

//	std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter again:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<std::endl;
	Index_t numZone = FRAME(numZone);
	Index_t numRegMinus1  = FRAME(numRegMinus1);
	std::cout<<"CalcCourantConstraintForElemsTP, numZone:"<<numZone<<",numRegMinus1:"<<numRegMinus1<<std::endl;
	if (numZone < numRegMinus1){
//		++FRAME(numZone);
	//	FRAME(regElemSize) = domain->regElemSize(numZone);
	//	FRAME(regElemlist) = domain->regElemlist(numZone);
////		for(size_t i=0;i<N_CORES;++i){
////			SYNC(calcCourantConstraintForElemsP1Iter[i] );
////		}

//		std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter again2:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<",numZone:"<<numZone<<std::endl;
		//size_t tree = MIN(g_treeBarrier,N_CORES);
		for ( size_t i = 0; i < N_TREE; ++i ) {
			SYNC(calcCourantConstraintForElemsP1Iter[i] );
		}

//		std::cout <<"CalcCourantConstraintForElemsTP IterCD.counter again3:"<<FRAME(calcCourantConstraintForElemsP1Iter[0]).getCounter()<<",numZone"<<numZone<<std::endl;
	
	}else{
		SIGNAL(signalUp);

		std::cout <<"CalcCourantConstraintForElemsTP finished!"<<std::endl;
	}

	++FRAME(numZone);
//	for (size_t i=0;i<N_TREE;++i){
//		std::cout <<"CalcCourantConstraintForElemsTP IterCD["<<i<<"].counter final:"<<FRAME(calcCourantConstraintForElemsP1Iter[i]).getCounter()<< ",numZone :" <<numZone<<std::endl;
//	}
	EXIT_TP();
}


