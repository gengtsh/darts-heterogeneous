#include "CalcPressureForElemsTP.h"


void  CalcPressureForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcPressureForElemsTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcPressureForElemsP1Iter[IdL])=CalcPressureForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdL};
//		ADD( calcPressureForElemsP1Iter+ IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(calcPressureForElemsP1Iter[IdR])=CalcPressureForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdR};
//		ADD( calcPressureForElemsP1Iter+ IdR);
//	}

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);

	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcPressureForElemsP1Iter[i])=CalcPressureForElemsP1IterCD{0,1,getTP(),SHORTWAIT,i};
		ADD( calcPressureForElemsP1Iter+ i);
	}

	Real_t  *bvc		 = FRAME(bvc		);
	Real_t  *pbvc		 = FRAME(pbvc		);
	Real_t  *compression = FRAME(compression);
	Index_t  numElemReg	 = FRAME(numElemReg	);
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
	
	CalcPressureForElemsP1_darts(bvc,pbvc,compression,lw,hi);
	
////	SYNC(calcPressureForElemsP2Iter[Id]);

	FRAME(calcPressureForElemsP2Iter[Id])=CalcPressureForElemsP2IterCD{0,1,getTP(),SHORTWAIT,Id};
	ADD( calcPressureForElemsP2Iter+ Id);
	//std::cout<<"CalcPressureForElemsP1Iter["<<Id<<"],parent's TP address="<<getTP()->parentTP_<<",ref_="<<getTP()->parentTP_->getRef()<<std::endl;
	EXIT_TP();
}

void  CalcPressureForElemsP2IterCD ::fire(void)
{
	LOAD_FRAME( CalcPressureForElemsTP);
	
	Real_t  *p_new		 = FRAME(p_new		);
	Real_t  *bvc		 = FRAME(bvc		);
	Real_t  *e_old		 = FRAME(e_old		);
	Real_t  *vnewc		 = FRAME(vnewc		);
	Real_t   pmin		 = FRAME(pmin		);
	Real_t   p_cut		 = FRAME(p_cut		);
	Real_t   eosvmax	 = FRAME(eosvmax	);
	Index_t  numElemReg	 = FRAME(numElemReg	);
	Index_t *regElemList = FRAME(regElemList); 
	
	size_t	Chunk = numElemReg/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
	
	CalcPressureForElemsP2_darts(p_new,bvc,e_old,vnewc,pmin,p_cut,eosvmax,lw,hi, regElemList);
	
	SYNC(calcPressureForElemsSync);

//	std::cout<<"CalcPressureForElemsP2Iter["<<Id<<"],parent's TP address="<<getTP()->parentTP_<<",ref_="<<getTP()->parentTP_->getRef()<<std::endl;
	EXIT_TP();
}


void  CalcPressureForElemsSyncCD ::fire(void)
{
	LOAD_FRAME( CalcPressureForElemsTP);

//	std::cout<<"CalcPressureForeElems Parents'TP address="<<getTP()->parentTP_<<", CalcPressureForeElemsTP.signalUp,address="<<FRAME(signalUp)<<",counter="<<FRAME(signalUp)->getCounter()<<std::endl;	 
//	std::cout<<"CalcPressureForeElemsSync.signalUp's TP.ref_="<<FRAME(signalUp)->getTP()->getRef()<<std::endl;	
	SIGNAL(signalUp);
//	std::cout<<"0"<<std::endl;
	EXIT_TP();
}
