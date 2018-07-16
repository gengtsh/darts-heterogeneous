#include "CalcPressureForElemsTP.h"


void  CalcPressureForElemsP1IterCD ::fire(void)
{
	LOAD_FRAME( CalcPressureForElemsTP);
	
	Real_t  *bvc		 = FRAME(bvc		);
	Real_t  *pbvc		 = FRAME(pbvc		);
	Real_t  *compression = FRAME(compression);
	Index_t  numElemReg	 = FRAME(numElemReg	);
	
	size_t	Chunk = numElemReg/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
	
	CalcPressureForElemsP1_darts(bvc,pbvc,compression,lw,hi);
	
	SYNC(calcPressureForElemsP2Iter[Id]);

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
	
	//SIGNAL(signalUp);
	SYNC(calcPressureForElemsSync);

	EXIT_TP();
}


void  CalcPressureForElemsSyncCD ::fire(void)
{

	LOAD_FRAME( CalcPressureForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
