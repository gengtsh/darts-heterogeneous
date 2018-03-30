#include "CalcSoundSpeedForElemsTP.h"


void  CalcSoundSpeedForElemsIterCD ::fire(void)
{
	LOAD_FRAME( CalcSoundSpeedForElemsTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcSoundSpeedForElemsIter[IdL])=CalcSoundSpeedForElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
//		ADD( calcSoundSpeedForElemsIter+ IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(calcSoundSpeedForElemsIter[IdR])=CalcSoundSpeedForElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
//		ADD( calcSoundSpeedForElemsIter+ IdR);
//	}


	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcSoundSpeedForElemsIter[i])=CalcSoundSpeedForElemsIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD( calcSoundSpeedForElemsIter+ i);
	}


	Domain *domain	=FRAME(domain);
	Real_t *vnewc   = FRAME(vnewc);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 

	Real_t rho0			=FRAME( rho0		)	;
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *bvc			=FRAME( bvc			)	;
	Real_t *pbvc		=FRAME( pbvc		)	;


	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
	
	CalcSoundSpeedForElems_darts(*domain,vnewc, rho0, e_new, p_new,pbvc, bvc,lw,hi, regElemList) ;

	SIGNAL(signalUp);
	EXIT_TP();
}
