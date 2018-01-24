#include "CalcSoundSpeedForElemsTP.h"


void  CalcSoundSpeedForElemsIterCD ::fire(void)
{
	LOAD_FRAME( CalcSoundSpeedForElemsTP);

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
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
	
	CalcSoundSpeedForElems_darts(*domain,vnewc, rho0, e_new, p_new,pbvc, bvc,lw,hi, regElemList) ;
	
	SIGNAL(signalUp);

	EXIT_TP();
}
