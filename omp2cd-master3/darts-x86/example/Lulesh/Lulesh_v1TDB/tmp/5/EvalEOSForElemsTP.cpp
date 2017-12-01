#include "EvalEOSForElemsTP.h"
#include "CalcSoundSpeedForElemsTP.h"
#include "EvalEOSForElemsP1TP.h"

void EvalEOSForElemsP1CD ::fire(void)
{
	LOAD_FRAME(EvalEOSForElemsTP);
	
	Domain *domain = FRAME(domain);
	Real_t *vnewc   = FRAME(vnew);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 
	Int_t rep			= FRAME(rep);
	                                         
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *q_new		=FRAME( q_new		)	;
	Real_t *bvc			=FRAME( bvc			)	;
	Real_t *pbvc		=FRAME( pbvc		)	;


	EvalEOSForElemsP1_darts(*domain,vnewc,p_new,e_new,q_new,pbvc,bvc,numElemReg,regElemList,rep);
	//loop to add load imbalance based on region number 
	for(Int_t j = 0; j < rep; j++) {
	   /* compress data, minimal set */
	//	INVOKE(EvalEOSForElemsP1TP,domain,vnewc,p_new,e_new,q_new,pbvc,bvc,numElemReg,regElemList ,&FRAME(evalEOSForElemsP1Sync));
		
		SYNC(evalEOSForElemsP1Sync);
	}


	EXIT_TP();
}

void EvalEOSForElemsP1SyncCD::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsTP);
	for(size_t i=0;i<N_CORES;++i){
		SYNC(evalEOSForElemsP2Iter[i]  );
	}
	EXIT_TP();
}

void EvalEOSForElemsP2IterCD ::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsTP);
	Domain *domain=FRAME(domain);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *q_new		=FRAME( q_new		)	;
	
	size_t	Chunk = numElemReg/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
    
	EvalEOSForElemsP2_darts(*domain,regElemList,p_new,e_new,q_new,lw,hi) ;
	
	SYNC( calcSoundSpeedForElems);

	EXIT_TP();
}



void  CalcSoundSpeedForElemsCD ::fire(void)
{
	
	LOAD_FRAME(EvalEOSForElemsTP);

	Domain *domain	=FRAME(domain);
	Real_t *vnewc   = FRAME(vnew);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 

	Real_t rho0			=FRAME( rho0		)	;
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *bvc			=FRAME( bvc			)	;
	Real_t *pbvc		=FRAME( pbvc		)	;
	
	INVOKE(CalcSoundSpeedForElemsTP,domain,vnewc, rho0, e_new, p_new,pbvc, bvc,numElemReg, regElemList,&FRAME( evalEOSForElemsSync )) ;

	EXIT_TP();
}
void EvalEOSForElemsSyncCD::fire(){
	LOAD_FRAME(EvalEOSForElemsTP);
	SIGNAL(signalUp);
	
	EXIT_TP();

}
	
