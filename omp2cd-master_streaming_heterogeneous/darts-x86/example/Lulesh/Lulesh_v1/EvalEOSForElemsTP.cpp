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
	                                       
	Real_t *e_old		=FRAME( e_old		)	;
	Real_t *delvc		=FRAME( delvc		)	;
	Real_t *p_old		=FRAME( p_old		)	;
	Real_t *q_old		=FRAME( q_old		)	;
	Real_t *compression	=FRAME( compression)	;
	Real_t *compHalfStep=FRAME( compHalfStep)	;
	Real_t *qq_old		=FRAME( qq_old		)	;
	Real_t *ql_old		=FRAME( ql_old		)	;
	Real_t *work		=FRAME( work		)	;
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *q_new		=FRAME( q_new		)	;
	Real_t *bvc			=FRAME( bvc			)	;
	Real_t *pbvc		=FRAME( pbvc		)	;

//	//EvalEOSForElemsP1_darts(*domain,vnewc,p_new,e_new,q_new,pbvc,bvc,numElemReg,regElemList,rep);
//	//loop to add load imbalance based on region number 
//	for(Int_t j = 0; j < rep; j++) {
//	   /* compress data, minimal set */
//		INVOKE(EvalEOSForElemsP1TP,domain,vnewc,p_new,e_new,q_new,pbvc,bvc,numElemReg,regElemList ,&FRAME(evalEOSForElemsP1Sync));
//		
//	//	SYNC(evalEOSForElemsP1Sync);
//	}
	INVOKE(EvalEOSForElemsP1TP,rep,domain,vnewc,e_old,delvc,p_old,q_old,compression,compHalfStep,qq_old,ql_old,work,p_new,e_new,q_new,bvc,pbvc,numElemReg,regElemList ,&FRAME(evalEOSForElemsP1Sync));
//	std::cout<<"EvalEOSForElemsP1CD,invoke EvalEOSForElemsP1TP"<<std::endl;
	//Int_t numZone = FRAME(numZone);
	//if(numZone<domain->numReg()-1){
	//	RESET(evalEOSForElemsP1);
	//	getTP()->incRef();
	//}
//	std::cout<<"EvalEOSForElemsTP.ref="<<getTP()->getRef()<<",EvalEOSForElemsP1CD,rep="<<rep<<",numZone="<<FRAME(numZone)<<std::endl;

	RESET(evalEOSForElemsP1);
	EXIT_TP();
}

void EvalEOSForElemsP1SyncCD::fire(void)
{
	LOAD_FRAME(EvalEOSForElemsTP);
////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(evalEOSForElemsP2Iter[i]  );
////	}

	FRAME(evalEOSForElemsP2Iter[0])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,0};	
	ADD ( evalEOSForElemsP2Iter+0 );
	if(N_CORES>1){
		FRAME(evalEOSForElemsP2Iter[1])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,1};	
		ADD ( evalEOSForElemsP2Iter+1 );
	}

	RESET(evalEOSForElemsP1Sync);
//	std::cout<<"EvalEOSForElemsP1SyncCD!"<<std::endl;
	EXIT_TP();
}

void EvalEOSForElemsP2IterCD ::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(evalEOSForElemsP2Iter[IdL])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,IdL};	
		ADD ( evalEOSForElemsP2Iter+IdL );
	}
	if(IdR<N_CORES){
		FRAME(evalEOSForElemsP2Iter[IdR])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,IdR};	
		ADD ( evalEOSForElemsP2Iter+IdR );
	}
	Domain *domain=FRAME(domain);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 
	Real_t *p_new		=FRAME( p_new		)	;
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *q_new		=FRAME( q_new		)	;
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
    
	EvalEOSForElemsP2_darts(*domain,regElemList,p_new,e_new,q_new,lw,hi) ;
	
//	std::cout<<"EvalEOSForElemsP2IterCD, numElemReg="<<FRAME(numElemReg)<<std::endl;
	SYNC( calcSoundSpeedForElems);
	
	RESET(evalEOSForElemsP2Iter[Id]);
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
	//std::cout<<"invoke CalcSoundSpeedForElemsTP"<<std::endl;	
	
	RESET(calcSoundSpeedForElems);
	EXIT_TP();
}
void EvalEOSForElemsSyncCD::fire(){
	LOAD_FRAME(EvalEOSForElemsTP);
	Domain *domain = FRAME(domain);
	Int_t numZone  = ++FRAME(numZone);	
	if(numZone == domain->numReg()){
		SIGNAL(signalUp);
		//unsigned int ref=getTP()->getRef();
		//for(size_t i=0;i<ref-1;++i){
		//	getTP()->decRef();
		//}
	}else if(numZone < domain->numReg()){
		FRAME(numElemReg) = domain->regElemSize(numZone);
		FRAME(regElemList)= domain->regElemlist(numZone);
		if (numZone<domain->numReg()/2){
			FRAME(rep) = 1;
		}else if(numZone<(domain->numReg()-(domain->numReg()+15)/20)){
			FRAME(rep) = 1 +domain->cost();
		}else{
			FRAME(rep) = 10 * (1+domain->cost());
		}
		SYNC(evalEOSForElemsP1);
		RESET(evalEOSForElemsSync);

	//	std::cout<<"EvalEOSForElemsTP.ref="<<getTP()->getRef()<<",EvalEOSForElemsSyncCD,numZone="<<numZone<<std::endl;
	}
	
	EXIT_TP();
}
	
