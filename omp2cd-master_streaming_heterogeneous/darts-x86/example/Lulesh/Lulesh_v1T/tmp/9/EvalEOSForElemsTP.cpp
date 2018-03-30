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

	INVOKE(EvalEOSForElemsP1TP,rep,domain,vnewc,e_old,delvc,p_old,q_old,compression,compHalfStep,qq_old,ql_old,work,p_new,e_new,q_new,bvc,pbvc,numElemReg,regElemList ,&FRAME(evalEOSForElemsP1Sync));
//	std::cout<<"EvalEOSForElemsP1CD,invoke EvalEOSForElemsP1TP"<<std::endl;
//	std::cout<<"EvalEOSForElemsTP.ref="<<getTP()->getRef()<<",EvalEOSForElemsP1CD,rep="<<rep<<",numZone="<<FRAME(numZone)<<std::endl;

	RESET(evalEOSForElemsP1);
	EXIT_TP();
}

void EvalEOSForElemsP1SyncCD::fire(void)
{
//	std::cout<<"EvalEOSForElemsP1Sync begin!"<<std::endl;	
	LOAD_FRAME(EvalEOSForElemsTP);

	//size_t tree = MIN(g_treeBarrier,N_CORES);
	for ( size_t i = 0; i <N_TREE; ++i ) {
		FRAME(evalEOSForElemsP2Iter[i])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,i};	
		ADD ( evalEOSForElemsP2Iter+i );
	}

	RESET(evalEOSForElemsP1Sync);

	//	std::cout<<"EvalEOSForElemsP1SyncCD!"<<std::endl;
	EXIT_TP();
}

void EvalEOSForElemsP2IterCD ::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsTP);
	size_t	Id	= getID();

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(evalEOSForElemsP2Iter[i])=EvalEOSForElemsP2IterCD{0,1,getTP(),SHORTWAIT,i};	
		ADD ( evalEOSForElemsP2Iter+i );
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
//	std::cout<<"EvalEOSForElemsP2 Iter finish!"<<std::endl;
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

	RESET(evalEOSForElemsSync);
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

	//	std::cout<<"EvalEOSForElemsTP.ref="<<getTP()->getRef()<<",EvalEOSForElemsSyncCD,numZone="<<numZone<<std::endl;
	}
	
	EXIT_TP();
}
	
