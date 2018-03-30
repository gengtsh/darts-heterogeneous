#include "EvalEOSForElemsP1TP.h"
#include "CalcEnergyForElemsTP.h"

void EvalEOSForElemsP1TPP1IterCD ::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsP1TP);
	
	Int_t rep=FRAME(rep);
	size_t	Id	= getID();

	std::cout<<"EvalEOSForElemsP1TPP1Iter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);

	if(rep == FRAME(repInit)){
		for (size_t i=IdC0;i<tree;++i){
			FRAME(evalEOSForElemsP1TPP1Iter[i])=EvalEOSForElemsP1TPP1IterCD{0,1,getTP(),SHORTWAIT,i} ;
			ADD(evalEOSForElemsP1TPP1Iter +i );
		}
	}else{
		for (size_t i=IdC0;i<tree;++i){
			SYNC(evalEOSForElemsP1TPP1Iter[i] );
		}
	}

	Domain *domain = FRAME(domain);
	Real_t *vnewc   = FRAME(vnewc);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 
	                                         
	Real_t *e_old		=FRAME( e_old		)	;
	Real_t *delvc		=FRAME( delvc		)	;
	Real_t *p_old		=FRAME( p_old		)	;
	Real_t *q_old		=FRAME( q_old		)	;
	Real_t *compression	=FRAME( compression)	;
	Real_t *compHalfStep=FRAME( compHalfStep)	;
	Real_t *qq_old		=FRAME( qq_old		)	;
	Real_t *ql_old		=FRAME( ql_old		)	;

	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	SYNC(evalEOSForElemsP1TPP1Sync );	
	EvalEOSForElemsP1TPP1_darts(*domain,vnewc, lw,hi,regElemList,e_old,delvc,p_old,q_old,compression,compHalfStep,qq_old,ql_old	);

	RESET( evalEOSForElemsP1TPP1Iter[Id] );
	

//	std::cout<<"EvalEOSForElemsP1TPP1Iter["<<Id<<"]"<<",ref_="<<getTP()->getRef()<<",rep="<<rep << ",P1TPP1SyncCD.counter="<<FRAME(evalEOSForElemsP1TPP1Sync).getCounter()<<std::endl;	

	EXIT_TP();
}

void EvalEOSForElemsP1TPP1SyncCD ::fire(void)
{
	
	std::cout<<"EvalEOSForElemsP1TPP1Sync is running!"<<std::endl;
	LOAD_FRAME(EvalEOSForElemsP1TP);
	
	Int_t rep=FRAME(rep);
	
	if(rep==FRAME(repInit)){
		for ( size_t i = 0; i < N_TREE; ++i ) {
			FRAME(evalEOSForElemsP1TPP2Iter[i])=EvalEOSForElemsP1TPP2IterCD{0,1,getTP(),SHORTWAIT,i} ;
			ADD(evalEOSForElemsP1TPP2Iter+i );
		}
	}else{
		for ( size_t i = 0; i < N_TREE; ++i ) {
			SYNC(evalEOSForElemsP1TPP2Iter[i] );
		}
	}
		
	
	RESET( evalEOSForElemsP1TPP1Sync);
	
//	std::cout<<"P1TPP1SyncCD!!"<<std::endl;
	EXIT_TP();
}

void EvalEOSForElemsP1TPP2IterCD ::fire(void)
{

	LOAD_FRAME(EvalEOSForElemsP1TP);
	
	Int_t rep=FRAME(rep);
	size_t	Id	= getID();
	
	std::cout<<"EvalEOSForElemsP1TPP2Iter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);

	if(rep==FRAME(repInit)){
		for (size_t i=IdC0;i<tree;++i){
			FRAME(evalEOSForElemsP1TPP2Iter[i])=EvalEOSForElemsP1TPP2IterCD{0,1,getTP(),SHORTWAIT,i} ;
			ADD(evalEOSForElemsP1TPP2Iter+i );
		}
	}else{
		for (size_t i=IdC0;i<tree;++i){
			SYNC(evalEOSForElemsP1TPP2Iter[i] );
		}
	}

	Real_t *vnewc   = FRAME(vnewc);
	
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 

	Real_t eosvmax		=FRAME( eosvmax		)	;
	Real_t eosvmin		=FRAME( eosvmin		)	;
	                                         
	Real_t *p_old		=FRAME( p_old		)	;
	Real_t *compression	=FRAME( compression)	;
	Real_t *compHalfStep=FRAME( compHalfStep)	;
	Real_t *work		=FRAME( work		)	;

	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;
   
	EvalEOSForElemsP1TPP2_darts(vnewc, lw,hi,regElemList,eosvmin,eosvmax,p_old,compression,compHalfStep,work);
	
	SYNC( calcEnergyForElems);

	
	RESET( evalEOSForElemsP1TPP2Iter[Id]);
	

//	std::cout<<"EvalEOSForElemsP1TPP2Iter["<<Id<<"]"<<",ref_="<<getTP()->getRef()<<",rep="<<rep << ",CalcEnergyForElemsCD.counter="<<FRAME(calcEnergyForElems).getCounter()<<std::endl;	

	EXIT_TP();
}



void  CalcEnergyForElemsCD ::fire(void)
{
	
	std::cout<<"CalcEnergyForElemsTP begin!"<<std::endl;
	LOAD_FRAME(EvalEOSForElemsP1TP);

	Real_t *vnewc   = FRAME(vnewc);
	Index_t numElemReg	= FRAME(numElemReg);
	Index_t *regElemList= FRAME(regElemList); 

	Real_t  e_cut		=FRAME(  e_cut		)	;
	Real_t  p_cut		=FRAME(  p_cut		)	;
	// Real_t  ss4o3	=FRAME(  ss4o3		)	;
	Real_t  q_cut		=FRAME(  q_cut		)	;
	                                         
	Real_t eosvmax		=FRAME( eosvmax		)	;
	Real_t pmin			=FRAME( pmin		)	;
	Real_t emin			=FRAME( emin		)	;
	Real_t rho0			=FRAME( rho0		)	;
	                                         
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


//    CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,p_old, e_old,  q_old, compression, compHalfStep,vnewc, work,  delvc, pmin,p_cut, e_cut, q_cut, emin,qq_old, ql_old, rho0, eosvmax,numElemReg, regElemList);
//	SYNC( evalEOSForElemsP1TPSync);

    INVOKE(CalcEnergyForElemsTP,p_new, e_new, q_new, bvc, pbvc,p_old, e_old,  q_old, compression, compHalfStep,vnewc, work,  delvc, pmin,p_cut, e_cut, q_cut, emin,qq_old, ql_old, rho0, eosvmax,numElemReg, regElemList,&FRAME(evalEOSForElemsP1TPSync));
	//std::cout<<"EvalEOSForElemsP1TP,invoke CalcEnergyForElemsTP"<<std::endl;
	
	RESET( calcEnergyForElems);

	EXIT_TP();
}
void EvalEOSForElemsP1TPSyncCD::fire(){
	LOAD_FRAME(EvalEOSForElemsP1TP);
	std::cout<<"EvalEOSForElemsP1TPSync in rep:"<<FRAME(rep)<<",is running!"<<std::endl;

	RESET( evalEOSForElemsP1TPSync);
	Int_t rep= --FRAME(rep);
	if(rep==0){
		SIGNAL(signalUp);
	}else{
	////	for(size_t i=0;i<N_CORES;++i){
	////		SYNC(evalEOSForElemsP1TPP1Iter[i]  );
	////	}
			
			size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < tree; ++i ) {
				SYNC(evalEOSForElemsP1TPP1Iter[i]  );
			}
	}
//	std::cout<<"P1TPSyncCD"<<std::endl;
	EXIT_TP();

}
	
