#ifndef DARTS_EVALEOSFORELEMSP1TP_H
#define DARTS_EVALEOSFORELEMSP1TP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(EvalEOSForElemsP1TPP1IterCD,0,SHORTWAIT);
DEF_CODELET(EvalEOSForElemsP1TPP1SyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(EvalEOSForElemsP1TPP2IterCD,1,SHORTWAIT);

DEF_CODELET(CalcEnergyForElemsCD,1,SHORTWAIT );
DEF_CODELET(EvalEOSForElemsP1TPSyncCD,1,SHORTWAIT );

DEF_TP( EvalEOSForElemsP1TP)
{
	Int_t rep;
	Domain *domain;
    Real_t *vnewc ;

	
	// These temporaries will be of different size for 
	// each call (due to different sized region element
	// lists)
	Real_t *e_old		;
	Real_t *delvc		;
	Real_t *p_old		;
	Real_t *q_old		;
	Real_t *compression ;
	Real_t *compHalfStep;
	Real_t *qq_old		;
	Real_t *ql_old		;
	Real_t *work		;
	Real_t *p_new		;
	Real_t *e_new		;
	Real_t *q_new		;
	Real_t *pbvc		;
	Real_t *bvc			;

	Index_t numElemReg;
	Index_t *regElemList; 
	
	Real_t  e_cut ;
	Real_t  p_cut ;
	// Real_t  ss4o3 ;
	Real_t  q_cut  ;
	
	Real_t eosvmax ;
	Real_t eosvmin ;
	Real_t pmin    ;
	Real_t emin    ;
	Real_t rho0    ;
	Int_t repInit	;
	
	EvalEOSForElemsP1TPP1IterCD *evalEOSForElemsP1TPP1Iter ;
	EvalEOSForElemsP1TPP1SyncCD evalEOSForElemsP1TPP1Sync ;
	EvalEOSForElemsP1TPP2IterCD *evalEOSForElemsP1TPP2Iter ;
	
	CalcEnergyForElemsCD	calcEnergyForElems	;
	EvalEOSForElemsP1TPSyncCD evalEOSForElemsP1TPSync;
	
	
	Codelet *signalUp;


	EvalEOSForElemsP1TP(Int_t rep, Domain *domain,Real_t *vnewc,Real_t *e_old, Real_t *delvc,Real_t *p_old,Real_t *q_old,Real_t *compression,Real_t *compHalfStep,Real_t *qq_old,Real_t *ql_old,Real_t *work,Real_t *p_new,Real_t *e_new,Real_t *q_new,Real_t *bvc,Real_t *pbvc ,Int_t numElemReg, Index_t *regElemList,Codelet *up)
		:rep(rep)
		,domain(domain)
		,vnewc(vnewc)
		,e_old		 ( e_old		)						
		,delvc		 ( delvc		)	
		,p_old		 ( p_old		)	
		,q_old		 ( q_old		)	
		,compression ( compression	)	
		,compHalfStep( compHalfStep )	
		,qq_old		 ( qq_old		)	
		,ql_old		 ( ql_old		)	
		,work		 ( work			)	
		,p_new		 ( p_new		)	
		,e_new		 ( e_new		)	
		,q_new		 ( q_new		)	
		,pbvc		 ( pbvc			)		
		,bvc		 ( bvc			) 
		,numElemReg  (numElemReg	)
		,regElemList (regElemList	)
		,evalEOSForElemsP1TPP1Iter(new EvalEOSForElemsP1TPP1IterCD[N_CORES] )
		,evalEOSForElemsP1TPP1Sync(N_CORES,N_CORES,this,SHORTWAIT ) 
		,evalEOSForElemsP1TPP2Iter(new EvalEOSForElemsP1TPP2IterCD[N_CORES] )
		,calcEnergyForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,evalEOSForElemsP1TPSync(1,1,this,SHORTWAIT)
		,signalUp(up)
		{
			e_cut = domain->e_cut() ;
			p_cut = domain->p_cut() ;
			//ss4o3 = domain->ss4o3() ;
			q_cut = domain->q_cut() ;
			
			eosvmax = domain->eosvmax() ;
			eosvmin = domain->eosvmin() ;
			pmin    = domain->pmin() ;
			emin    = domain->emin() ;
			rho0    = domain->refdens() ;
			
			repInit	= rep;
			// These temporaries will be of different size for 
			// each call (due to different sized region element
			// lists)
			
////			for (size_t i=0;i<N_CORES;++i){
////				evalEOSForElemsP1TPP1Iter[i]=EvalEOSForElemsP1TPP1IterCD{0,1,this,SHORTWAIT,i} ;
////				evalEOSForElemsP1TPP2Iter[i]=EvalEOSForElemsP1TPP2IterCD{1,1,this,SHORTWAIT,i} ;
////				add (evalEOSForElemsP1TPP1Iter +i );
////				add (evalEOSForElemsP1TPP2Iter +i );
////			}
////			add (&evalEOSForElemsP1TPP1Sync);
////			add (&calcEnergyForElems);
////			add (&evalEOSForElemsP1TPSync);


			//size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < N_TREE; ++i ) {
				evalEOSForElemsP1TPP1Iter[i]=EvalEOSForElemsP1TPP1IterCD{0,1,this,SHORTWAIT,i} ;
				add (evalEOSForElemsP1TPP1Iter +i );
			}	
		}

	virtual ~EvalEOSForElemsP1TP(){
//		std::cout<<"~EvalEOSForElemsP1TP--"<<std::endl;
		delete []evalEOSForElemsP1TPP1Iter;
		delete []evalEOSForElemsP1TPP2Iter;
//		std::cout<<"~EvalEOSForElemsP1TP!!"<<std::endl;
	}

};


#endif
