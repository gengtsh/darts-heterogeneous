#ifndef DARTS_EVALEOSFORELEMSTP_H
#define DARTS_EVALEOSFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET(EvalEOSForElemsP1CD,0,SHORTWAIT );
DEF_CODELET(EvalEOSForElemsP1SyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(EvalEOSForElemsP2IterCD,1,SHORTWAIT);
DEF_CODELET(CalcSoundSpeedForElemsCD,1,SHORTWAIT );
DEF_CODELET(EvalEOSForElemsSyncCD,1,SHORTWAIT );

DEF_TP( EvalEOSForElemsTP)
{
	Domain *domain;
    Real_t *vnew ;
	Index_t numElemReg;
	Index_t *regElemList; 
	Int_t rep;

	EvalEOSForElemsP1CD evalEOSForElemsP1;
	EvalEOSForElemsP1SyncCD evalEOSForElemsP1Sync;
	EvalEOSForElemsP2IterCD *evalEOSForElemsP2Iter;
	CalcSoundSpeedForElemsCD calcSoundSpeedForElems;  
	EvalEOSForElemsSyncCD evalEOSForElemsSync;
	Codelet *signalUp;

	
	// These temporaries will be of different size for 
	// each call (due to different sized region element
	// lists)
	Real_t rho0;
	Real_t *p_new	;
	Real_t *e_new	;
	Real_t *q_new	;
	Real_t *bvc		;
	Real_t *pbvc	;

	EvalEOSForElemsTP(Domain *domain,Real_t *vnew,Int_t numElemReg, Index_t *regElemList, Int_t rep,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,numElemReg(numElemReg)
		,regElemList(regElemList)
		,rep(rep)
		,evalEOSForElemsP1(0,1,this,SHORTWAIT)
		,evalEOSForElemsP1Sync(rep,rep,this,SHORTWAIT)
		,evalEOSForElemsP2Iter(new EvalEOSForElemsP2IterCD[N_CORES])
		,calcSoundSpeedForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,evalEOSForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			rho0    = domain->refdens() ;
			
			// These temporaries will be of different size for 
			// each call (due to different sized region element
			// lists)
			p_new = Allocate<Real_t>(numElemReg) ;
			e_new = Allocate<Real_t>(numElemReg) ;
			q_new = Allocate<Real_t>(numElemReg) ;
			bvc = Allocate<Real_t>(numElemReg) ;
			pbvc = Allocate<Real_t>(numElemReg) ;
		
			add (& evalEOSForElemsP1);
			add (& evalEOSForElemsP1Sync);
			for (size_t i=0;i<N_CORES;++i){
				evalEOSForElemsP2Iter[i]=EvalEOSForElemsP2IterCD{1,1,this,SHORTWAIT,i};	
				add ( evalEOSForElemsP2Iter+i );
			}
			add (&calcSoundSpeedForElems);
			add (&evalEOSForElemsSync);
		}

	virtual ~EvalEOSForElemsTP(){

		Release(&pbvc) ;
		Release(&bvc) ;
		Release(&q_new) ;
		Release(&e_new) ;
		Release(&p_new) ;
		delete []evalEOSForElemsP2Iter;
	}

};


#endif
