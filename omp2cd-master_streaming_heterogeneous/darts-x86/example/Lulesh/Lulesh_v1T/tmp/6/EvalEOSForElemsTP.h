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
	Int_t numElemRegMax;
	Int_t numZone;
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

	Real_t *e_old ;
	Real_t *delvc ;
	Real_t *p_old ;
	Real_t *q_old ;
	Real_t *compression ;
	Real_t *compHalfStep ;
	Real_t *qq_old	;
	Real_t *ql_old	;
	Real_t *work	;
	Real_t *p_new	;
	Real_t *e_new	;
	Real_t *q_new	;
	Real_t *bvc		;
	Real_t *pbvc	;
	
	EvalEOSForElemsTP(Domain *domain,Real_t *vnew,Index_t numElemRegMax ,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,numElemRegMax(numElemRegMax)
		//,numElemReg(numElemReg)
		//,regElemList(regElemList)
		//,rep(rep)
		,evalEOSForElemsP1(0,1,this,SHORTWAIT)
		,evalEOSForElemsP1Sync(1,1,this,SHORTWAIT)
		,evalEOSForElemsP2Iter(new EvalEOSForElemsP2IterCD[N_CORES])
		,calcSoundSpeedForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,evalEOSForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			rho0    = domain->refdens() ;
			numZone = 0;	
			numElemReg = domain->regElemSize(0);
			regElemList= domain->regElemlist(0);
			if (numZone<domain->numReg()/2){
				rep = 1;
			}else if(numZone<(domain->numReg()-(domain->numReg()+15)/20)){
				rep = 1 +domain->cost();
			}else{
				rep = 10 * (1+domain->cost());
			}
			// These temporaries will be of different size for 
			// each call (due to different sized region element
			// lists)
	
			
			e_old = Allocate<Real_t>(numElemRegMax) ;
			delvc = Allocate<Real_t>(numElemRegMax) ;
			p_old = Allocate<Real_t>(numElemRegMax) ;
			q_old = Allocate<Real_t>(numElemRegMax) ;
			compression = Allocate<Real_t>(numElemRegMax) ;
			compHalfStep = Allocate<Real_t>(numElemRegMax) ;
			qq_old = Allocate<Real_t>(numElemRegMax) ;
			ql_old = Allocate<Real_t>(numElemRegMax) ;
			work = Allocate<Real_t>(numElemRegMax) ;
			
			p_new = Allocate<Real_t>(numElemRegMax) ;
			e_new = Allocate<Real_t>(numElemRegMax) ;
			q_new = Allocate<Real_t>(numElemRegMax) ;
			bvc = Allocate<Real_t>(numElemRegMax) ;
			pbvc = Allocate<Real_t>(numElemRegMax) ;
		
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
//		std::cout<<"~EvalEOSForElemsTP--"<<std::endl;
		Release(&work) ;
		Release(&ql_old) ;
		Release(&qq_old) ;
		Release(&compHalfStep) ;
		Release(&compression) ;
		Release(&q_old) ;
		Release(&p_old) ;
		Release(&delvc) ;
		Release(&e_old) ;
		
		Release(&pbvc) ;
		Release(&bvc) ;
		Release(&q_new) ;
		Release(&e_new) ;
		Release(&p_new) ;
		delete []evalEOSForElemsP2Iter;
//		std::cout<<"~EvalEOSForElemsTP!!"<<std::endl;
	}

};


#endif
