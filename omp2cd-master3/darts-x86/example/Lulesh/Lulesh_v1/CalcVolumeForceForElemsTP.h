#ifndef CALCVOLUMEFORCEFORELEMS_H
#define CALCVOLUMEFORCEFORELEMS_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET(InitStressTermsForElemsCD,0,SHORTWAIT);
DEF_CODELET(IntegrateStressForElemsCD,1,SHORTWAIT);
DEF_CODELET(IntegrateStressForElemsTPSyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(CheckForNegativeVolumeElemsIterCD,1,SHORTWAIT);
DEF_CODELET(CalcHourglassControlForElemsCD,1,SHORTWAIT);

DEF_TP(CalcVolumeForceForElemsTP)
{
	Domain *domain;
	Index_t numElem;	
	Real_t  hgcoef ;
    Real_t *sigxx  ;
	Real_t *sigyy  ;
	Real_t *sigzz  ;
    Real_t *determ ;

	InitStressTermsForElemsCD initStressTermsForElems;
	IntegrateStressForElemsCD integrateStressForElems;
	IntegrateStressForElemsTPSyncCD  integrateStressForElemsTPSync;  
	CheckForNegativeVolumeElemsIterCD *checkForNegativeVolumeElemsIter;
	CalcHourglassControlForElemsCD calcHourglassControlForElems;
	Codelet *signalUp;
	CalcVolumeForceForElemsTP(Domain *domain,Codelet *up)
		:domain(domain)
		,initStressTermsForElems(0,1,this,SHORTWAIT)
		,integrateStressForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,integrateStressForElemsTPSync(N_CORES,N_CORES,this,SHORTWAIT)
		,checkForNegativeVolumeElemsIter(new CheckForNegativeVolumeElemsIterCD[N_CORES])
		,calcHourglassControlForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			numElem = domain->numElem() ;
			hgcoef = domain->hgcoef() ;
			sigxx  = Allocate<Real_t>(numElem) ;
			sigyy  = Allocate<Real_t>(numElem) ;
			sigzz  = Allocate<Real_t>(numElem) ;
			determ = Allocate<Real_t>(numElem) ;

			add (& initStressTermsForElems);
////			add (& integrateStressForElems);
////			for ( size_t i = 0; i < N_CORES; ++i ) {
////				checkForNegativeVolumeElemsIter[i]=CheckForNegativeVolumeElemsIterCD{1,1,this,SHORTWAIT,i};
////				add ( checkForNegativeVolumeElemsIter+ i);
////			}
////
////			add (& calcHourglassControlForElems);
		}
	virtual ~CalcVolumeForceForElemsTP(){
		delete [] checkForNegativeVolumeElemsIter;
	
		Release(&determ) ;
		Release(&sigzz) ;
		Release(&sigyy) ;
		Release(&sigxx) ;
	}

};


#endif
