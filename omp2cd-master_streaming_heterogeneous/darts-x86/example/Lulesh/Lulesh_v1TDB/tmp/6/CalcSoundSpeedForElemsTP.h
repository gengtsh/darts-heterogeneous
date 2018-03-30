#ifndef DARTS_CALCSOUNDSPEEDFORELEMSTP_H
#define DARTS_CALCSOUNDSPEEDFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

DEF_CODELET_ITER(CalcSoundSpeedForElemsIterCD,0,SHORTWAIT);

DEF_TP(CalcSoundSpeedForElemsTP)
{
	Domain *domain;
	Real_t *vnewc   ;
	Real_t rho0			;
	Real_t *e_new		;
	Real_t *p_new		;
	Real_t *pbvc		;
	Real_t *bvc			;
	Index_t numElemReg	;
	Index_t *regElemList; 
	
	CalcSoundSpeedForElemsIterCD *calcSoundSpeedForElemsIter;
	Codelet *signalUp;
	
	CalcSoundSpeedForElemsTP(Domain *domain,Real_t * vnewc, Real_t rho0, Real_t *e_new, Real_t *p_new,Real_t * pbvc, Real_t *bvc,Index_t numElemReg, Index_t *regElemList, Codelet *up) 
		:domain(domain)
		,vnewc(vnewc)
		,rho0	(rho0	)	
		,e_new	(e_new	)	
		,p_new	(p_new	)	
		,pbvc	(pbvc	)	
		,bvc	(bvc	)	
		,numElemReg(numElemReg)	
		,regElemList(regElemList)
		,calcSoundSpeedForElemsIter(new  CalcSoundSpeedForElemsIterCD[N_CORES])
		,signalUp(up)
		{
			for ( size_t i = 0; i < N_CORES; ++i ) {
				calcSoundSpeedForElemsIter[i]=CalcSoundSpeedForElemsIterCD{0,1,this,SHORTWAIT,i};
				add ( calcSoundSpeedForElemsIter+ i);
			}
		}
	virtual ~CalcSoundSpeedForElemsTP(){
		//std::cout<<"~CalcSoundSpeedForElemsTP--"<<std::endl;
		delete [] calcSoundSpeedForElemsIter;
		//std::cout<<"~CalcSoundSpeedForElemsTP!!"<<std::endl;
	}

};


#endif
