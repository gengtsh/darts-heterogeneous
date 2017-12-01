#ifndef DARTS_CALCQFORELEMSTP_H
#define DARTS_CALCQFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcMonotonicQGradientsForElemsIterCD,0,SHORTWAIT);
DEF_CODELET(CalcMonotonicQForElemsCD,1,SHORTWAIT );
DEF_CODELET(CalcQForElemsSyncCD,1,SHORTWAIT );

DEF_TP(CalcQForElemsTP)
{
	Domain *domain;
	Real_t *vnew;
	Index_t numElem;	
	Int_t allElem;
	Index_t numReg;
	const Real_t ptiny = Real_t(1.e-36) ;
	
	CalcMonotonicQGradientsForElemsIterCD *calcMonotonicQGradientsForElemsIter;
	CalcMonotonicQForElemsCD calcMonotonicQForElems; 
	CalcQForElemsSyncCD  calcQForElemsSync;  
	Codelet *signalUp;

	CalcQForElemsTP(Domain *domain,Real_t *vnew,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,calcMonotonicQGradientsForElemsIter(new CalcMonotonicQGradientsForElemsIterCD[N_CORES])
		,calcMonotonicQForElems(N_CORES,N_CORES,this,SHORTWAIT)
		//,calcQForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			numElem = domain->numElem();
			numReg  = domain->numReg();	
			allElem = numElem +  /* local elem */
					2*domain->sizeX()*domain->sizeY() + /* plane ghosts */
					2*domain->sizeX()*domain->sizeZ() + /* row ghosts */
					2*domain->sizeY()*domain->sizeZ() ; /* col ghosts */

			domain->AllocateGradients(numElem, allElem);
			for (size_t i=0;i<N_CORES;++i){
				calcMonotonicQGradientsForElemsIter[i]=CalcMonotonicQGradientsForElemsIterCD{0,1,this,SHORTWAIT,i};
				add ( calcMonotonicQGradientsForElemsIter +i);
			}
		
			uint32_t numLargeZero =0;
			for (Index_t r=0;r<numReg;++r){
				if(domain->regElemSize(r)>0){
					++numLargeZero;
				}
			}
			calcQForElemsSync = CalcQForElemsSyncCD{numLargeZero,numLargeZero,this,SHORTWAIT};
			add (&calcMonotonicQForElems);
			add (&calcQForElemsSync);
		}

	virtual ~CalcQForElemsTP(){

		delete []calcMonotonicQGradientsForElemsIter; 
		domain->DeallocateGradients();
	}

};


#endif
