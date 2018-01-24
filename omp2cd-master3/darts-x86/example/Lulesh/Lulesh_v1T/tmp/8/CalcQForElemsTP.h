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
	
	CalcMonotonicQGradientsForElemsIterCD *calcMonotonicQGradientsForElemsIter;
	CalcMonotonicQForElemsCD calcMonotonicQForElems; 
	CalcQForElemsSyncCD  calcQForElemsSync;  
	Codelet *signalUp;

	CalcQForElemsTP(Domain *domain,Real_t *vnew,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,calcMonotonicQGradientsForElemsIter(new CalcMonotonicQGradientsForElemsIterCD[N_CORES])
		,calcMonotonicQForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,calcQForElemsSync(1,1,this,SHORTWAIT)
		,signalUp(up)
		{
			numElem = domain->numElem();
			numReg  = domain->numReg();	
			allElem = numElem +  /* local elem */
					2*domain->sizeX()*domain->sizeY() + /* plane ghosts */
					2*domain->sizeX()*domain->sizeZ() + /* row ghosts */
					2*domain->sizeY()*domain->sizeZ() ; /* col ghosts */

			domain->AllocateGradients(numElem, allElem);
////			for (size_t i=0;i<N_CORES;++i){
////				calcMonotonicQGradientsForElemsIter[i]=CalcMonotonicQGradientsForElemsIterCD{0,1,this,SHORTWAIT,i};
////				add ( calcMonotonicQGradientsForElemsIter +i);
////			}
////		
////			add (&calcMonotonicQForElems);
////			add (&calcQForElemsSync);

			calcMonotonicQGradientsForElemsIter[0]=CalcMonotonicQGradientsForElemsIterCD{0,1,this,SHORTWAIT,0};
			add ( calcMonotonicQGradientsForElemsIter +0);
			if(N_CORES>1){
				calcMonotonicQGradientsForElemsIter[1]=CalcMonotonicQGradientsForElemsIterCD{0,1,this,SHORTWAIT,1};
				add ( calcMonotonicQGradientsForElemsIter +1);
			}
		}

	virtual ~CalcQForElemsTP(){

		delete []calcMonotonicQGradientsForElemsIter; 
		domain->DeallocateGradients();
	}

};


#endif
