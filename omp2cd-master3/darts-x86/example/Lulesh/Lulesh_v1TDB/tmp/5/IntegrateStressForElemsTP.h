#ifndef INTEGRATESTRESSFORELEMSTP_H
#define INTEGRATESTRESSFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(IntegrateStressForElemsP1IterCD,0,SHORTWAIT);
DEF_CODELET(IntegrateStressForElemsP2CD,1,SHORTWAIT );

DEF_TP( IntegrateStressForElemsTP)
{
	Domain *domain;
    Real_t *sigxx  ;
	Real_t *sigyy  ;
	Real_t *sigzz  ;
    Real_t *determ ;
	Index_t numElem;	
	Index_t numNode;

	Real_t *fx_elem;
	Real_t *fy_elem;
	Real_t *fz_elem;
	
	IntegrateStressForElemsP1IterCD *integrateStressForElemsP1Iter;
	IntegrateStressForElemsP2CD  integrateStressForElemsP2;  
	Codelet *signalUp;

	IntegrateStressForElemsTP(Domain *domain,Real_t *sigxx,Real_t *sigyy, Real_t *sigzz,Real_t *determ,Index_t numElem, Index_t numNode,Codelet *up)
		:domain(domain)
		,sigxx(sigxx)
		,sigyy(sigyy)
		,sigzz(sigzz)
		,determ(determ)
		,numElem(numElem)
		,numNode(numNode)
		,integrateStressForElemsP1Iter(new IntegrateStressForElemsP1IterCD[N_CORES])
		,integrateStressForElemsP2(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			Index_t numElem8 = numElem*8 ;

			fx_elem = Allocate<Real_t>(numElem8) ;
			fy_elem = Allocate<Real_t>(numElem8) ;
			fz_elem = Allocate<Real_t>(numElem8) ;
			for (size_t i=0;i<N_CORES;++i){
				integrateStressForElemsP1Iter[i]=IntegrateStressForElemsP1IterCD{0,1,this,SHORTWAIT,i};	
				add ( integrateStressForElemsP1Iter +i);
			}
			add (& integrateStressForElemsP2);
		}

	virtual ~IntegrateStressForElemsTP(){

		delete []integrateStressForElemsP1Iter; 
		Release(&fz_elem) ;
		Release(&fy_elem) ;
		Release(&fx_elem) ;
	}

};


#endif
