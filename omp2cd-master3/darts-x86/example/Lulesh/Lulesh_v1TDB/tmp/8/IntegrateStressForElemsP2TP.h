#ifndef INTEGRATESTRESSFORELEMSP2TP_H
#define INTEGRATESTRESSFORELEMSP2TP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(IntegrateStressForElemsP2IterCD,0,SHORTWAIT);

DEF_TP( IntegrateStressForElemsP2TP)
{
	Domain *domain;
	Index_t numNode;

	Real_t *fx_elem;
	Real_t *fy_elem;
	Real_t *fz_elem;
	
	IntegrateStressForElemsP2IterCD *integrateStressForElemsP2Iter;
	Codelet *signalUp;

	IntegrateStressForElemsP2TP(Domain *domain,Index_t numNode,Real_t *fx_elem,Real_t *fy_elem,Real_t *fz_elem,Codelet *up)
		:domain(domain)
		,numNode(numNode)
		,fx_elem(fx_elem)
		,fy_elem(fy_elem)
		,fz_elem(fz_elem)
		,integrateStressForElemsP2Iter(new IntegrateStressForElemsP2IterCD[N_CORES])
		,signalUp(up)
		{
////			for (size_t i=0;i<N_CORES;++i){
////				integrateStressForElemsP2Iter[i]=IntegrateStressForElemsP2IterCD{0,1,this,SHORTWAIT,i};	
////				add ( integrateStressForElemsP2Iter +i);
////			}

			integrateStressForElemsP2Iter[0]=IntegrateStressForElemsP2IterCD{0,1,this,SHORTWAIT,0};	
			add ( integrateStressForElemsP2Iter +0);
			if(N_CORES>1){
				integrateStressForElemsP2Iter[1]=IntegrateStressForElemsP2IterCD{0,1,this,SHORTWAIT,1};	
				add ( integrateStressForElemsP2Iter +1);
			}
		
		}

	virtual ~IntegrateStressForElemsP2TP(){

		delete []integrateStressForElemsP2Iter; 
	}

};


#endif
