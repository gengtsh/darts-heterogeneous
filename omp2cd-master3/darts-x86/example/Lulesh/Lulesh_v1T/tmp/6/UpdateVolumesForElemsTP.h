#ifndef INITSTRESSTERMSFORELEMSTP_H
#define INITSTRESSTERMSFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(UpdateVolumesForElemsIterCD,0,SHORTWAIT);

DEF_TP(UpdateVolumesForElemsTP)
{
	Domain *domain;
	Real_t *vnew  ;
	Real_t  v_cut;
	Index_t numElem;	
	
	UpdateVolumesForElemsIterCD *upadteVolumesForElemsIter;
	Codelet *signalUp;
	
	UpdateVolumesForElemsTP(Domain *domain,Real_t *vnew,Real_t v_cut,Index_t numElem,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,v_cut(v_cut)
		,numElem(numElem)
		,upadteVolumesForElemsIter(new UpdateVolumesForElemsIterCD[N_CORES])
		,signalUp(up)
		{
			for ( size_t i = 0; i < N_CORES; ++i ) {
				upadteVolumesForElemsIter[i]= UpdateVolumesForElemsIterCD{0,1,this,SHORTWAIT,i};
				add ( upadteVolumesForElemsIter+ i);
			}
		}
	virtual ~UpdateVolumesForElemsTP(){
		delete [] upadteVolumesForElemsIter;
	}

};


#endif
