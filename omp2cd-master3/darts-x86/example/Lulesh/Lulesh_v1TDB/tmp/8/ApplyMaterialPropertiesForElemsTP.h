#ifndef DARTS_APPLYMATERIALPROPERTIESFORELEMSTP_H
#define DARTS_APPLYMATERIALPROPERTIESFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(ApplyMaterialPropertiesForElemsP1IterCD,0,SHORTWAIT);

DEF_CODELET(ApplyMaterialPropertiesForElemsP2CD,1,SHORTWAIT );

DEF_CODELET(ApplyMaterialPropertiesForElemsP2SyncCD,1,SHORTWAIT );

DEF_TP( ApplyMaterialPropertiesForElemsTP)
{
	Domain *domain;
    Real_t *vnew ;

	ApplyMaterialPropertiesForElemsP1IterCD *applyMaterialPropertiesForElemsP1Iter;

	ApplyMaterialPropertiesForElemsP2CD  applyMaterialPropertiesForElemsP2;  
	ApplyMaterialPropertiesForElemsP2SyncCD  applyMaterialPropertiesForElemsP2Sync;  
	Codelet *signalUp;

	ApplyMaterialPropertiesForElemsTP(Domain *domain,Real_t *vnew,Codelet *up)
		:domain(domain)
		,vnew(vnew)
		,applyMaterialPropertiesForElemsP1Iter(new ApplyMaterialPropertiesForElemsP1IterCD[N_CORES])
		,applyMaterialPropertiesForElemsP2(N_CORES,N_CORES,this,SHORTWAIT)
		,applyMaterialPropertiesForElemsP2Sync(1,1,this,SHORTWAIT)
		,signalUp(up)
		{
////			for (size_t i=0;i<N_CORES;++i){
////				applyMaterialPropertiesForElemsP1Iter[i]=ApplyMaterialPropertiesForElemsP1IterCD{0,1,this,SHORTWAIT,i};	
////				add ( applyMaterialPropertiesForElemsP1Iter +i);
////			}
////			add (& applyMaterialPropertiesForElemsP2);
////			add(&applyMaterialPropertiesForElemsP2Sync);	

			applyMaterialPropertiesForElemsP1Iter[0]=ApplyMaterialPropertiesForElemsP1IterCD{0,1,this,SHORTWAIT,0};	
			add ( applyMaterialPropertiesForElemsP1Iter +0);
			if(N_CORES>1){
				applyMaterialPropertiesForElemsP1Iter[1]=ApplyMaterialPropertiesForElemsP1IterCD{0,1,this,SHORTWAIT,1};	
				add ( applyMaterialPropertiesForElemsP1Iter +1);
			}
		}

	virtual ~ApplyMaterialPropertiesForElemsTP(){

		delete []applyMaterialPropertiesForElemsP1Iter; 
	}

};


#endif
