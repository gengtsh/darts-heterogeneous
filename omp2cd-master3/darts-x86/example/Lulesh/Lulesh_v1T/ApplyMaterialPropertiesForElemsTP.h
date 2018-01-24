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


			//size_t tree = MIN(g_treeBarrier,N_CORES);
			for ( size_t i = 0; i < N_TREE; ++i ) {
				applyMaterialPropertiesForElemsP1Iter[i]=ApplyMaterialPropertiesForElemsP1IterCD{0,1,this,SHORTWAIT,i};	
				add ( applyMaterialPropertiesForElemsP1Iter +i);
			}
		}

	virtual ~ApplyMaterialPropertiesForElemsTP(){

		delete []applyMaterialPropertiesForElemsP1Iter; 
	}

};


#endif
