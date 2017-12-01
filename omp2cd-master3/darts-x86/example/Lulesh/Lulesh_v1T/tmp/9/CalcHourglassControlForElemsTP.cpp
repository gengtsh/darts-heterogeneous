#define __STDC_FORMAT_MACROS
#include <cstdint>
//#include <cmath>
#include <inttypes.h>

#include "CalcHourglassControlForElemsTP.h"
#include "CalcFBHourglassForceForElemsTP.h"

void 
CalcHourglassControlForElemsIterCD::fire(void) 
{
	LOAD_FRAME(CalcHourglassControlForElemsTP);
	
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcHourglassControlForElemsIter[IdL]) = CalcHourglassControlForElemsIterCD {0,1,getTP(),SHORTWAIT,IdL};
//		ADD( calcHourglassControlForElemsIter + IdL );
//	}
//	if(IdR<N_CORES){
//		FRAME(calcHourglassControlForElemsIter[IdR]) = CalcHourglassControlForElemsIterCD {0,1,getTP(),SHORTWAIT,IdR};
//		ADD( calcHourglassControlForElemsIter + IdR );
//	}
	
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcHourglassControlForElemsIter[i]) = CalcHourglassControlForElemsIterCD {0,1,getTP(),SHORTWAIT,i};
		ADD( calcHourglassControlForElemsIter + i );
	}	

	Domain *locDom = FRAME(locDom);
	Real_t *determ = FRAME(determ);

	Index_t numElem = FRAME(numElem);
		
	Real_t *dvdx	= FRAME(dvdx) ;
	Real_t *dvdy	= FRAME(dvdy) ;
	Real_t *dvdz	= FRAME(dvdz) ;
	Real_t *x8n 	= FRAME(x8n ) ;
	Real_t *y8n 	= FRAME(y8n ) ;
	Real_t *z8n 	= FRAME(z8n ) ;
	
	Index_t	lw;
	Index_t	hi;
	size_t	Chunk = numElem/N_CORES;

	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;

	CalcHourglassControlForElems_p1(*locDom,determ,dvdx,dvdy,dvdz,x8n,y8n,z8n ,lw, hi );
	
	SYNC(calcHourglassControlForElemsSync);

	//std::cout<<"cpeSyncSub["<<Id/nGrain<<"]:="<<FRAME(cpeSyncSub)->getCounter()<<std::endl;

	EXIT_TP();
}


void 
CalcHourglassControlForElemsSyncCD::fire(void) 
{
	LOAD_FRAME(CalcHourglassControlForElemsTP);

	Domain *locDom = FRAME(locDom);
	Real_t *determ = FRAME(determ);
	Real_t hgcoef  = FRAME(hgcoef);
	
	Index_t numElem = FRAME(numElem);
		
	Real_t *dvdx	= FRAME(dvdx) ;
	Real_t *dvdy	= FRAME(dvdy) ;
	Real_t *dvdz	= FRAME(dvdz) ;
	Real_t *x8n 	= FRAME(x8n ) ;
	Real_t *y8n 	= FRAME(y8n ) ;
	Real_t *z8n 	= FRAME(z8n ) ;

	Index_t numNode=locDom->numNode();

	INVOKE(CalcFBHourglassForceForElemsTP, locDom,determ, x8n, y8n, z8n, dvdx, dvdy, dvdz,hgcoef,numElem, numNode,FRAME(signalUp) ) ;

	EXIT_TP();
}


