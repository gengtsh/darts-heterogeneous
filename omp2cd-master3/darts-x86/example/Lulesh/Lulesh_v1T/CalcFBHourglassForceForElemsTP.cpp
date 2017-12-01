#define __STDC_FORMAT_MACROS
#include <cstdint>
//#include <cmath>
#include <inttypes.h>
#include "CalcFBHourglassForceForElemsTP.h"


void 
CalcFBHourglassForceForElemsP1IterCD::fire(void) 
{
	LOAD_FRAME(CalcFBHourglassForceForElemsTP);
	
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcFBHourglassForceForElemsP1Iter[IdL]) = CalcFBHourglassForceForElemsP1IterCD {0,1,getTP(),SHORTWAIT,IdL};
//		ADD( calcFBHourglassForceForElemsP1Iter+IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(calcFBHourglassForceForElemsP1Iter[IdR]) = CalcFBHourglassForceForElemsP1IterCD {0,1,getTP(),SHORTWAIT,IdR};
//		ADD(calcFBHourglassForceForElemsP1Iter+IdR);
//	}

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcFBHourglassForceForElemsP1Iter[i]) = CalcFBHourglassForceForElemsP1IterCD {0,1,getTP(),SHORTWAIT,i};
		ADD(calcFBHourglassForceForElemsP1Iter+i);
	}
	
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
	
	Real_t *fx_elem = FRAME(fx_elem); 
	Real_t *fy_elem = FRAME(fy_elem); 
	Real_t *fz_elem = FRAME(fz_elem); 

	Real_t (*gamma)[8]  = FRAME(gamma);

	Index_t	lw;
	Index_t	hi;
	size_t	Chunk = numElem/N_CORES;

	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;

   CalcFBHourglassForceForElems_p1( *locDom,determ, x8n, y8n, z8n, dvdx, dvdy, dvdz,hgcoef, fx_elem,fy_elem,fz_elem,gamma,lw,hi) ;
	
	SYNC(calcFBHourglassForceForElemsP2);

	//std::cout<<"cpeSyncSub["<<Id/nGrain<<"]:="<<FRAME(cpeSyncSub)->getCounter()<<std::endl;

	EXIT_TP();
}

void 
CalcFBHourglassForceForElemsP2CD::fire(void) 
{
	LOAD_FRAME(CalcFBHourglassForceForElemsTP);
////	for (size_t i=0;i<N_CORES;++i){
////		SYNC(calcFBHourglassForceForElemsP2Iter[i]);
////	}


	//size_t tree = MIN(g_treeBarrier,N_CORES);
	for ( size_t i = 0; i < N_TREE; ++i ) {
		FRAME(calcFBHourglassForceForElemsP2Iter[i]) = CalcFBHourglassForceForElemsP2IterCD {0,1,getTP(),SHORTWAIT,i};
		ADD(calcFBHourglassForceForElemsP2Iter+i);
	}
	
	EXIT_TP();
}

void
CalcFBHourglassForceForElemsP2IterCD::fire(void) 
{
	LOAD_FRAME(CalcFBHourglassForceForElemsTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		//SYNC(calcFBHourglassForceForElemsP2Iter[IdL]);
//		FRAME(calcFBHourglassForceForElemsP2Iter[IdL]) = CalcFBHourglassForceForElemsP2IterCD {0,1,getTP(),SHORTWAIT,IdL};
//		ADD(calcFBHourglassForceForElemsP2Iter+IdL);
//	}
//	if(IdR<N_CORES){
//		//SYNC(calcFBHourglassForceForElemsP2Iter[IdR]);
//		FRAME(calcFBHourglassForceForElemsP2Iter[IdR]) = CalcFBHourglassForceForElemsP2IterCD {0,1,getTP(),SHORTWAIT,IdR};
//		ADD(calcFBHourglassForceForElemsP2Iter+IdR);
//	}

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcFBHourglassForceForElemsP2Iter[i]) = CalcFBHourglassForceForElemsP2IterCD {0,1,getTP(),SHORTWAIT,i};
		ADD(calcFBHourglassForceForElemsP2Iter+i);
	}

	Domain *locDom = FRAME(locDom);
	
	Real_t *fx_elem = FRAME(fx_elem); 
	Real_t *fy_elem = FRAME(fy_elem); 
	Real_t *fz_elem = FRAME(fz_elem); 

	Index_t	lw;
	Index_t	hi;
	Index_t numNode=locDom->numNode();
	size_t	Chunk = numNode/N_CORES;

	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;
	CalcFBHourglassForceForElems_p2( *locDom,fx_elem,fy_elem,fz_elem,lw,hi) ;

	SYNC(calcFBHourglassForceForElemsSync);	
	//std::cout<<"chceSyncSub["<<Id/nGrain<<"]:="<<FRAME(chceSyncSub)->getCounter()<<std::endl;
	
	EXIT_TP();
}

void CalcFBHourglassForceForElemsSyncCD::fire(void)
{

	LOAD_FRAME(CalcFBHourglassForceForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();
}
