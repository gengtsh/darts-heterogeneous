#include "CalcQForElemsTP.h"
#include "CalcMonotonicQRegionForElemsTP.h"

void CalcMonotonicQGradientsForElemsIterCD::fire(void)
{

	LOAD_FRAME(CalcQForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);
	Index_t numElem=FRAME(numElem); 		

	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
    
    CalcMonotonicQGradientsForElems_darts(*domain, vnew, lw,hi) ;
	
	SYNC(calcMonotonicQForElems);

	EXIT_TP();
}

void  CalcMonotonicQForElemsCD ::fire(void)
{
	LOAD_FRAME(CalcQForElemsTP);

	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);
	Real_t ptiny=FRAME(ptiny);
	Index_t numReg=FRAME(numReg);
	
	for(Index_t r=0;r<numReg;++r){
		if(domain->regElemSize(r)>0){
			INVOKE( CalcMonotonicQRegionForElemsTP, domain,r,vnew,ptiny,&FRAME(calcQForElemsSync)) ;
		}
	}
	EXIT_TP();
}

void  CalcQForElemsSyncCD ::fire(void)
{
	LOAD_FRAME(CalcQForElemsTP);
	Domain *domain=FRAME(domain);
	Real_t numElem=FRAME(numElem);

	/* Don't allow excessive artificial viscosity */
	Index_t idx = -1; 
	for (Index_t i=0; i<numElem; ++i) {
	   if ( domain->q(i) > domain->qstop() ) {
	      idx = i ;
	      break ;
	   }
	}
	
	if(idx >= 0) {
	   exit(QStopError);
	}
	
	SIGNAL(signalUp);
	EXIT_TP();
}
	
