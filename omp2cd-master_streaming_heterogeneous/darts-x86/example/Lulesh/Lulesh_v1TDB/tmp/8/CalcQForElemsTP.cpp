#include "CalcQForElemsTP.h"
#include "CalcMonotonicQRegionForElemsTP.h"

void CalcMonotonicQGradientsForElemsIterCD::fire(void)
{

	LOAD_FRAME(CalcQForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcMonotonicQGradientsForElemsIter[IdL])=CalcMonotonicQGradientsForElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
		ADD ( calcMonotonicQGradientsForElemsIter +IdL);
	}
	if(IdR<N_CORES){
		FRAME(calcMonotonicQGradientsForElemsIter[IdR])=CalcMonotonicQGradientsForElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
		ADD ( calcMonotonicQGradientsForElemsIter +IdR);
	}
	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);
	Index_t numElem=FRAME(numElem); 		

	size_t	Chunk = numElem/N_CORES;
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
	Index_t numReg=FRAME(numReg);
	
	for(Index_t	r=0;r<numReg;++r){
		if(domain->regElemSize(r)>0){
			INVOKE( CalcMonotonicQRegionForElemsTP, domain,vnew,r,&FRAME(calcQForElemsSync)) ;
			break;
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
	
