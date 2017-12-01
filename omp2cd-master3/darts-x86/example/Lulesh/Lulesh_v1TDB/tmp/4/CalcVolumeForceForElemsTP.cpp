#include "CalcVolumeForceForElemsTP.h"
#include "InitStressTermsForElemsTP.h"
#include "IntegrateStressForElemsTP.h"
#include "CalcHourglassControlForElemsTP.h"

/* Sum contributions to total stress tensor */
void InitStressTermsForElemsCD::fire(void)
{

	LOAD_FRAME(CalcVolumeForceForElemsTP);
	Domain *domain=FRAME(domain);
	Real_t *sigxx =FRAME(sigxx); 
	Real_t *sigyy =FRAME(sigyy);
	Real_t *sigzz =FRAME(sigzz);
	Index_t numElem=FRAME(numElem); 		

    //InitStressTermsForElems(*domain, sigxx, sigyy, sigzz, numElem);
	
	//SYNC(integrateStressForElems);

	INVOKE(InitStressTermsForElemsTP,domain, sigxx, sigyy, sigzz, numElem, &FRAME(integrateStressForElems)) ;
	
	EXIT_TP();
}

// call elemlib stress integration loop to produce nodal forces from material stresses.
void IntegrateStressForElemsCD ::fire(void)
{

	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *sigxx =FRAME(sigxx); 
	Real_t *sigyy =FRAME(sigyy);
	Real_t *sigzz =FRAME(sigzz);
	Real_t *determ =FRAME(determ);
	Index_t numElem=FRAME(numElem); 		
	
//    IntegrateStressForElems(*domain,sigxx, sigyy, sigzz, determ, numElem,domain->numNode()) ;
//	for (size_t i=0;i<N_CORES;i++){
//		SYNC(checkForNegativeVolumeElemsIter[i]);
//	}

	INVOKE(IntegrateStressForElemsTP,domain,sigxx, sigyy, sigzz, determ, numElem,domain->numNode(),&FRAME(integrateStressForElemsTPSync)) ;

	EXIT_TP();
}

void  IntegrateStressForElemsTPSyncCD ::fire(void)
{

	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	for (size_t i=0;i<N_CORES;i++){
		SYNC(checkForNegativeVolumeElemsIter[i]);
	}

	EXIT_TP();
}


void CheckForNegativeVolumeElemsIterCD ::fire(void)
{
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	Real_t *determ =FRAME(determ);
	Index_t numElem=FRAME(numElem); 		
	
	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
	
	for ( Index_t k=lw ; k<hi ; ++k ) {
		if (determ[k] <= Real_t(0.0)) {
            exit(VolumeError);
		}
	}
	
	SYNC(calcHourglassControlForElems);

	//INVOKE(TP, domain,&FRAME(signalUp)) ;

	EXIT_TP();
}
void CalcHourglassControlForElemsCD ::fire(void)
{

	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	Domain *domain =FRAME(domain);
	Real_t *determ =FRAME(determ);
	Real_t hgcoef  =FRAME(hgcoef);

	//CalcHourglassControlForElems(*domain, determ, hgcoef) ;
	//SIGNAL(signalUp);

	INVOKE(CalcHourglassControlForElemsTP,domain,determ,hgcoef,FRAME(signalUp)) ;

	EXIT_TP();
}
