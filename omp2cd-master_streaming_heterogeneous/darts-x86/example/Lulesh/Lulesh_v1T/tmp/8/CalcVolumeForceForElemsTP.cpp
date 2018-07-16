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
	
////	for (size_t i=0;i<N_CORES;i++){
////		SYNC(checkForNegativeVolumeElemsIter[i]);
////	}

	FRAME(checkForNegativeVolumeElemsIter[0])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,0};
	ADD ( checkForNegativeVolumeElemsIter+ 0);
	if(N_CORES>1){
		FRAME(checkForNegativeVolumeElemsIter[1])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,1};
		ADD ( checkForNegativeVolumeElemsIter+ 1);
	}
	EXIT_TP();
}


void CheckForNegativeVolumeElemsIterCD ::fire(void)
{
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(checkForNegativeVolumeElemsIter[IdL])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
		ADD ( checkForNegativeVolumeElemsIter+ IdL);
	}
	if(IdR<N_CORES){
		FRAME(checkForNegativeVolumeElemsIter[IdR])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
		ADD ( checkForNegativeVolumeElemsIter+ IdR);
	}
	Real_t *determ =FRAME(determ);
	Index_t numElem=FRAME(numElem); 		
	
	size_t	Chunk = numElem/N_CORES;
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
