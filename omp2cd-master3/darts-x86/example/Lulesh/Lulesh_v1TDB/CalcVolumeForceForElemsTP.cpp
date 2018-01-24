#include "CalcVolumeForceForElemsTP.h"
#include "InitStressTermsForElemsTP.h"
#include "IntegrateStressForElemsTP.h"
#include "CalcHourglassControlForElemsTP.h"

/* Sum contributions to total stress tensor */
void InitStressTermsForElemsCD::fire(void)
{
	std::cout<<"InitStressTermsForElemsTP begin!"<<std::endl;
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	Domain *domain=FRAME(domain);
	Real_t *sigxx =FRAME(sigxx); 
	Real_t *sigyy =FRAME(sigyy);
	Real_t *sigzz =FRAME(sigzz);
	Index_t numElem=FRAME(numElem); 		

	INVOKE(InitStressTermsForElemsTP,domain, sigxx, sigyy, sigzz, numElem, &FRAME(integrateStressForElems)) ;
	
	EXIT_TP();
}

// call elemlib stress integration loop to produce nodal forces from material stresses.
void IntegrateStressForElemsCD ::fire(void)
{

	std::cout<<"IntegrateStressForElemsTP begin!"<<std::endl;

	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *sigxx =FRAME(sigxx); 
	Real_t *sigyy =FRAME(sigyy);
	Real_t *sigzz =FRAME(sigzz);
	Real_t *determ =FRAME(determ);
	Index_t numElem=FRAME(numElem); 		

	INVOKE(IntegrateStressForElemsTP,domain,sigxx, sigyy, sigzz, determ, numElem,domain->numNode(),&FRAME(integrateStressForElemsTPSync)) ;

	EXIT_TP();
}

void  IntegrateStressForElemsTPSyncCD ::fire(void)
{
	std::cout<<"IntegrateStressForElemsTPSync is running!"<<std::endl;
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	//size_t tree = MIN(g_treeBarrier,N_CORES);
	for ( size_t i = 0; i < N_TREE; ++i ) {
		FRAME(checkForNegativeVolumeElemsIter[i])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( checkForNegativeVolumeElemsIter+ i);
	}

	EXIT_TP();
}


void CheckForNegativeVolumeElemsIterCD ::fire(void)
{
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	size_t	Id	= getID();
	
	std::cout<<"CheckForNegativeVolumeElemsIter["<<Id<<"] is running!"<<std::endl;

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(checkForNegativeVolumeElemsIter[i])=CheckForNegativeVolumeElemsIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( checkForNegativeVolumeElemsIter+ i);
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
	std::cout<<"CalcHourglassControlForElemsTP begins!"<<std::endl;
	LOAD_FRAME(CalcVolumeForceForElemsTP);
	
	Domain *domain =FRAME(domain);
	Real_t *determ =FRAME(determ);
	Real_t hgcoef  =FRAME(hgcoef);

	INVOKE(CalcHourglassControlForElemsTP,domain,determ,hgcoef,FRAME(signalUp)) ;

	EXIT_TP();
}
