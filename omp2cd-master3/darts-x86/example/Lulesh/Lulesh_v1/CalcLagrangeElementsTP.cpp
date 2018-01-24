#include "CalcLagrangeElementsTP.h"
#include "CalcLagrangeElementsP2TP.h"

void CalcKinematicsForElemsIterCD::fire(void)
{

	LOAD_FRAME(CalcLagrangeElementsTP);
	
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcKinematicsForElemsIter[IdL])=CalcKinematicsForElemsIterCD{0,1,getTP(),SHORTWAIT,IdL};
		ADD ( calcKinematicsForElemsIter +IdL);
	}
	if(IdR<N_CORES){
		FRAME(calcKinematicsForElemsIter[IdR])=CalcKinematicsForElemsIterCD{0,1,getTP(),SHORTWAIT,IdR};
		ADD ( calcKinematicsForElemsIter +IdR);
	}
	
	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);
	Index_t numElem=FRAME(numElem); 		
	const Real_t deltatime = domain->deltatime() ;

	size_t	Chunk = numElem/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
    
    CalcKinematicsForElems_darts(*domain, vnew, deltatime, lw,hi) ;
	
	SYNC(calcLagrangeElementsP2);


	EXIT_TP();
}

void  CalcLagrangeElementsP2CD ::fire(void)
{
	
	LOAD_FRAME(CalcLagrangeElementsTP);

	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);

	INVOKE(CalcLagrangeElementsP2TP, domain,vnew,FRAME(signalUp)) ;
	EXIT_TP();
}

	
