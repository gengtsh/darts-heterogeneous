#include "CalcLagrangeElementsTP.h"
#include "CalcLagrangeElementsP2TP.h"

void CalcKinematicsForElemsIterCD::fire(void)
{

	LOAD_FRAME(CalcLagrangeElementsTP);
	
	size_t	Id	= getID();
	std::cout<<" CalcKinematicsForElemsIter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcKinematicsForElemsIter[i])=CalcKinematicsForElemsIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( calcKinematicsForElemsIter +i);
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
	std::cout<<" CalcLagrangeElementsP2TP begin!"<<std::endl;	
	LOAD_FRAME(CalcLagrangeElementsTP);

	Domain *domain=FRAME(domain);
	Real_t *vnew=FRAME(vnew);

	INVOKE(CalcLagrangeElementsP2TP, domain,vnew,FRAME(signalUp)) ;
	EXIT_TP();
}

	
