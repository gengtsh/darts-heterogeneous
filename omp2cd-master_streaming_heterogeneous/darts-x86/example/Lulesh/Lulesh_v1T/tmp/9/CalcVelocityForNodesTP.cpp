#include "CalcVelocityForNodesTP.h"


void  CalcVelocityForNodesIterCD ::fire(void)
{
	LOAD_FRAME( CalcVelocityForNodesTP);
	size_t	Id	= getID();
//	size_t  IdL = 2*Id+2;
//	size_t  IdR = 2*Id+3;
//	if(IdL<N_CORES){
//		FRAME(calcVelocityForNodesIter[IdL])= CalcVelocityForNodesIterCD{0,1,getTP(),SHORTWAIT,IdL};
//		ADD ( calcVelocityForNodesIter+ IdL);
//	}
//	if(IdR<N_CORES){
//		FRAME(calcVelocityForNodesIter[IdR])= CalcVelocityForNodesIterCD{0,1,getTP(),SHORTWAIT,IdR};
//		ADD ( calcVelocityForNodesIter+ IdR);
//	}

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcVelocityForNodesIter[i])= CalcVelocityForNodesIterCD{0,1,getTP(),SHORTWAIT,i};
		ADD ( calcVelocityForNodesIter+ i);
	}

	Domain *domain =FRAME(domain);	
	Index_t numNode=domain->numNode();	
	
	const Real_t delt = domain->deltatime() ;
	Real_t u_cut = domain->u_cut() ;
	
	size_t	Chunk = numNode/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numNode%N_CORES):0) + (Id+1)* Chunk;

	CalcVelocityForNodes_darts( *domain, delt, u_cut, lw,hi) ;

	SIGNAL(signalUp);

	EXIT_TP();
}
