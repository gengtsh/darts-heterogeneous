#include "ApplyMaterialPropertiesForElemsTP.h"
#include "EvalEOSForElemsTP.h"

void ApplyMaterialPropertiesForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	size_t	Id	= getID();
	std::cout<<" ApplyMaterialPropertiesForElemsP1Iter["<<Id<<"] is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(applyMaterialPropertiesForElemsP1Iter[i])=ApplyMaterialPropertiesForElemsP1IterCD{0,1,getTP(),SHORTWAIT,i};	
		ADD ( applyMaterialPropertiesForElemsP1Iter +i);
	}

	Domain *domain=FRAME(domain);
	Real_t *vnew  =FRAME(vnew);
	Index_t numElem=domain->numElem(); 		

	size_t	Chunk = numElem/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElem%N_CORES):0) + (Id+1)* Chunk;
    
	ApplyMaterialPropertiesForElemsP1_darts(*domain,vnew,lw,hi) ;
	
	SYNC(applyMaterialPropertiesForElemsP2);

	EXIT_TP();
}


void  ApplyMaterialPropertiesForElemsP2CD ::fire(void)
{
	std::cout<<"EvalEOSForElemsTP begin!"<<std::endl;
	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);

	Domain *domain	=FRAME(domain);
	Real_t *vnew	=FRAME(vnew);

	Index_t numElemRegMax=0;
	for (Int_t r=0 ; r<domain->numReg() ; r++) {
		Index_t numElemReg = domain->regElemSize(r);
		if (numElemReg>numElemRegMax){
			numElemRegMax = numElemReg;
		}
	}
	INVOKE(EvalEOSForElemsTP,domain, vnew, numElemRegMax, &FRAME(applyMaterialPropertiesForElemsP2Sync));
	EXIT_TP();
}
void ApplyMaterialPropertiesForElemsP2SyncCD::fire(){
	std::cout<<" ApplyMaterialPropertiesForElemsP2Sync is running!"<<std::endl;
	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();

}
	
