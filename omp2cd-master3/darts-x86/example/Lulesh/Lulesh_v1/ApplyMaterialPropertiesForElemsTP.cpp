#include "ApplyMaterialPropertiesForElemsTP.h"
#include "EvalEOSForElemsTP.h"

void ApplyMaterialPropertiesForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(applyMaterialPropertiesForElemsP1Iter[IdL])=ApplyMaterialPropertiesForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdL};	
		ADD ( applyMaterialPropertiesForElemsP1Iter +IdL);
	}
	if(IdR<N_CORES){
		FRAME(applyMaterialPropertiesForElemsP1Iter[IdR])=ApplyMaterialPropertiesForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdR};	
		ADD ( applyMaterialPropertiesForElemsP1Iter +IdR);
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
//    ApplyMaterialPropertiesForElemsP2_darts(*domain,vnew) ;
	EXIT_TP();
}
void ApplyMaterialPropertiesForElemsP2SyncCD::fire(){
	//std::cout<<"EvalEOSForElemsTP finish!"<<std::endl;
	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();

}
	
