#include "ApplyMaterialPropertiesForElemsTP.h"
#include "EvalEOSForElemsTP.h"

void ApplyMaterialPropertiesForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	
	Domain *domain=FRAME(domain);
	Real_t *vnew  =FRAME(vnew);
	Index_t numElem=domain->numElem(); 		

	size_t	Chunk = numElem/N_CORES;
	size_t	Id	= getID();
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

//	for (Int_t r=0 ; r<domain->numReg() ; r++) {
//		Index_t numElemReg = domain->regElemSize(r);
//		Index_t *regElemList = domain->regElemlist(r);
//		Int_t rep;
//		//Determine load imbalance for this region
//		//round down the number with lowest cost
//		if(r < domain->numReg()/2)
//			rep = 1;
//		//you don't get an expensive region unless you at least have 5 regions
//		else if(r < (domain->numReg() - (domain->numReg()+15)/20))
//			rep = 1 + domain->cost();
//		//very expensive regions
//		else
//			rep = 10 * (1+ domain->cost());
//	
//		INVOKE(EvalEOSForElemsTP,domain, vnew, numElemReg, regElemList, rep,&FRAME(applyMaterialPropertiesForElemsP2Sync));
//	}
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

	LOAD_FRAME(ApplyMaterialPropertiesForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();

}
	
