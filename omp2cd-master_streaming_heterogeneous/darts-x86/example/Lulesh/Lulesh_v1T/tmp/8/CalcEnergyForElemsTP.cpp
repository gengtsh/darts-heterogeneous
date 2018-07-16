#include "CalcEnergyForElemsTP.h"
#include "CalcPressureForElemsTP.h"

void CalcEnergyForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcEnergyForElemsP1Iter[IdL])=CalcEnergyForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdL} ;
		ADD (calcEnergyForElemsP1Iter +IdL );
	}
	if(IdR<N_CORES){
		FRAME(calcEnergyForElemsP1Iter[IdR])=CalcEnergyForElemsP1IterCD{0,1,getTP(),SHORTWAIT,IdR} ;
		ADD (calcEnergyForElemsP1Iter +IdR );
	}

	Index_t numElemReg	= FRAME(numElemReg);
	                                         
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *e_old		=FRAME( e_old		)	;
	Real_t *delvc		=FRAME( delvc		)	;
	Real_t *p_old		=FRAME( p_old		)	;
	Real_t *q_old		=FRAME( q_old		)	;
	Real_t *work		=FRAME( work		)	;
	Real_t  emin		=FRAME( emin		)	;

	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP1_darts(e_new, lw,hi,e_old,delvc,p_old,q_old,work,emin);
	
	SYNC(calcEnergyForElemsP2 );	
	
//	std::cout<<"CalcEnergyForElemsP1Iter["<<Id<<"],TP address="<<getTP()<<",ref_="<<getTP()->getRef()<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP2CD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	Real_t  *pHalfStep	  =  FRAME(pHalfStep	 );
	Real_t  *bvc		  =  FRAME(bvc			 );
	Real_t  *pbvc		  =  FRAME(pbvc			 );
	Real_t  *e_new		  =  FRAME(e_new		 );
	Real_t  *compHalfStep =  FRAME(compHalfStep  );
	Real_t  *vnewc		  =  FRAME(vnewc		 );
	Real_t   pmin		  =  FRAME(pmin			 );
	Real_t   p_cut		  =  FRAME(p_cut		 );
	Real_t   eosvmax	  =  FRAME(eosvmax		 );
	Index_t  numElemReg	  =  FRAME(numElemReg	 );
	Index_t *regElemList  =  FRAME(regElemList	 ); 
	
	INVOKE(CalcPressureForElemsTP,pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,pmin, p_cut, eosvmax, numElemReg, regElemList,&FRAME(calcEnergyForElemsP2Sync));
//	std::cout<<"*2,invoke CalcPressureForElemsTP"<<std::endl;
//	std::cout<<"*2,"<<"CalcEnergyForElemsTP address="<<getTP()<<",ref_="<<getTP()->getRef()<<",CalcEnergyForElemsP2SyncCD.address="<<&FRAME(calcEnergyForElemsP2Sync)<<",counter="<<FRAME(calcEnergyForElemsP2Sync).getCounter()<<std::endl;	
	EXIT_TP();
}

void CalcEnergyForElemsP2SyncCD ::fire(void)
{
	
//	std::cout<<"CalcEnergyForElemsP2Sync"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	//std::cout<<"CalcEnergyForElemsTP.ref_="<<getTP()->getRef()<<",P3Iter[0].counter="<<FRAME(calcEnergyForElemsP3Iter[0]).getCounter()<<std::endl;
////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP3Iter[i]  );
////	}
	FRAME(calcEnergyForElemsP3Iter[0])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,0} ;
	ADD (calcEnergyForElemsP3Iter +0 );
	if(N_CORES>1){
		FRAME(calcEnergyForElemsP3Iter[1])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,1} ;
		ADD(calcEnergyForElemsP3Iter +1 );
	}
	
//	std::cout<<"2"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP3IterCD ::fire(void)
{
//	std::cout<<"CalcEnergyForElemsP3Iter"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcEnergyForElemsP3Iter[IdL])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,IdL} ;
		ADD (calcEnergyForElemsP3Iter +IdL );
	}
	if(IdR<N_CORES){
		FRAME(calcEnergyForElemsP3Iter[IdR])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,IdR} ;
		ADD (calcEnergyForElemsP3Iter +IdR );
	}

	Index_t numElemReg	= FRAME(numElemReg);

	Real_t  *compHalfStep =  FRAME(compHalfStep );
	Real_t  *delvc		  =  FRAME(delvc		);
	Real_t  *q_new		  =  FRAME(q_new		);
	Real_t  *pbvc		  =  FRAME(pbvc			);
	Real_t  *e_new		  =  FRAME(e_new		);
	Real_t  *pHalfStep	  =  FRAME(pHalfStep	);
	Real_t  *bvc		  =  FRAME(bvc			);
	Real_t  *qq_old		  =  FRAME(qq_old		);
	Real_t  *ql_old		  =  FRAME(ql_old		);
	Real_t  *p_old		  =  FRAME(p_old		);
	Real_t  *q_old		  =  FRAME(q_old		);
	Real_t   rho0		  =  FRAME(rho0			);
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP3_darts(e_new,q_new, lw,hi,delvc,p_old,q_old,pHalfStep, compHalfStep,bvc,pbvc,ql_old,qq_old,rho0);
	
////	SYNC(calcEnergyForElemsP4Iter[Id] );	
	FRAME(calcEnergyForElemsP4Iter[Id])=CalcEnergyForElemsP4IterCD{0,1,getTP(),SHORTWAIT,Id} ;
	ADD(calcEnergyForElemsP4Iter+Id);
	
	EXIT_TP();
}


void CalcEnergyForElemsP4IterCD ::fire(void)
{
	//std::cout<<"CalcEnergyForElemsP4Iter"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();
	
	Index_t numElemReg	= FRAME(numElemReg);

	Real_t  *e_new		  =  FRAME(e_new		);
	Real_t  *work		  =  FRAME(work			);
	Real_t   e_cut		  =  FRAME(e_cut		);
	Real_t   emin		  =  FRAME(emin			);
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP4_darts(e_new,lw,hi,work,e_cut,emin);
	
	SYNC(calcEnergyForElemsP5 );	
	EXIT_TP();
}


void CalcEnergyForElemsP5CD ::fire(void)
{
	//std::cout<<"CalcEnergyForElemsP5"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	Real_t  *p_new		  =  FRAME(p_new	 	 );
	Real_t  *bvc		  =  FRAME(bvc			 );
	Real_t  *pbvc		  =  FRAME(pbvc			 );
	Real_t  *e_new		  =  FRAME(e_new		 );
	Real_t  *compression  =  FRAME(compression   );
	Real_t  *vnewc		  =  FRAME(vnewc		 );
	Real_t   pmin		  =  FRAME(pmin			 );
	Real_t   p_cut		  =  FRAME(p_cut		 );
	Real_t   eosvmax	  =  FRAME(eosvmax		 );
	Index_t  numElemReg	  =  FRAME(numElemReg	 );
	Index_t *regElemList  =  FRAME(regElemList	 ); 
	
	INVOKE(CalcPressureForElemsTP,p_new, bvc, pbvc, e_new, compression, vnewc,pmin, p_cut, eosvmax, numElemReg, regElemList,&FRAME( calcEnergyForElemsP5Sync));

//	std::cout<<"*5,invoke CalcPressureForElemsTP"<<std::endl;
//	std::cout<<"*5,"<<"CalcEnergyForElemsTP.ref="<<getTP()->getRef()<<",CalcEnergyForElemsP5SyncCD.address="<<&FRAME(calcEnergyForElemsP5Sync)<<",counter="<<FRAME(calcEnergyForElemsP5Sync).getCounter()<<std::endl;	
	EXIT_TP();
}

void CalcEnergyForElemsP5SyncCD ::fire(void)
{
	//std::cout<<"CalcEnergyForElemsP5Sync"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);

////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP6Iter[i]  );
////	}

	FRAME(calcEnergyForElemsP6Iter[0])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,0} ;
	ADD (calcEnergyForElemsP6Iter +0 );
	if(N_CORES>1){
		FRAME(calcEnergyForElemsP6Iter[1])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,1} ;
		ADD(calcEnergyForElemsP6Iter +1 );
	}


//	std::cout<<"5"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP6IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);

	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcEnergyForElemsP6Iter[IdL])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,IdL} ;
		ADD (calcEnergyForElemsP6Iter +IdL );
	}
	if(IdR<N_CORES){
		FRAME(calcEnergyForElemsP6Iter[IdR])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,IdR} ;
		ADD (calcEnergyForElemsP6Iter +IdR );
	}

	Index_t numElemReg	= FRAME(numElemReg);

	Real_t  *delvc		  =  FRAME(delvc		);
	Real_t  *q_new		  =  FRAME(q_new		);
	Real_t  *p_new		  =  FRAME(p_new		);
	Real_t  *e_new		  =  FRAME(e_new		);
	Real_t  *pHalfStep	  =  FRAME(pHalfStep	);
	Real_t  *pbvc		  =  FRAME(pbvc			);
	Real_t  *bvc		  =  FRAME(bvc			);
	Real_t  *qq_old		  =  FRAME(qq_old		);
	Real_t  *ql_old		  =  FRAME(ql_old		);
	Real_t  *p_old		  =  FRAME(p_old		);
	Real_t  *q_old		  =  FRAME(q_old		);
	Real_t   e_cut		  =  FRAME(e_cut		);
	Real_t   emin		  =  FRAME(emin			);
	Real_t   rho0		  =  FRAME(rho0			);
	Real_t  *vnewc		  =  FRAME(vnewc		);
	Index_t *regElemList  =  FRAME(regElemList	 ); 
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP6_darts(e_new,p_new,q_new,vnewc,regElemList,lw,hi,delvc,p_old,q_old,pHalfStep, bvc,pbvc,ql_old,qq_old,e_cut,emin,rho0);
	
	SYNC(calcEnergyForElemsP7 );	
	EXIT_TP();
}


void CalcEnergyForElemsP7CD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	
	Real_t  *p_new		  =  FRAME(p_new	 	 );
	Real_t  *bvc		  =  FRAME(bvc			 );
	Real_t  *pbvc		  =  FRAME(pbvc			 );
	Real_t  *e_new		  =  FRAME(e_new		 );
	Real_t  *compression  =  FRAME(compression   );
	Real_t  *vnewc		  =  FRAME(vnewc		 );
	Real_t   pmin		  =  FRAME(pmin			 );
	Real_t   p_cut		  =  FRAME(p_cut		 );
	Real_t   eosvmax	  =  FRAME(eosvmax		 );
	Index_t  numElemReg	  =  FRAME(numElemReg	 );
	Index_t *regElemList  =  FRAME(regElemList	 ); 
	
	INVOKE(CalcPressureForElemsTP,p_new, bvc, pbvc, e_new, compression, vnewc,pmin, p_cut, eosvmax, numElemReg, regElemList,&FRAME(calcEnergyForElemsP7Sync ));

//	std::cout<<"*7,invoke CalcPressureForElemsTP"<<std::endl;
//	std::cout<<"*7,"<<"CalcEnergyForElemsTP.ref="<<getTP()->getRef()<<",CalcEnergyForElemsP7SyncCD.address="<<&FRAME(calcEnergyForElemsP7Sync)<<",counter="<<FRAME(calcEnergyForElemsP7Sync).getCounter()<<std::endl;	
	
	EXIT_TP();
}

void CalcEnergyForElemsP7SyncCD ::fire(void)
{
	//std::cout<<"CalcEnergyForElemsP7Sync"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	
////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP8Iter[i]  );
////	}

	FRAME(calcEnergyForElemsP8Iter[0])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,0} ;
	ADD (calcEnergyForElemsP8Iter +0 );
	if(N_CORES>1){
		FRAME(calcEnergyForElemsP8Iter[1])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,1} ;
		ADD(calcEnergyForElemsP8Iter +1 );
	}

//	std::cout<<"7"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP8IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	size_t	Id	= getID();
	size_t  IdL = 2*Id+2;
	size_t  IdR = 2*Id+3;
	if(IdL<N_CORES){
		FRAME(calcEnergyForElemsP8Iter[IdL])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,IdL} ;
		ADD (calcEnergyForElemsP8Iter +IdL );
	}
	if(IdR<N_CORES){
		FRAME(calcEnergyForElemsP8Iter[IdR])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,IdR} ;
		ADD (calcEnergyForElemsP8Iter +IdR );
	}
	
	Index_t numElemReg	= FRAME(numElemReg);

	Real_t  *q_new		  =  FRAME(q_new		);
	Real_t  *p_new		  =  FRAME(p_new		);
	Real_t  *e_new		  =  FRAME(e_new		);
	Real_t  *qq_old		  =  FRAME(qq_old		);
	Real_t  *ql_old		  =  FRAME(ql_old		);
	Real_t   q_cut		  =  FRAME(q_cut		);
	Real_t   rho0		  =  FRAME(rho0			);
	Real_t  *vnewc		  =  FRAME(vnewc		);
	Index_t *regElemList  =  FRAME(regElemList	 ); 
	Real_t  *pbvc		  =  FRAME(pbvc			);
	Real_t  *bvc		  =  FRAME(bvc			);
	Real_t  *delvc		  =  FRAME(delvc		);
	
	size_t	Chunk = numElemReg/N_CORES;
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP8_darts(q_new,e_new,p_new,vnewc,regElemList,lw,hi,bvc,pbvc,delvc,ql_old,qq_old,q_cut,rho0);
	
	SYNC(calcEnergyForElemsSync );	
	EXIT_TP();
}

void CalcEnergyForElemsSyncCD::fire(){
	LOAD_FRAME(CalcEnergyForElemsTP);
//	std::cout<<"****CalcEnergyForElemsSyncCD****"<<std::endl;
	SIGNAL(signalUp);
	EXIT_TP();

}
	
