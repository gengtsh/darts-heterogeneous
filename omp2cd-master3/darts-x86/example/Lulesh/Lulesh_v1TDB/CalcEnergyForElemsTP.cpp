#include "CalcEnergyForElemsTP.h"
#include "CalcPressureForElemsTP.h"

void CalcEnergyForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();
	std::cout<<"CalcEnergyForElemsP1Iter["<<Id<<" is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcEnergyForElemsP1Iter[i])=CalcEnergyForElemsP1IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP1Iter +i);
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
	std::cout<<"CalcEnergyForElemsP2,CalcPressureForElemsTP begin!"<<std::endl;
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
	
	std::cout<<"CalcEnergyForElemsP2Sync is running!"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	//std::cout<<"CalcEnergyForElemsTP.ref_="<<getTP()->getRef()<<",P3Iter[0].counter="<<FRAME(calcEnergyForElemsP3Iter[0]).getCounter()<<std::endl;
////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP3Iter[i]  );
////	}

	//size_t tree = MIN(g_treeBarrier,N_CORES);
	for (size_t i=0;i<N_TREE;++i){
		FRAME(calcEnergyForElemsP3Iter[i])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP3Iter +i );
	}
	
//	std::cout<<"2"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP3IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();

	std::cout<<"CalcEnergyForElemsP3Iter["<<Id<<"] is running!"<<std::endl;

	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcEnergyForElemsP3Iter[i])=CalcEnergyForElemsP3IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP3Iter +i );
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
	LOAD_FRAME(CalcEnergyForElemsTP);
	size_t	Id	= getID();
	std::cout<<"CalcEnergyForElemsP4Iter["<<Id<<"] is running!"<<std::endl;
	
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
	std::cout<<"CalcEnergyForElemsP5, CalcEnergyForElemsTP begin!"<<std::endl;
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
	std::cout<<"CalcEnergyForElemsP5Sync is running!"<<std::endl;
	LOAD_FRAME(CalcEnergyForElemsTP);

////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP6Iter[i]  );
////	}


	//size_t tree = MIN(g_treeBarrier,N_CORES);
	for (size_t i=0;i<N_TREE;++i){
		FRAME(calcEnergyForElemsP6Iter[i])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP6Iter +i );
	}

//	std::cout<<"5"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP6IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);

	size_t	Id	= getID();

	std::cout<<"CalcEnergyForElemsP6Iter["<<Id<<" is running!"<<std::endl;
	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcEnergyForElemsP6Iter[i])=CalcEnergyForElemsP6IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP6Iter +i );
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

	std::cout<<"CalcEnergyForElemsP7, CalcEnergyForElemsTP begin!"<<std::endl;
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

	std::cout<<"CalcEnergyForElemsP7Sync is running!"<<std::endl;

	LOAD_FRAME(CalcEnergyForElemsTP);
	
////	for(size_t i=0;i<N_CORES;++i){
////		SYNC(calcEnergyForElemsP8Iter[i]  );
////	}


//	size_t tree = MIN(g_treeBarrier,N_CORES);
	for (size_t i=0;i<N_TREE;++i){
		FRAME(calcEnergyForElemsP8Iter[i])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP8Iter +i );
	}

//	std::cout<<"7"<<std::endl;
	EXIT_TP();
}

void CalcEnergyForElemsP8IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	
	size_t	Id	= getID();

	std::cout<<"CalcEnergyForElemsP8Iter["<<Id<<" is running!"<<std::endl;


	size_t IdC0 = g_treeBarrier*(Id+1);
	size_t IdCE = IdC0+g_treeBarrier;
	size_t tree = MIN(IdCE,N_CORES);
	for (size_t i=IdC0;i<tree;++i){
		FRAME(calcEnergyForElemsP8Iter[i])=CalcEnergyForElemsP8IterCD{0,1,getTP(),SHORTWAIT,i} ;
		ADD (calcEnergyForElemsP8Iter +i );
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
	
	std::cout<<"CalcEnergyForElemsSync is running!"<<std::endl;
	
	LOAD_FRAME(CalcEnergyForElemsTP);
	SIGNAL(signalUp);
	EXIT_TP();

}
	
