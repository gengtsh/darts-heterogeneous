#include "CalcEnergyForElemsTP.h"
#include "CalcPressureForElemsTP.h"

void CalcEnergyForElemsP1IterCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	
	Index_t numElemReg	= FRAME(numElemReg);
	                                         
	Real_t *e_new		=FRAME( e_new		)	;
	Real_t *e_old		=FRAME( e_old		)	;
	Real_t *delvc		=FRAME( delvc		)	;
	Real_t *p_old		=FRAME( p_old		)	;
	Real_t *q_old		=FRAME( q_old		)	;
	Real_t *work		=FRAME( work		)	;
	Real_t  emin		=FRAME( emin		)	;

	size_t	Chunk = numElemReg/N_CORES;
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP1_darts(e_new, lw,hi,e_old,delvc,p_old,q_old,work,emin);
	
	SYNC(calcEnergyForElemsP2 );	
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

	EXIT_TP();
}

void CalcEnergyForElemsP2SyncCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);

	for(size_t i=0;i<N_CORES;++i){
		SYNC(calcEnergyForElemsP3Iter[i]  );
	}
	EXIT_TP();
}

void CalcEnergyForElemsP3IterCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	
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
	size_t	Id	= getID();
	Index_t	lw;
	Index_t	hi;
	lw	= Id*Chunk;
	hi	= ((Id==(N_CORES-1))? (numElemReg%N_CORES):0) + (Id+1)* Chunk;

	CalcEnergyForElemsP3_darts(e_new,q_new, lw,hi,delvc,p_old,q_old,pHalfStep, compHalfStep,bvc,pbvc,ql_old,qq_old,rho0);
	
	SYNC(calcEnergyForElemsP4Iter[Id] );	
	EXIT_TP();
}


void CalcEnergyForElemsP4IterCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	
	Index_t numElemReg	= FRAME(numElemReg);

	Real_t  *e_new		  =  FRAME(e_new		);
	Real_t  *work		  =  FRAME(work			);
	Real_t   e_cut		  =  FRAME(e_cut		);
	Real_t   emin		  =  FRAME(emin			);
	
	size_t	Chunk = numElemReg/N_CORES;
	size_t	Id	= getID();
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

	EXIT_TP();
}

void CalcEnergyForElemsP5SyncCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);

	for(size_t i=0;i<N_CORES;++i){
		SYNC(calcEnergyForElemsP6Iter[i]  );
	}
	EXIT_TP();
}

void CalcEnergyForElemsP6IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	
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
	size_t	Id	= getID();
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

	
	EXIT_TP();
}

void CalcEnergyForElemsP7SyncCD ::fire(void)
{

	LOAD_FRAME(CalcEnergyForElemsTP);
	
	for(size_t i=0;i<N_CORES;++i){
		SYNC(calcEnergyForElemsP8Iter[i]  );
	}
	EXIT_TP();
}

void CalcEnergyForElemsP8IterCD ::fire(void)
{
	LOAD_FRAME(CalcEnergyForElemsTP);
	
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
	size_t	Id	= getID();
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
	SIGNAL(signalUp);
	EXIT_TP();

}
	
