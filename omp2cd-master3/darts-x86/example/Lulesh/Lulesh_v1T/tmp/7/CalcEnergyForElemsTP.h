#ifndef DARTS_CALCENERGYFORELEMSTP_H
#define DARTS_CALCENERGYFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"


DEF_CODELET_ITER(CalcEnergyForElemsP1IterCD,0,SHORTWAIT);
DEF_CODELET(CalcEnergyForElemsP2CD,1,SHORTWAIT );
DEF_CODELET(CalcEnergyForElemsP2SyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(CalcEnergyForElemsP3IterCD,1,SHORTWAIT);
DEF_CODELET_ITER(CalcEnergyForElemsP4IterCD,1,SHORTWAIT);
DEF_CODELET(CalcEnergyForElemsP5CD,1,SHORTWAIT );
DEF_CODELET(CalcEnergyForElemsP5SyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(CalcEnergyForElemsP6IterCD,1,SHORTWAIT);
DEF_CODELET(CalcEnergyForElemsP7CD,1,SHORTWAIT );
DEF_CODELET(CalcEnergyForElemsP7SyncCD,1,SHORTWAIT );
DEF_CODELET_ITER(CalcEnergyForElemsP8IterCD,1,SHORTWAIT);
DEF_CODELET(CalcEnergyForElemsSyncCD,1,SHORTWAIT );


DEF_TP( CalcEnergyForElemsTP)
{
    
	Real_t  *p_new		 ;
	Real_t  *e_new		 ;
	Real_t  *q_new		 ;
	Real_t  *bvc		 ;
	Real_t  *pbvc		 ;
	Real_t  *p_old		 ;
	Real_t  *e_old		 ;
	Real_t  *q_old		 ;
	Real_t  *compression ;
	Real_t  *compHalfStep;
	Real_t  *vnewc		 ;
	Real_t  *work		 ;
	Real_t  *delvc		 ;
	Real_t   pmin		 ;
	Real_t   p_cut		 ;
	Real_t   e_cut		 ;
	// Rea l_t  ss4o3	 ;
	Real_t   q_cut		 ;
	Real_t   emin		 ;
	Real_t  *qq_old		 ;
	Real_t  *ql_old		 ;
	Real_t   rho0		 ;
	Real_t   eosvmax	 ;
	Index_t  numElemReg	 ;
	Index_t *regElemList ; 
	Real_t  *pHalfStep	 ;


	CalcEnergyForElemsP1IterCD  *calcEnergyForElemsP1Iter  	;
	CalcEnergyForElemsP2CD		 calcEnergyForElemsP2		;
	CalcEnergyForElemsP2SyncCD	 calcEnergyForElemsP2Sync	;
	CalcEnergyForElemsP3IterCD  *calcEnergyForElemsP3Iter  	;
	CalcEnergyForElemsP4IterCD  *calcEnergyForElemsP4Iter  	;
	CalcEnergyForElemsP5CD		 calcEnergyForElemsP5	  	;
	CalcEnergyForElemsP5SyncCD	 calcEnergyForElemsP5Sync  	;
	CalcEnergyForElemsP6IterCD  *calcEnergyForElemsP6Iter  	;
	CalcEnergyForElemsP7CD		 calcEnergyForElemsP7	  	;
	CalcEnergyForElemsP7SyncCD	 calcEnergyForElemsP7Sync  	;
	CalcEnergyForElemsP8IterCD  *calcEnergyForElemsP8Iter  	;
	CalcEnergyForElemsSyncCD	 calcEnergyForElemsSync		;
	
	Codelet *signalUp;

	CalcEnergyForElemsTP(Real_t* p_new, Real_t* e_new, Real_t* q_new,Real_t* bvc, Real_t* pbvc,Real_t* p_old, Real_t* e_old, Real_t* q_old,Real_t* compression, Real_t* compHalfStep,Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,Real_t* qq_old, Real_t* ql_old,Real_t rho0,Real_t eosvmax,Index_t numElemReg, Index_t *regElemList, Codelet *up)
		:p_new		 (p_new		  )		  
		,e_new		 (e_new		  )  
		,q_new		 (q_new		  )  
		,bvc		 (bvc		  )  
		,pbvc		 (pbvc		  )  
		,p_old		 (p_old		  )  
		,e_old		 (e_old		  )  
		,q_old		 (q_old		  )  
		,compression (compression )  
		,compHalfStep(compHalfStep) 
		,vnewc		 (vnewc		  )  
		,work		 (work		  )  
		,delvc		 (delvc		  )  
		,pmin		 (pmin		  )  
		,p_cut		 (p_cut		  )  
		,e_cut		 (e_cut		  )  
		,q_cut		 (q_cut		  )  
		,emin		 (emin		  )
		,qq_old		 (qq_old	  )  
		,ql_old		 (ql_old	  )  
		,rho0		 (rho0		  )  
		,eosvmax	 (eosvmax	  )  
		,numElemReg	 (numElemReg  )  
		,regElemList (regElemList )  
		,calcEnergyForElemsP1Iter (new CalcEnergyForElemsP1IterCD[N_CORES]) 	
		,calcEnergyForElemsP2	  (N_CORES,N_CORES,this,SHORTWAIT)  
		,calcEnergyForElemsP2Sync (1,1,this,SHORTWAIT)  
		 ,calcEnergyForElemsP3Iter (new CalcEnergyForElemsP3IterCD[N_CORES]) 	
		,calcEnergyForElemsP4Iter (new CalcEnergyForElemsP4IterCD[N_CORES]) 	
		,calcEnergyForElemsP5	  (N_CORES,N_CORES,this,SHORTWAIT)	
		,calcEnergyForElemsP5Sync (1,1,this,SHORTWAIT)	
		,calcEnergyForElemsP6Iter (new CalcEnergyForElemsP6IterCD[N_CORES]) 	
		,calcEnergyForElemsP7	  (N_CORES,N_CORES,this,SHORTWAIT)	
		,calcEnergyForElemsP7Sync (1,1,this,SHORTWAIT)	
		,calcEnergyForElemsP8Iter (new CalcEnergyForElemsP8IterCD[N_CORES]) 	
		,calcEnergyForElemsSync	  (N_CORES,N_CORES,this,SHORTWAIT)	
		,signalUp(up)
		{
			pHalfStep = Allocate<Real_t>(numElemReg); 
			for (size_t i=0;i<N_CORES;++i){
				calcEnergyForElemsP1Iter[i]=CalcEnergyForElemsP1IterCD{0,1,this,SHORTWAIT,i} ;
				calcEnergyForElemsP3Iter[i]=CalcEnergyForElemsP3IterCD{1,1,this,SHORTWAIT,i} ;
				calcEnergyForElemsP4Iter[i]=CalcEnergyForElemsP4IterCD{1,1,this,SHORTWAIT,i} ;
				calcEnergyForElemsP6Iter[i]=CalcEnergyForElemsP6IterCD{1,1,this,SHORTWAIT,i} ;
				calcEnergyForElemsP8Iter[i]=CalcEnergyForElemsP8IterCD{1,1,this,SHORTWAIT,i} ;
				
				add (calcEnergyForElemsP1Iter +i );
				add (calcEnergyForElemsP3Iter +i );
				add (calcEnergyForElemsP4Iter +i );
				add (calcEnergyForElemsP6Iter +i );
				add (calcEnergyForElemsP8Iter +i );

				
			}
			add (&calcEnergyForElemsP2);
			add (&calcEnergyForElemsP5);
			add (&calcEnergyForElemsP7);
			add (&calcEnergyForElemsP2Sync);
			add (&calcEnergyForElemsP5Sync);
			add (&calcEnergyForElemsP7Sync);
			add (&calcEnergyForElemsSync);
		}

	virtual ~CalcEnergyForElemsTP(){
//		std::cout<<"~CalcEnergyForElemsTP--"<<std::endl;
		Release(&pHalfStep) ;
		delete []calcEnergyForElemsP1Iter;
		delete []calcEnergyForElemsP3Iter;
		delete []calcEnergyForElemsP4Iter;
		delete []calcEnergyForElemsP6Iter;
		delete []calcEnergyForElemsP8Iter;

//		std::cout<<"~CalcEnergyForElemsTP!!"<<std::endl;
	}

};


#endif
