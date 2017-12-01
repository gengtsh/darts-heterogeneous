#ifndef DARTS_CALCPRESSUREFORELEMSTP_H
#define DARTS_CALCPRESSUREFORELEMSTP_H

#include <stdint.h>

#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

DEF_CODELET_ITER(CalcPressureForElemsP1IterCD,0,SHORTWAIT);
DEF_CODELET_ITER(CalcPressureForElemsP2IterCD,0,SHORTWAIT);
DEF_CODELET(CalcPressureForElemsSyncCD,1,SHORTWAIT );

DEF_TP(CalcPressureForElemsTP)
{
	Real_t  *p_new		 ;
	Real_t  *bvc		 ;
	Real_t  *pbvc		 ;
	Real_t  *e_old		 ;
	Real_t  *compression ;
	Real_t  *vnewc		 ;
	Real_t   pmin		 ;
	Real_t   p_cut		 ;
	Real_t   eosvmax	 ;
	Index_t  numElemReg	 ;
	Index_t *regElemList ; 
	
	CalcPressureForElemsP1IterCD *calcPressureForElemsP1Iter;
	CalcPressureForElemsP2IterCD *calcPressureForElemsP2Iter;
	CalcPressureForElemsSyncCD calcPressureForElemsSync;
	Codelet *signalUp;
	
	CalcPressureForElemsTP(Real_t* p_new, Real_t* bvc,Real_t* pbvc, Real_t* e_old,Real_t* compression, Real_t *vnewc,Real_t pmin,Real_t p_cut, Real_t eosvmax,Index_t numElemReg, Index_t *regElemList,Codelet *up)
		:p_new		 (p_new		  )		  
		,bvc		 (bvc		  )  
		,pbvc		 (pbvc		  )  
		,e_old		 (e_old		  )  
		,compression (compression )  
		,vnewc		 (vnewc		  )  
		,pmin		 (pmin		  )  
		,p_cut		 (p_cut		  )  
		,eosvmax	 (eosvmax	  )  
		,numElemReg	 (numElemReg  )  
		,regElemList (regElemList )  
		,calcPressureForElemsP1Iter(new  CalcPressureForElemsP1IterCD[N_CORES])
		,calcPressureForElemsP2Iter(new  CalcPressureForElemsP2IterCD[N_CORES])
		,calcPressureForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
		{
			for ( size_t i = 0; i < N_CORES; ++i ) {
				calcPressureForElemsP1Iter[i]=CalcPressureForElemsP1IterCD{0,1,this,SHORTWAIT,i};
				calcPressureForElemsP2Iter[i]=CalcPressureForElemsP2IterCD{1,1,this,SHORTWAIT,i};
				add ( calcPressureForElemsP1Iter+ i);
				add ( calcPressureForElemsP2Iter+ i);
			}
			add (&calcPressureForElemsSync);
		}
	virtual ~CalcPressureForElemsTP(){
	//	std::cout<<"~CalcPressureForElemsTP--"<<std::endl;	
		delete [] calcPressureForElemsP1Iter;
		delete [] calcPressureForElemsP2Iter;
	//	std::cout<<"~CalcPressureForElemsTP!!"<<std::endl;	
	}

};


#endif
