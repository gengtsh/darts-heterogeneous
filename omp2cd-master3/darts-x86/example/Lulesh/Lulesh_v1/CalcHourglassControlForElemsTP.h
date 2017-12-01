#ifndef DARTS_CalcHourglassControlForElemsTP_H
#define DARTS_CalcHourglassControlForElemsTP_H

#include <stdint.h>
#include "lulesh.h"
#include "luleshFunction.h"
#include "luleshMain.h"

DEF_CODELET_ITER(CalcHourglassControlForElemsIterCD,0,SHORTWAIT);
DEF_CODELET(CalcHourglassControlForElemsSyncCD,2,LONGWAIT);

DEF_TP(CalcHourglassControlForElemsTP)
{
	
	Domain *locDom;
	Real_t *determ;
	Real_t hgcoef;

	Index_t numElem;
	Index_t numElem8;

	Real_t *dvdx ;
	Real_t *dvdy ;
	Real_t *dvdz ;
	Real_t *x8n  ;
	Real_t *y8n  ;
	Real_t *z8n  ;

	CalcHourglassControlForElemsIterCD *calcHourglassControlForElemsIter;
	CalcHourglassControlForElemsSyncCD calcHourglassControlForElemsSync;	

//	CalcHourglassControlForElemsSyncFinalCD calcHourglassControlForElemsSyncFinal;
	Codelet *signalUp;
	
	CalcHourglassControlForElemsTP(Domain *locDom,Real_t *determ,Real_t hgcoef,Codelet *up)
	:locDom(locDom)
	,determ(determ)
	,hgcoef(hgcoef)
	,calcHourglassControlForElemsIter(new CalcHourglassControlForElemsIterCD[N_CORES])
    ,calcHourglassControlForElemsSync(N_CORES,N_CORES,this,SHORTWAIT)	
	,signalUp(up)
	{
		numElem		= locDom->numElem();
		numElem8	= numElem * 8;
		dvdx = Allocate<Real_t>(numElem8) ;
		dvdy = Allocate<Real_t>(numElem8) ;
		dvdz = Allocate<Real_t>(numElem8) ;
		x8n  = Allocate<Real_t>(numElem8) ;
		y8n  = Allocate<Real_t>(numElem8) ;
		z8n  = Allocate<Real_t>(numElem8) ;

		if (hgcoef>0){

////			for ( size_t i = 0; i < N_CORES; ++i ) {
////        	    calcHourglassControlForElemsIter[i] = CalcHourglassControlForElemsIterCD {0,1,this,SHORTWAIT,i};
////        	    add( calcHourglassControlForElemsIter + i );
////        	}
////			add (&calcHourglassControlForElemsSync);

			calcHourglassControlForElemsIter[0] = CalcHourglassControlForElemsIterCD {0,1,this,SHORTWAIT,0};
			add( calcHourglassControlForElemsIter + 0 );
			if(N_CORES>1){
				calcHourglassControlForElemsIter[1] = CalcHourglassControlForElemsIterCD {0,1,this,SHORTWAIT,1};
				add( calcHourglassControlForElemsIter + 1 );
			}
		}else {
			signalUp->decDep();
		}
	}

	~CalcHourglassControlForElemsTP(){
		
		delete []calcHourglassControlForElemsIter;

		Release(&z8n) ;
		Release(&y8n) ;
		Release(&x8n) ;
		Release(&dvdz) ;
		Release(&dvdy) ;
		Release(&dvdx) ;

	}
};

#endif
