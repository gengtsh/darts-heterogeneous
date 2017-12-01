#ifndef DARTS_CALCFBHOURGLASSFORCEFORELEMSTP_H
#define DARTS_CALCFBHOURGLASSFORCEFORELEMSTP_H

#include <stdint.h>
#include "lulesh.h"
#include "luleshFunction.h"
#include "luleshMain.h"


DEF_CODELET_ITER(CalcFBHourglassForceForElemsP1IterCD,1,SHORTWAIT);
DEF_CODELET_ITER(CalcFBHourglassForceForElemsP2IterCD,1,SHORTWAIT);

DEF_CODELET(CalcFBHourglassForceForElemsP2CD,1,LONGWAIT);
DEF_CODELET(CalcFBHourglassForceForElemsSyncCD,1,LONGWAIT);

DEF_TP(CalcFBHourglassForceForElemsTP)
{
	
	Domain *locDom;
	Real_t *determ;


	Real_t *x8n  ;
	Real_t *y8n  ;
	Real_t *z8n  ;
	Real_t *dvdx ;
	Real_t *dvdy ;
	Real_t *dvdz ;
	
	Real_t hgcoef;
	
	Index_t numElem;
	Index_t numNode;	
	Index_t numElem8;
	
	Real_t *fx_elem; 
	Real_t *fy_elem; 
	Real_t *fz_elem; 

	Real_t  gamma[4][8];


	CalcFBHourglassForceForElemsP1IterCD *calcFBHourglassForceForElemsP1Iter;
	CalcFBHourglassForceForElemsP2CD	  calcFBHourglassForceForElemsP2;
	CalcFBHourglassForceForElemsP2IterCD *calcFBHourglassForceForElemsP2Iter;
	CalcFBHourglassForceForElemsSyncCD calcFBHourglassForceForElemsSync;
	Codelet *signalUp;
		
	CalcFBHourglassForceForElemsTP(Domain *locDom,Real_t *determ, Real_t *x8n, Real_t *y8n, Real_t *z8n, Real_t *dvdx, Real_t *dvdy, Real_t *dvdz, Real_t hgcoef,Index_t numElem, Index_t numNode,Codelet *up ) 
	:locDom(locDom)
	,determ(determ)
	,x8n(x8n)
	,y8n(y8n)
	,z8n(z8n)
	,dvdx(dvdx)
	,dvdy(dvdy)
	,dvdz(dvdz)
	,hgcoef(hgcoef)
	,numElem(numElem)
	,numNode(numNode)
	,calcFBHourglassForceForElemsP1Iter(new CalcFBHourglassForceForElemsP1IterCD[N_CORES])
	,calcFBHourglassForceForElemsP2(N_CORES,N_CORES,this,SHORTWAIT)
	,calcFBHourglassForceForElemsP2Iter(new CalcFBHourglassForceForElemsP2IterCD[N_CORES])
	,calcFBHourglassForceForElemsSync(N_CORES,N_CORES,this,LONGWAIT)
	,signalUp(up)
	{
		numElem8	= numElem * 8;

		fx_elem = Allocate<Real_t>(numElem8) ;
		fy_elem = Allocate<Real_t>(numElem8) ;
		fz_elem = Allocate<Real_t>(numElem8) ;


		gamma[0][0] = Real_t( 1.);
		gamma[0][1] = Real_t( 1.);
		gamma[0][2] = Real_t(-1.);
		gamma[0][3] = Real_t(-1.);
		gamma[0][4] = Real_t(-1.);
		gamma[0][5] = Real_t(-1.);
		gamma[0][6] = Real_t( 1.);
		gamma[0][7] = Real_t( 1.);
		gamma[1][0] = Real_t( 1.);
		gamma[1][1] = Real_t(-1.);
		gamma[1][2] = Real_t(-1.);
		gamma[1][3] = Real_t( 1.);
		gamma[1][4] = Real_t(-1.);
		gamma[1][5] = Real_t( 1.);
		gamma[1][6] = Real_t( 1.);
		gamma[1][7] = Real_t(-1.);
		gamma[2][0] = Real_t( 1.);
		gamma[2][1] = Real_t(-1.);
		gamma[2][2] = Real_t( 1.);
		gamma[2][3] = Real_t(-1.);
		gamma[2][4] = Real_t( 1.);
		gamma[2][5] = Real_t(-1.);
		gamma[2][6] = Real_t( 1.);
		gamma[2][7] = Real_t(-1.);
		gamma[3][0] = Real_t(-1.);
		gamma[3][1] = Real_t( 1.);
		gamma[3][2] = Real_t(-1.);
		gamma[3][3] = Real_t( 1.);
		gamma[3][4] = Real_t( 1.);
		gamma[3][5] = Real_t(-1.);
		gamma[3][6] = Real_t( 1.);
		gamma[3][7] = Real_t(-1.);


////		for ( size_t i = 0; i < N_CORES; ++i ) {
////            calcFBHourglassForceForElemsP1Iter[i] = CalcFBHourglassForceForElemsP1IterCD {0,1,this,SHORTWAIT,i};
////            calcFBHourglassForceForElemsP2Iter[i] = CalcFBHourglassForceForElemsP2IterCD {1,1,this,SHORTWAIT,i};
////			add( calcFBHourglassForceForElemsP1Iter + i);
////			add( calcFBHourglassForceForElemsP2Iter + i);
////        }
////		add (&calcFBHourglassForceForElemsP2);
////		add (&calcFBHourglassForceForElemsSync);
	

////		for ( size_t i = 0; i < N_CORES; ++i ) {
////            calcFBHourglassForceForElemsP1Iter[i] = CalcFBHourglassForceForElemsP1IterCD {0,1,this,SHORTWAIT,i};
////            calcFBHourglassForceForElemsP2Iter[i] = CalcFBHourglassForceForElemsP2IterCD {1,1,this,SHORTWAIT,i};
////			add( calcFBHourglassForceForElemsP1Iter + i);
////        }

        calcFBHourglassForceForElemsP1Iter[0] = CalcFBHourglassForceForElemsP1IterCD {0,1,this,SHORTWAIT,0};
		add( calcFBHourglassForceForElemsP1Iter + 0);
		if(N_CORES>1){
			calcFBHourglassForceForElemsP1Iter[1] = CalcFBHourglassForceForElemsP1IterCD {0,1,this,SHORTWAIT,1};
			add( calcFBHourglassForceForElemsP1Iter + 1);
		}


	}

	~CalcFBHourglassForceForElemsTP(){
		
		delete []calcFBHourglassForceForElemsP1Iter;
		delete []calcFBHourglassForceForElemsP2Iter;

		Release(&fz_elem) ;
		Release(&fy_elem) ;
		Release(&fx_elem) ;
	}
};

#endif
