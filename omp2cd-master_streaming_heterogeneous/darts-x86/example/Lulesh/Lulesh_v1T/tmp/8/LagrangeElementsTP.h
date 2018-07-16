#ifndef DARTS_LAGRANGEELEMENTSTP_H
#define DARTS_LAGRANGEELEMENTSTP_H

#include <stdint.h>
#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

using namespace darts;

DEF_CODELET(CalcLagrangeElementsCD,0,SHORTWAIT);
DEF_CODELET(CalcQForElemsCD,0,SHORTWAIT);
DEF_CODELET(ApplyMaterialPropertiesForElemsCD,0,SHORTWAIT);
DEF_CODELET(UpdateVolumesForElemsCD,0,SHORTWAIT);
DEF_CODELET(LagrangeElementsSyncCD,1,SHORTWAIT);

DEF_TP(LagrangeElementsTP)
{
	Domain *domain;
	Index_t numElem;
	Real_t *vnew;
	CalcLagrangeElementsCD calcLagrangeElements;
	CalcQForElemsCD calcQForElems;
	ApplyMaterialPropertiesForElemsCD applyMaterialPropertiesForElems;

	UpdateVolumesForElemsCD updateVolumesForElems;
	LagrangeElementsSyncCD lagrangeElementsSync;
	Codelet *signalUp;
	LagrangeElementsTP(Domain *domain,Index_t numElem,Codelet *up)
		:domain(domain)
		,numElem(numElem)
		,calcLagrangeElements(0,1,this,SHORTWAIT)
		,calcQForElems(N_CORES,N_CORES,this,SHORTWAIT)
		,applyMaterialPropertiesForElems(1,1,this,SHORTWAIT)
		,updateVolumesForElems(1,1,this,SHORTWAIT)
		,lagrangeElementsSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
	{
		vnew = Allocate<Real_t>(numElem) ;  /* new relative vol -- temp */
		add(&calcLagrangeElements);	
//		add(&calcQForElems);
//		add(&applyMaterialPropertiesForElems);
//		add(&updateVolumesForElems);
//		add(&lagrangeElementsSync);

	}
	virtual ~LagrangeElementsTP(){
		Release(&vnew);
	}
};



#endif
