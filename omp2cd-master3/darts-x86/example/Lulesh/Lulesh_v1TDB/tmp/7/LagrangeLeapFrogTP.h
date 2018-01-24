#ifndef LAGRANGELEAPFROG_H
#define LAGRANGELEAPFROG_H

#include <stdint.h>
#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

using namespace darts;



DEF_CODELET(LagrangeNodalCD,0,SHORTWAIT);
DEF_CODELET(LagrangeElementsCD,1,SHORTWAIT);
DEF_CODELET(LagrangeCalcTimeConstraintsForElemsCD,1,SHORTWAIT);

DEF_CODELET(LeapFrogSyncCD,1,LONGWAIT);

DEF_TP(LagrangeLeapFrogTP)
{
	Domain *domain;
		
	LagrangeNodalCD lagrangeNodal;
	LagrangeElementsCD lagrangeElements;
	LagrangeCalcTimeConstraintsForElemsCD lagrangeCalcTimeConstraintsForElems;
	LeapFrogSyncCD sync;
	Codelet *signalUp;
	LagrangeLeapFrogTP(Domain *domain,Codelet *up)
		:domain(domain)
		,lagrangeNodal(0,1,this,SHORTWAIT)
		,lagrangeElements(1,1,this,SHORTWAIT)
		,lagrangeCalcTimeConstraintsForElems(1,1,this,SHORTWAIT)
		,sync(1,1,this, SHORTWAIT)
		,signalUp(up)
	{
		add(&lagrangeNodal);
		add(&lagrangeElements);
		add(&lagrangeCalcTimeConstraintsForElems);
		add(&sync);

	}
	virtual ~LagrangeLeapFrogTP(){}
};

#endif
