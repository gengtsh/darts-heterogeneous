#ifndef LAGRANGELEAPFROG_H
#define LAGRANGELEAPFROG_H

#include <stdint.h>
#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

using namespace darts;


DEF_CODELET(LeapFrogSyncCD,1,LONGWAIT);
DEF_CODELET(NodalCD,0,SHORTWAIT);
DEF_CODELET(ElementsCD,1,SHORTWAIT);
DEF_CODELET(TimeConstraintsCD,1,SHORTWAIT);

DEF_TP(LagrangeLeapFrogTP)
{
	Domain *domain;
	NodalCD nodal;
	ElementsCD elements;
	TimeConstraintsCD timeConstraints;
	LeapFrogSyncCD sync;
	Codelet *signalUp;
	LagrangeLeapFrogTP(Domain *domain,Codelet *up)
		:domain(domain)
		,nodal(0,0,this,SHORTWAIT)
		,elements(1,1,this,SHORTWAIT)
		,timeConstraints(1,1,this,SHORTWAIT)
		,sync(1,1,this, SHORTWAIT)
		,signalUp(up)
	{
	
		add(&nodal);
		add(&elements);
		add(&timeConstraints);
		add(&sync);

	}
	virtual ~LagrangeLeapFrogTP(){}
};

#endif
