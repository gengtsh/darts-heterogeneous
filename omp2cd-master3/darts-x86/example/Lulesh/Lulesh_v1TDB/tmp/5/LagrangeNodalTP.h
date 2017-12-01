#ifndef DARTS_LAGRANGENONAL_H
#define DARTS_LAGRANGENODAL_H

#include <stdint.h>
#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

using namespace darts;


DEF_CODELET(CalcForceForNodesCD,0,SHORTWAIT);
DEF_CODELET(CalcAccelerationForNodesCD,1,SHORTWAIT);
DEF_CODELET(ApplyAccelerationBoundaryConditionsForNodesCD,1,SHORTWAIT);
DEF_CODELET(CalcVelocityForNodesCD,1,SHORTWAIT);
DEF_CODELET(CalcPositionForNodesCD,1,SHORTWAIT);

DEF_CODELET(LagrangeNodalSyncCD,0,SHORTWAIT);

DEF_TP(LagrangeNodalTP)
{
	Domain *domain;
	CalcForceForNodesCD calcForceForNodes;
	CalcAccelerationForNodesCD calcAccelerationForNodes;
	ApplyAccelerationBoundaryConditionsForNodesCD applyAccelerationBoundaryConditionsForNodes;
 	CalcVelocityForNodesCD calcVelocityForNodes;
	CalcPositionForNodesCD calcPositionForNodes;
	LagrangeNodalSyncCD lagrangeNodalSync;
	Codelet *signalUp;
	LagrangeNodalTP(Domain *domain,Codelet *up)
		:domain(domain)
		,calcForceForNodes(0,0,this,SHORTWAIT)
		,calcAccelerationForNodes(1,1,this,SHORTWAIT)
		,applyAccelerationBoundaryConditionsForNodes(N_CORES,N_CORES,this,SHORTWAIT)
		,calcVelocityForNodes(N_CORES,N_CORES,this,SHORTWAIT)
		,calcPositionForNodes(N_CORES,N_CORES,this,SHORTWAIT)
		,lagrangeNodalSync(N_CORES,N_CORES,this,SHORTWAIT)
		,signalUp(up)
	{
		add(&calcForceForNodes);
		add(&calcAccelerationForNodes);
		add(&applyAccelerationBoundaryConditionsForNodes);
		add(&calcVelocityForNodes);
		add(&calcPositionForNodes);
		add(&lagrangeNodalSync);

	}
	virtual ~LagrangeNodalTP(){}
};

#endif
