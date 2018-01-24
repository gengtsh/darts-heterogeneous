#ifndef DARTS_CALCTIMECONSTRAINTSFORELEMSTP_H
#define DARTS_CALCTIMECONSTRAINTSFORELEMSTP_H 

#include <stdint.h>
#include "luleshMain.h"
#include "lulesh.h"
#include "luleshFunction.h"

using namespace darts;


DEF_CODELET(CalcTimeConstraintsForElemsCD,1,SHORTWAIT);
DEF_CODELET(CalcTimeConstraintsForElemsSyncCD,1,SHORTWAIT);

DEF_TP(CalcTimeConstraintsForElemsTP)
{
	Domain *domain;

	Index_t *courant_elem_per_thread;
	Real_t *dtcourant_per_thread;
	Index_t *hydro_elem_per_thread;
	Real_t *dthydro_per_thread;
	
	CalcTimeConstraintsForElemsCD calcTimeConstraintsForElems;
	CalcTimeConstraintsForElemsSyncCD calcTimeConstraintsForElemsSync;
	Codelet *signalUp;
	CalcTimeConstraintsForElemsTP(Domain *domain,Codelet *up)
		:domain(domain)
		,calcTimeConstraintsForElems(0,1,this,SHORTWAIT)
		,signalUp(up)
	{
		courant_elem_per_thread = new Index_t[N_CORES];
		dtcourant_per_thread = new Real_t[N_CORES];
		
		hydro_elem_per_thread = new Index_t[N_CORES];
		dthydro_per_thread = new Real_t[N_CORES];
		
		domain->dtcourant() = 1.0e+20;
		domain->dthydro() = 1.0e+20;
		
		uint32_t numReg=domain->numReg();	
		calcTimeConstraintsForElemsSync =CalcTimeConstraintsForElemsSyncCD{numReg*2,numReg*2,this,SHORTWAIT}; 	
		add(&calcTimeConstraintsForElems);
		add(&calcTimeConstraintsForElemsSync);
	}
	virtual ~CalcTimeConstraintsForElemsTP(){
	
		delete []courant_elem_per_thread ;
		delete []dtcourant_per_thread ;
	
		delete []hydro_elem_per_thread ;
		delete []dthydro_per_thread ;
	}
};



#endif
