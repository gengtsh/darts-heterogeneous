#include "LagrangeElementsTP.h"
#include "CalcLagrangeElementsTP.h"
#include "CalcQForElemsTP.h"
#include "ApplyMaterialPropertiesForElemsTP.h"
#include "UpdateVolumesForElemsTP.h"

void CalcLagrangeElementsCD::fire(void)
{
	
	std::cout<<"CalcLagrangeElements begin!"<<std::endl;	
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	
	INVOKE(CalcLagrangeElementsTP, domain,vnew,&FRAME(calcQForElems)) ;

	EXIT_TP();
}

void CalcQForElemsCD::fire(void)
{
	std::cout<<"CalcQForElems begin!"<<std::endl;	
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	
	INVOKE(CalcQForElemsTP, domain,vnew,&FRAME(applyMaterialPropertiesForElems )) ;

	EXIT_TP();
}

void ApplyMaterialPropertiesForElemsCD::fire(void)
{

	std::cout<<"ApplyMaterialPropertiesForElems begin!"<<std::endl;	
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	
	INVOKE(ApplyMaterialPropertiesForElemsTP, domain,vnew,&FRAME(updateVolumesForElems )) ;

	EXIT_TP();
}


void UpdateVolumesForElemsCD::fire(void)
{
	std::cout<<"UpdateVolumesForElems begin!"<<std::endl;
	LOAD_FRAME(LagrangeElementsTP);
	Domain *domain=FRAME(domain);
	Real_t *vnew = FRAME(vnew) ; 
	Index_t numElem=domain->numElem();
	
	INVOKE( UpdateVolumesForElemsTP, domain,vnew,domain->v_cut(), numElem, &FRAME(lagrangeElementsSync)) ;

	EXIT_TP();
}


void LagrangeElementsSyncCD::fire(void)
{
	std::cout<<"LagrangeElemsSync begin!"<<std::endl;
	LOAD_FRAME(LagrangeElementsTP);
//	Domain *domain=FRAME(domain);
//	Index_t numElem=domain->numElem();
//	Real_t *vnew = FRAME(vnew) ; 
	
	//CalcLagrangeElements(*domain, vnew) ;
	
	/* Calculate Q.  (Monotonic q option requires communication) */
//	CalcQForElems(*domain, vnew) ;
	
//	ApplyMaterialPropertiesForElems(*domain, vnew) ;
	
//	UpdateVolumesForElems(*domain, vnew,domain->v_cut(), numElem) ;
	
	SIGNAL(signalUp);	
	EXIT_TP();
}
