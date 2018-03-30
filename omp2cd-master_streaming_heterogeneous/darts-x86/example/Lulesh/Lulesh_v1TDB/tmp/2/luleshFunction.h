#ifndef LULESHFUNCTION_H
#define LULESHFUNCTION_H


#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <unistd.h>

#if _OPENMP
# include <omp.h>
#endif

#include "lulesh.h"


/*********************************/
/* Data structure implementation */
/*********************************/

/* might want to add access methods so that memory can be */
/* better managed, as in luleshFT */

template <typename T>
T *Allocate(size_t size)
{
   return static_cast<T *>(malloc(sizeof(T)*size)) ;
}

template <typename T>
void Release(T **ptr)
{
   if (*ptr != NULL) {
      free(*ptr) ;
      *ptr = NULL ;
   }
}



/******************************************/

/* Work Routines */

//static inline
void TimeIncrement(Domain& domain);
/******************************************/

//static inline
void CollectDomainNodesToElemNodes(Domain &domain,
                                   const Index_t* elemToNode,
                                   Real_t elemX[8],
                                   Real_t elemY[8],
                                   Real_t elemZ[8]);
/******************************************/

//static inline
void InitStressTermsForElems(Domain &domain,
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                             Index_t numElem);
/******************************************/

//static inline
void CalcElemShapeFunctionDerivatives( Real_t const x[],
                                       Real_t const y[],
                                       Real_t const z[],
                                       Real_t b[][8],
                                       Real_t* const volume );
/******************************************/

//static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3);
/******************************************/

//static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8]);
/******************************************/

//static inline
void SumElemStressesToNodeForces( const Real_t B[][8],
                                  const Real_t stress_xx,
                                  const Real_t stress_yy,
                                  const Real_t stress_zz,
                                  Real_t fx[], Real_t fy[], Real_t fz[] );
/******************************************/

//static inline
void IntegrateStressForElems( Domain &domain,
                              Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                              Real_t *determ, Index_t numElem, Index_t numNode);


void IntegrateStressForElemsP1( Domain &domain,Real_t *sigxx, Real_t *sigyy, Real_t *sigzz, Real_t *determ, Real_t *fx_elem,Real_t *fy_elem, Real_t *fz_elem,Index_t lw,Index_t hi);
void IntegrateStressForElemsP2( Domain &domain,Real_t *fx_elem,Real_t *fy_elem,Real_t *fz_elem, Index_t lw, Index_t hi);

/******************************************/

//static inline
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t* dvdx, Real_t* dvdy, Real_t* dvdz);
/******************************************/

//static inline
void CalcElemVolumeDerivative(Real_t dvdx[8],
                              Real_t dvdy[8],
                              Real_t dvdz[8],
                              const Real_t x[8],
                              const Real_t y[8],
                              const Real_t z[8]);
/******************************************/

//static inline
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t hourgam[][4],
                              Real_t coefficient,
                              Real_t *hgfx, Real_t *hgfy, Real_t *hgfz );

/******************************************/

//static inline
void CalcFBHourglassForceForElems( Domain &domain,
                                   Real_t *determ,
                                   Real_t *x8n, Real_t *y8n, Real_t *z8n,
                                   Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
                                   Real_t hourg, Index_t numElem,
                                   Index_t numNode);


void CalcFBHourglassForceForElems_p1( Domain &domain,
                                   Real_t *determ,
                                   Real_t *x8n, Real_t *y8n, Real_t *z8n,
                                   Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
                                   Real_t hourg, 
								   Real_t *fx_elem,Real_t *fy_elem, Real_t *fz_elem, Real_t (*gamma)[8],Index_t lw, Index_t hi );


void CalcFBHourglassForceForElems_p2( Domain &domain,Real_t *fx_elem,Real_t *fy_elem,Real_t* fz_elem,Index_t lw,Index_t hi);
/******************************************/

//static inline
void CalcHourglassControlForElems(Domain& domain,
                                  Real_t determ[], Real_t hgcoef);

void CalcHourglassControlForElems_p1(Domain& domain,Real_t determ[], Real_t *dvdx,Real_t *dvdy,Real_t *dvdz,Real_t *x8n,Real_t *y8n,Real_t *z8n ,Index_t lw, Index_t hi );

/******************************************/

//static inline
void CalcVolumeForceForElems(Domain& domain);

/******************************************/

//static inline 
void CalcForceForNodes(Domain& domain);
void CalcForceForNodesP1(Domain& domain,Index_t lw, Index_t hi);
/******************************************/

//static inline
void CalcAccelerationForNodes(Domain &domain, Index_t numNode);
/******************************************/

//static inline
void ApplyAccelerationBoundaryConditionsForNodes(Domain& domain);

void ApplyAccelerationBoundaryConditionsForNodes_darts(Domain& domain,Index_t lw,Index_t hi);

/******************************************/

//static inline
void CalcVelocityForNodes(Domain &domain, const Real_t dt, const Real_t u_cut,Index_t numNode);
void CalcVelocityForNodes_darts(Domain &domain, const Real_t dt, const Real_t u_cut,Index_t lw,Index_t hi);

/******************************************/

//static inline
void CalcPositionForNodes(Domain &domain, const Real_t dt, Index_t numNode);

void CalcPositionForNodes_darts(Domain &domain, const Real_t dt, Index_t lw, Index_t hi);
/******************************************/

//static inline
void LagrangeNodal(Domain& domain);
/******************************************/

//static inline
Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
               const Real_t x2, const Real_t x3,
               const Real_t x4, const Real_t x5,
               const Real_t x6, const Real_t x7,
               const Real_t y0, const Real_t y1,
               const Real_t y2, const Real_t y3,
               const Real_t y4, const Real_t y5,
               const Real_t y6, const Real_t y7,
               const Real_t z0, const Real_t z1,
               const Real_t z2, const Real_t z3,
               const Real_t z4, const Real_t z5,
               const Real_t z6, const Real_t z7 );
/******************************************/

//inline
Real_t CalcElemVolume3( const Real_t x[8], const Real_t y[8], const Real_t z[8] );

/******************************************/

//static inline
Real_t AreaFace( const Real_t x0, const Real_t x1,
                 const Real_t x2, const Real_t x3,
                 const Real_t y0, const Real_t y1,
                 const Real_t y2, const Real_t y3,
                 const Real_t z0, const Real_t z1,
                 const Real_t z2, const Real_t z3);

/******************************************/

//static inline
Real_t CalcElemCharacteristicLength( const Real_t x[8],
                                     const Real_t y[8],
                                     const Real_t z[8],
                                     const Real_t volume);

/******************************************/

//static inline
void CalcElemVelocityGradient( const Real_t* const xvel,
                                const Real_t* const yvel,
                                const Real_t* const zvel,
                                const Real_t b[][8],
                                const Real_t detJ,
                                Real_t* const d );

/******************************************/

////static inline
void CalcKinematicsForElems( Domain &domain, Real_t *vnew, 
                             Real_t deltaTime, Index_t numElem );

void CalcKinematicsForElems_darts( Domain &domain, Real_t *vnew,Real_t deltaTime, Index_t lw, Index_t hi );

/******************************************/

//static inline
void CalcLagrangeElements(Domain& domain, Real_t* vnew);

void CalcLagrangeElementsP2_darts(Domain& domain, Real_t* vnew,Index_t lw,Index_t hi);
/******************************************/

//static inline
void CalcMonotonicQGradientsForElems(Domain& domain, Real_t vnew[]);
/******************************************/

//static inline
void CalcMonotonicQRegionForElems(Domain &domain, Int_t r,
                                  Real_t vnew[], Real_t ptiny);
/******************************************/

//static inline
void CalcMonotonicQForElems(Domain& domain, Real_t vnew[]);
/******************************************/

//static inline
void CalcQForElems(Domain& domain, Real_t vnew[]);
/******************************************/

//static inline
void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
                          Real_t* pbvc, Real_t* e_old,
                          Real_t* compression, Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t length, Index_t *regElemList);
/******************************************/

//static inline
void CalcEnergyForElems(Real_t* p_new, Real_t* e_new, Real_t* q_new,
                        Real_t* bvc, Real_t* pbvc,
                        Real_t* p_old, Real_t* e_old, Real_t* q_old,
                        Real_t* compression, Real_t* compHalfStep,
                        Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,
                        Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
                        Real_t* qq_old, Real_t* ql_old,
                        Real_t rho0,
                        Real_t eosvmax,
                        Index_t length, Index_t *regElemList);
/******************************************/

//static inline
//void CalcSoundSpeedForElems(Domain &domain,
//                            Real_t *vnewc, Real_t rho0, Real_t *enewc,
//                            Real_t *pnewc, Real_t *pbvc,
//                            Real_t *bvc, Real_t ss4o3,
//                            Index_t len, Index_t *regElemList);

void CalcSoundSpeedForElems(Domain &domain,
                            Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, 
                            Index_t len, Index_t *regElemList);


/******************************************/

//static inline
void EvalEOSForElems(Domain& domain, Real_t *vnewc,
                     Int_t numElemReg, Index_t *regElemList, Int_t rep);

/******************************************/

//static inline
void ApplyMaterialPropertiesForElems(Domain& domain, Real_t vnew[]);

/******************************************/

//static inline
void UpdateVolumesForElems(Domain &domain, Real_t *vnew,
                           Real_t v_cut, Index_t length);

/******************************************/

//static inline
void LagrangeElements(Domain& domain, Index_t numElem);

/******************************************/

//static inline
void CalcCourantConstraintForElems(Domain &domain, Index_t length,
                                   Index_t *regElemlist,
                                   Real_t qqc, Real_t& dtcourant);

/******************************************/

//static inline
void CalcHydroConstraintForElems(Domain &domain, Index_t length,
                                 Index_t *regElemlist, Real_t dvovmax, Real_t& dthydro);


//static inline
void CalcTimeConstraintsForElems(Domain& domain);


/******************************************/

//static inline
void LagrangeLeapFrog(Domain& domain);

#endif
