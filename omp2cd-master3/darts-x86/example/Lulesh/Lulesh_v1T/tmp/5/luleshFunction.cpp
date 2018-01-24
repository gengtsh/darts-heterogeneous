#include "luleshFunction.h"

void CalcForceForNodesP1(Domain& domain,Index_t lw, Index_t hi)
{
//  Index_t numNode = domain.numNode() ;

  for (Index_t i=lw; i<hi; ++i) {
     domain.fx(i) = Real_t(0.0) ;
     domain.fy(i) = Real_t(0.0) ;
     domain.fz(i) = Real_t(0.0) ;
  }

}


void IntegrateStressForElemsP1( Domain &domain,Real_t *sigxx, Real_t *sigyy, Real_t *sigzz, Real_t *determ, Real_t *fx_elem,Real_t *fy_elem, Real_t *fz_elem,Index_t lw, Index_t hi)
{

//#pragma omp parallel for firstprivate(numElem)
	for( Index_t k=lw ; k<hi ; ++k )
	{
		const Index_t* const elemToNode = domain.nodelist(k);
		Real_t B[3][8] ;// shape function derivatives
		Real_t x_local[8] ;
		Real_t y_local[8] ;
		Real_t z_local[8] ;
		
		// get nodal coordinates from global arrays and copy into local arrays.
		CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);
		
		// Volume calculation involves extra work for numerical consistency
		CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
		                                     B, &determ[k]);
		
		CalcElemNodeNormals( B[0] , B[1], B[2],
		                      x_local, y_local, z_local );
		
		// Eliminate thread writing conflicts at the nodes by giving
		// each element its own copy to write to
		SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],&fx_elem[k*8],&fy_elem[k*8],&fz_elem[k*8] ) ;
	}
		
}


void IntegrateStressForElemsP2( Domain &domain,Real_t *fx_elem,Real_t *fy_elem,Real_t *fz_elem,Index_t lw, Index_t hi)
{

     // If threaded, then we need to copy the data out of the temporary
     // arrays used above into the final forces field
//#pragma omp parallel for firstprivate(numNode)
     for( Index_t gnode=lw ; gnode<hi ; ++gnode )
     {
        Index_t count = domain.nodeElemCount(gnode) ;
        Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
        Real_t fx_tmp = Real_t(0.0) ;
        Real_t fy_tmp = Real_t(0.0) ;
        Real_t fz_tmp = Real_t(0.0) ;
        for (Index_t i=0 ; i < count ; ++i) {
           Index_t elem = cornerList[i] ;
           fx_tmp += fx_elem[elem] ;
           fy_tmp += fy_elem[elem] ;
           fz_tmp += fz_elem[elem] ;
        }
        domain.fx(gnode) = fx_tmp ;
        domain.fy(gnode) = fy_tmp ;
        domain.fz(gnode) = fz_tmp ;
     }
}

void CalcFBHourglassForceForElems_p1( Domain &domain,
                                   Real_t *determ,
                                   Real_t *x8n, Real_t *y8n, Real_t *z8n,
                                   Real_t *dvdx, Real_t *dvdy, Real_t *dvdz,
                                   Real_t hourg, 
								   Real_t *fx_elem,Real_t *fy_elem, Real_t *fz_elem, Real_t (*gamma)[8],Index_t lw, Index_t hi )
{

   /*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/
  

/*************************************************/
/*    compute the hourglass modes */


//#pragma omp parallel for firstprivate(numElem, hourg)
   for(Index_t i2=lw;i2<hi;++i2){
      Real_t *fx_local, *fy_local, *fz_local ;
      Real_t hgfx[8], hgfy[8], hgfz[8] ;

      Real_t coefficient;

      Real_t hourgam[8][4];
      Real_t xd1[8], yd1[8], zd1[8] ;

      const Index_t *elemToNode = domain.nodelist(i2);
      Index_t i3=8*i2;
      Real_t volinv=Real_t(1.0)/determ[i2];
      Real_t ss1, mass1, volume13 ;
      for(Index_t i1=0;i1<4;++i1){

         Real_t hourmodx =
            x8n[i3] * gamma[i1][0] + x8n[i3+1] * gamma[i1][1] +
            x8n[i3+2] * gamma[i1][2] + x8n[i3+3] * gamma[i1][3] +
            x8n[i3+4] * gamma[i1][4] + x8n[i3+5] * gamma[i1][5] +
            x8n[i3+6] * gamma[i1][6] + x8n[i3+7] * gamma[i1][7];

         Real_t hourmody =
            y8n[i3] * gamma[i1][0] + y8n[i3+1] * gamma[i1][1] +
            y8n[i3+2] * gamma[i1][2] + y8n[i3+3] * gamma[i1][3] +
            y8n[i3+4] * gamma[i1][4] + y8n[i3+5] * gamma[i1][5] +
            y8n[i3+6] * gamma[i1][6] + y8n[i3+7] * gamma[i1][7];

         Real_t hourmodz =
            z8n[i3] * gamma[i1][0] + z8n[i3+1] * gamma[i1][1] +
            z8n[i3+2] * gamma[i1][2] + z8n[i3+3] * gamma[i1][3] +
            z8n[i3+4] * gamma[i1][4] + z8n[i3+5] * gamma[i1][5] +
            z8n[i3+6] * gamma[i1][6] + z8n[i3+7] * gamma[i1][7];

         hourgam[0][i1] = gamma[i1][0] -  volinv*(dvdx[i3  ] * hourmodx +
                                                  dvdy[i3  ] * hourmody +
                                                  dvdz[i3  ] * hourmodz );

         hourgam[1][i1] = gamma[i1][1] -  volinv*(dvdx[i3+1] * hourmodx +
                                                  dvdy[i3+1] * hourmody +
                                                  dvdz[i3+1] * hourmodz );

         hourgam[2][i1] = gamma[i1][2] -  volinv*(dvdx[i3+2] * hourmodx +
                                                  dvdy[i3+2] * hourmody +
                                                  dvdz[i3+2] * hourmodz );

         hourgam[3][i1] = gamma[i1][3] -  volinv*(dvdx[i3+3] * hourmodx +
                                                  dvdy[i3+3] * hourmody +
                                                  dvdz[i3+3] * hourmodz );

         hourgam[4][i1] = gamma[i1][4] -  volinv*(dvdx[i3+4] * hourmodx +
                                                  dvdy[i3+4] * hourmody +
                                                  dvdz[i3+4] * hourmodz );

         hourgam[5][i1] = gamma[i1][5] -  volinv*(dvdx[i3+5] * hourmodx +
                                                  dvdy[i3+5] * hourmody +
                                                  dvdz[i3+5] * hourmodz );

         hourgam[6][i1] = gamma[i1][6] -  volinv*(dvdx[i3+6] * hourmodx +
                                                  dvdy[i3+6] * hourmody +
                                                  dvdz[i3+6] * hourmodz );

         hourgam[7][i1] = gamma[i1][7] -  volinv*(dvdx[i3+7] * hourmodx +
                                                  dvdy[i3+7] * hourmody +
                                                  dvdz[i3+7] * hourmodz );

      }

      /* compute forces */
      /* store forces into h arrays (force arrays) */

      ss1=domain.ss(i2);
      mass1=domain.elemMass(i2);
      volume13=CBRT(determ[i2]);

      Index_t n0si2 = elemToNode[0];
      Index_t n1si2 = elemToNode[1];
      Index_t n2si2 = elemToNode[2];
      Index_t n3si2 = elemToNode[3];
      Index_t n4si2 = elemToNode[4];
      Index_t n5si2 = elemToNode[5];
      Index_t n6si2 = elemToNode[6];
      Index_t n7si2 = elemToNode[7];

      xd1[0] = domain.xd(n0si2);
      xd1[1] = domain.xd(n1si2);
      xd1[2] = domain.xd(n2si2);
      xd1[3] = domain.xd(n3si2);
      xd1[4] = domain.xd(n4si2);
      xd1[5] = domain.xd(n5si2);
      xd1[6] = domain.xd(n6si2);
      xd1[7] = domain.xd(n7si2);

      yd1[0] = domain.yd(n0si2);
      yd1[1] = domain.yd(n1si2);
      yd1[2] = domain.yd(n2si2);
      yd1[3] = domain.yd(n3si2);
      yd1[4] = domain.yd(n4si2);
      yd1[5] = domain.yd(n5si2);
      yd1[6] = domain.yd(n6si2);
      yd1[7] = domain.yd(n7si2);

      zd1[0] = domain.zd(n0si2);
      zd1[1] = domain.zd(n1si2);
      zd1[2] = domain.zd(n2si2);
      zd1[3] = domain.zd(n3si2);
      zd1[4] = domain.zd(n4si2);
      zd1[5] = domain.zd(n5si2);
      zd1[6] = domain.zd(n6si2);
      zd1[7] = domain.zd(n7si2);

      coefficient = - hourg * Real_t(0.01) * ss1 * mass1 / volume13;

      CalcElemFBHourglassForce(xd1,yd1,zd1,
                      hourgam,
                      coefficient, hgfx, hgfy, hgfz);

      // With the threaded version, we write into local arrays per elem
      // so we don't have to worry about race conditions
         fx_local = &fx_elem[i3] ;
         fx_local[0] = hgfx[0];
         fx_local[1] = hgfx[1];
         fx_local[2] = hgfx[2];
         fx_local[3] = hgfx[3];
         fx_local[4] = hgfx[4];
         fx_local[5] = hgfx[5];
         fx_local[6] = hgfx[6];
         fx_local[7] = hgfx[7];

         fy_local = &fy_elem[i3] ;
         fy_local[0] = hgfy[0];
         fy_local[1] = hgfy[1];
         fy_local[2] = hgfy[2];
         fy_local[3] = hgfy[3];
         fy_local[4] = hgfy[4];
         fy_local[5] = hgfy[5];
         fy_local[6] = hgfy[6];
         fy_local[7] = hgfy[7];

         fz_local = &fz_elem[i3] ;
         fz_local[0] = hgfz[0];
         fz_local[1] = hgfz[1];
         fz_local[2] = hgfz[2];
         fz_local[3] = hgfz[3];
         fz_local[4] = hgfz[4];
         fz_local[5] = hgfz[5];
         fz_local[6] = hgfz[6];
         fz_local[7] = hgfz[7];
   }

   return ;
}


void CalcFBHourglassForceForElems_p2( Domain &domain,
								   Real_t *fx_elem,Real_t *fy_elem,Real_t *fz_elem,
								   Index_t lw,Index_t hi)
{

     // Collect the data from the local arrays into the final force arrays
//#pragma omp parallel for firstprivate(numNode)
      for( Index_t gnode=lw ; gnode<hi ; ++gnode )
      {
         Index_t count = domain.nodeElemCount(gnode) ;
         Index_t *cornerList = domain.nodeElemCornerList(gnode) ;
         Real_t fx_tmp = Real_t(0.0) ;
         Real_t fy_tmp = Real_t(0.0) ;
         Real_t fz_tmp = Real_t(0.0) ;
         for (Index_t i=0 ; i < count ; ++i) {
            Index_t elem = cornerList[i] ;
            fx_tmp += fx_elem[elem] ;
            fy_tmp += fy_elem[elem] ;
            fz_tmp += fz_elem[elem] ;
         }
         domain.fx(gnode) += fx_tmp ;
         domain.fy(gnode) += fy_tmp ;
         domain.fz(gnode) += fz_tmp ;
      }
	 
	  return;

}

void CalcHourglassControlForElems_p1(Domain& domain,Real_t determ[], Real_t *dvdx,Real_t *dvdy,Real_t *dvdz,Real_t *x8n,Real_t *y8n,Real_t *z8n ,Index_t lw, Index_t hi )
{

   /* start loop over elements */
//#pragma omp parallel for firstprivate(numElem)
   for (Index_t i=lw ; i<hi ; ++i){
      Real_t  x1[8],  y1[8],  z1[8] ;
      Real_t pfx[8], pfy[8], pfz[8] ;

      Index_t* elemToNode = domain.nodelist(i);
      CollectDomainNodesToElemNodes(domain, elemToNode, x1, y1, z1);

      CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

      /* load into temporary storage for FB Hour Glass control */
      for(Index_t ii=0;ii<8;++ii){
         Index_t jj=8*i+ii;

         dvdx[jj] = pfx[ii];
         dvdy[jj] = pfy[ii];
         dvdz[jj] = pfz[ii];
         x8n[jj]  = x1[ii];
         y8n[jj]  = y1[ii];
         z8n[jj]  = z1[ii];
      }

      determ[i] = domain.volo(i) * domain.v(i);

      /* Do a check for negative volumes */
      if ( domain.v(i) <= Real_t(0.0) ) {
//#if USE_MPI         
//         MPI_Abort(MPI_COMM_WORLD, VolumeError) ;
//#else
		  std::cout<<"CalcHourglassControlForElems_p1!"<<std::endl;
		  exit(VolumeError);
//#endif
      }
   }

   return ;
}

void ApplyAccelerationBoundaryConditionsForNodes_darts(Domain& domain,Index_t lw, Index_t hi)
{
//   Index_t size = domain.sizeX();
//   Index_t numNodeBC = (size+1)*(size+1) ;

//#pragma omp parallel
   {
      if (!domain.symmXempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=lw ; i<hi ; ++i)
            domain.xdd(domain.symmX(i)) = Real_t(0.0) ;
      }

      if (!domain.symmYempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=lw ; i<hi ; ++i)
            domain.ydd(domain.symmY(i)) = Real_t(0.0) ;
      }

      if (!domain.symmZempty() != 0) {
//#pragma omp for nowait firstprivate(numNodeBC)
         for(Index_t i=lw; i<hi ; ++i)
            domain.zdd(domain.symmZ(i)) = Real_t(0.0) ;
      }
   }
}


void CalcVelocityForNodes_darts(Domain &domain, const Real_t dt, const Real_t u_cut,Index_t lw,Index_t hi)
{

//#pragma omp parallel for firstprivate(numNode)
   for ( Index_t i = lw ; i < hi ; ++i )
   {
     Real_t xdtmp, ydtmp, zdtmp ;

     xdtmp = domain.xd(i) + domain.xdd(i) * dt ;
     if( FABS(xdtmp) < u_cut ) xdtmp = Real_t(0.0);
     domain.xd(i) = xdtmp ;

     ydtmp = domain.yd(i) + domain.ydd(i) * dt ;
     if( FABS(ydtmp) < u_cut ) ydtmp = Real_t(0.0);
     domain.yd(i) = ydtmp ;

     zdtmp = domain.zd(i) + domain.zdd(i) * dt ;
     if( FABS(zdtmp) < u_cut ) zdtmp = Real_t(0.0);
     domain.zd(i) = zdtmp ;
   }
}


void CalcPositionForNodes_darts(Domain &domain, const Real_t dt, Index_t lw, Index_t hi)
{
//#pragma omp parallel for firstprivate(numNode)
   for ( Index_t i = lw ; i < hi ; ++i )
   {
     domain.x(i) += domain.xd(i) * dt ;
     domain.y(i) += domain.yd(i) * dt ;
     domain.z(i) += domain.zd(i) * dt ;
   }
}

void CalcKinematicsForElems_darts( Domain &domain, Real_t *vnew,Real_t deltaTime, Index_t lw, Index_t hi )
{
  // loop over all elements
//#pragma omp parallel for firstprivate(numElem, deltaTime)
  for( Index_t k=lw ; k<hi ; ++k )
  {
    Real_t B[3][8] ; /** shape function derivatives */
    Real_t D[6] ;
    Real_t x_local[8] ;
    Real_t y_local[8] ;
    Real_t z_local[8] ;
    Real_t xd_local[8] ;
    Real_t yd_local[8] ;
    Real_t zd_local[8] ;
    Real_t detJ = Real_t(0.0) ;

    Real_t volume ;
    Real_t relativeVolume ;
    const Index_t* const elemToNode = domain.nodelist(k) ;

    // get nodal coordinates from global arrays and copy into local arrays.
    CollectDomainNodesToElemNodes(domain, elemToNode, x_local, y_local, z_local);

    // volume calculations
    volume = CalcElemVolume3(x_local, y_local, z_local );
    relativeVolume = volume / domain.volo(k) ;
    vnew[k] = relativeVolume ;
    domain.delv(k) = relativeVolume - domain.v(k) ;

    // set characteristic length
    domain.arealg(k) = CalcElemCharacteristicLength(x_local, y_local, z_local,
                                             volume);

    // get nodal velocities from global array and copy into local arrays.
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = elemToNode[lnode];
      xd_local[lnode] = domain.xd(gnode);
      yd_local[lnode] = domain.yd(gnode);
      zd_local[lnode] = domain.zd(gnode);
    }

    Real_t dt2 = Real_t(0.5) * deltaTime;
    for ( Index_t j=0 ; j<8 ; ++j )
    {
       x_local[j] -= dt2 * xd_local[j];
       y_local[j] -= dt2 * yd_local[j];
       z_local[j] -= dt2 * zd_local[j];
    }

    CalcElemShapeFunctionDerivatives( x_local, y_local, z_local,
                                      B, &detJ );

    CalcElemVelocityGradient( xd_local, yd_local, zd_local,
                               B, detJ, D );

    // put velocity gradient quantities into their global arrays.
    domain.dxx(k) = D[0];
    domain.dyy(k) = D[1];
    domain.dzz(k) = D[2];
  }
}


void CalcLagrangeElementsP2_darts(Domain& domain, Real_t* vnew,Index_t lw,Index_t hi)
{

      // element loop to do some stuff not included in the elemlib function.
//#pragma omp parallel for firstprivate(numElem)
      for ( Index_t k=lw ; k<hi ; ++k )
      {
         // calc strain rate and apply as constraint (only done in FB element)
         Real_t vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k) ;
         Real_t vdovthird = vdov/Real_t(3.0) ;

         // make the rate of deformation tensor deviatoric
         domain.vdov(k) = vdov ;
         domain.dxx(k) -= vdovthird ;
         domain.dyy(k) -= vdovthird ;
         domain.dzz(k) -= vdovthird ;

        // See if any volumes are negative, and take appropriate action.
         if (vnew[k] <= Real_t(0.0))
        {
			std::cout<<"CalcLagrangeElementsP2_darts!"<<std::endl;
			exit(VolumeError);
        }
      }
}

void CalcMonotonicQGradientsForElems_darts(Domain& domain, Real_t vnew[],Index_t lw,Index_t hi)
{
//#pragma omp parallel for firstprivate(numElem)
   for (Index_t i = lw ; i < hi ; ++i ) {
      const Real_t ptiny = Real_t(1.e-36) ;
      Real_t ax,ay,az ;
      Real_t dxv,dyv,dzv ;

      const Index_t *elemToNode = domain.nodelist(i);
      Index_t n0 = elemToNode[0] ;
      Index_t n1 = elemToNode[1] ;
      Index_t n2 = elemToNode[2] ;
      Index_t n3 = elemToNode[3] ;
      Index_t n4 = elemToNode[4] ;
      Index_t n5 = elemToNode[5] ;
      Index_t n6 = elemToNode[6] ;
      Index_t n7 = elemToNode[7] ;

      Real_t x0 = domain.x(n0) ;
      Real_t x1 = domain.x(n1) ;
      Real_t x2 = domain.x(n2) ;
      Real_t x3 = domain.x(n3) ;
      Real_t x4 = domain.x(n4) ;
      Real_t x5 = domain.x(n5) ;
      Real_t x6 = domain.x(n6) ;
      Real_t x7 = domain.x(n7) ;

      Real_t y0 = domain.y(n0) ;
      Real_t y1 = domain.y(n1) ;
      Real_t y2 = domain.y(n2) ;
      Real_t y3 = domain.y(n3) ;
      Real_t y4 = domain.y(n4) ;
      Real_t y5 = domain.y(n5) ;
      Real_t y6 = domain.y(n6) ;
      Real_t y7 = domain.y(n7) ;

      Real_t z0 = domain.z(n0) ;
      Real_t z1 = domain.z(n1) ;
      Real_t z2 = domain.z(n2) ;
      Real_t z3 = domain.z(n3) ;
      Real_t z4 = domain.z(n4) ;
      Real_t z5 = domain.z(n5) ;
      Real_t z6 = domain.z(n6) ;
      Real_t z7 = domain.z(n7) ;

      Real_t xv0 = domain.xd(n0) ;
      Real_t xv1 = domain.xd(n1) ;
      Real_t xv2 = domain.xd(n2) ;
      Real_t xv3 = domain.xd(n3) ;
      Real_t xv4 = domain.xd(n4) ;
      Real_t xv5 = domain.xd(n5) ;
      Real_t xv6 = domain.xd(n6) ;
      Real_t xv7 = domain.xd(n7) ;

      Real_t yv0 = domain.yd(n0) ;
      Real_t yv1 = domain.yd(n1) ;
      Real_t yv2 = domain.yd(n2) ;
      Real_t yv3 = domain.yd(n3) ;
      Real_t yv4 = domain.yd(n4) ;
      Real_t yv5 = domain.yd(n5) ;
      Real_t yv6 = domain.yd(n6) ;
      Real_t yv7 = domain.yd(n7) ;

      Real_t zv0 = domain.zd(n0) ;
      Real_t zv1 = domain.zd(n1) ;
      Real_t zv2 = domain.zd(n2) ;
      Real_t zv3 = domain.zd(n3) ;
      Real_t zv4 = domain.zd(n4) ;
      Real_t zv5 = domain.zd(n5) ;
      Real_t zv6 = domain.zd(n6) ;
      Real_t zv7 = domain.zd(n7) ;

      Real_t vol = domain.volo(i)*vnew[i] ;
      Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

      Real_t dxj = Real_t(-0.25)*((x0+x1+x5+x4) - (x3+x2+x6+x7)) ;
      Real_t dyj = Real_t(-0.25)*((y0+y1+y5+y4) - (y3+y2+y6+y7)) ;
      Real_t dzj = Real_t(-0.25)*((z0+z1+z5+z4) - (z3+z2+z6+z7)) ;

      Real_t dxi = Real_t( 0.25)*((x1+x2+x6+x5) - (x0+x3+x7+x4)) ;
      Real_t dyi = Real_t( 0.25)*((y1+y2+y6+y5) - (y0+y3+y7+y4)) ;
      Real_t dzi = Real_t( 0.25)*((z1+z2+z6+z5) - (z0+z3+z7+z4)) ;

      Real_t dxk = Real_t( 0.25)*((x4+x5+x6+x7) - (x0+x1+x2+x3)) ;
      Real_t dyk = Real_t( 0.25)*((y4+y5+y6+y7) - (y0+y1+y2+y3)) ;
      Real_t dzk = Real_t( 0.25)*((z4+z5+z6+z7) - (z0+z1+z2+z3)) ;

      /* find delvk and delxk ( i cross j ) */

      ax = dyi*dzj - dzi*dyj ;
      ay = dzi*dxj - dxi*dzj ;
      az = dxi*dyj - dyi*dxj ;

      domain.delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = Real_t(0.25)*((xv4+xv5+xv6+xv7) - (xv0+xv1+xv2+xv3)) ;
      dyv = Real_t(0.25)*((yv4+yv5+yv6+yv7) - (yv0+yv1+yv2+yv3)) ;
      dzv = Real_t(0.25)*((zv4+zv5+zv6+zv7) - (zv0+zv1+zv2+zv3)) ;

      domain.delv_zeta(i) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxi and delvi ( j cross k ) */

      ax = dyj*dzk - dzj*dyk ;
      ay = dzj*dxk - dxj*dzk ;
      az = dxj*dyk - dyj*dxk ;

      domain.delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = Real_t(0.25)*((xv1+xv2+xv6+xv5) - (xv0+xv3+xv7+xv4)) ;
      dyv = Real_t(0.25)*((yv1+yv2+yv6+yv5) - (yv0+yv3+yv7+yv4)) ;
      dzv = Real_t(0.25)*((zv1+zv2+zv6+zv5) - (zv0+zv3+zv7+zv4)) ;

      domain.delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxj and delvj ( k cross i ) */

      ax = dyk*dzi - dzk*dyi ;
      ay = dzk*dxi - dxk*dzi ;
      az = dxk*dyi - dyk*dxi ;

      domain.delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = Real_t(-0.25)*((xv0+xv1+xv5+xv4) - (xv3+xv2+xv6+xv7)) ;
      dyv = Real_t(-0.25)*((yv0+yv1+yv5+yv4) - (yv3+yv2+yv6+yv7)) ;
      dzv = Real_t(-0.25)*((zv0+zv1+zv5+zv4) - (zv3+zv2+zv6+zv7)) ;

      domain.delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
   }
}


void CalcMonotonicQRegionForElems_darts(Domain &domain, Int_t r,Real_t vnew[], Real_t ptiny,Index_t lw,Index_t hi)
{
   Real_t monoq_limiter_mult = domain.monoq_limiter_mult();
   Real_t monoq_max_slope = domain.monoq_max_slope();
   Real_t qlc_monoq = domain.qlc_monoq();
   Real_t qqc_monoq = domain.qqc_monoq();

//#pragma omp parallel for firstprivate(qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
   for ( Index_t ielem = lw ; ielem < hi; ++ielem ) {
      Index_t i = domain.regElemlist(r,ielem);
      Real_t qlin, qquad ;
      Real_t phixi, phieta, phizeta ;
      Int_t bcMask = domain.elemBC(i) ;
      Real_t delvm = 0.0, delvp =0.0;

      /*  phixi     */
      Real_t norm = Real_t(1.) / (domain.delv_xi(i)+ ptiny ) ;

      switch (bcMask & XI_M) {
         case XI_M_COMM: /* needs comm data */
         case 0:         delvm = domain.delv_xi(domain.lxim(i)); break ;
         case XI_M_SYMM: delvm = domain.delv_xi(i) ;       break ;
         case XI_M_FREE: delvm = Real_t(0.0) ;      break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & XI_P) {
         case XI_P_COMM: /* needs comm data */
         case 0:         delvp = domain.delv_xi(domain.lxip(i)) ; break ;
         case XI_P_SYMM: delvp = domain.delv_xi(i) ;       break ;
         case XI_P_FREE: delvp = Real_t(0.0) ;      break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phixi = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm < phixi ) phixi = delvm ;
      if ( delvp < phixi ) phixi = delvp ;
      if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
      if ( phixi > monoq_max_slope) phixi = monoq_max_slope;


      /*  phieta     */
      norm = Real_t(1.) / ( domain.delv_eta(i) + ptiny ) ;

      switch (bcMask & ETA_M) {
         case ETA_M_COMM: /* needs comm data */
         case 0:          delvm = domain.delv_eta(domain.letam(i)) ; break ;
         case ETA_M_SYMM: delvm = domain.delv_eta(i) ;        break ;
         case ETA_M_FREE: delvm = Real_t(0.0) ;        break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & ETA_P) {
         case ETA_P_COMM: /* needs comm data */
         case 0:          delvp = domain.delv_eta(domain.letap(i)) ; break ;
         case ETA_P_SYMM: delvp = domain.delv_eta(i) ;        break ;
         case ETA_P_FREE: delvp = Real_t(0.0) ;        break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phieta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm  < phieta ) phieta = delvm ;
      if ( delvp  < phieta ) phieta = delvp ;
      if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
      if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

      /*  phizeta     */
      norm = Real_t(1.) / ( domain.delv_zeta(i) + ptiny ) ;

      switch (bcMask & ZETA_M) {
         case ZETA_M_COMM: /* needs comm data */
         case 0:           delvm = domain.delv_zeta(domain.lzetam(i)) ; break ;
         case ZETA_M_SYMM: delvm = domain.delv_zeta(i) ;         break ;
         case ZETA_M_FREE: delvm = Real_t(0.0) ;          break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvm = 0; /* ERROR - but quiets the compiler */
            break;
      }
      switch (bcMask & ZETA_P) {
         case ZETA_P_COMM: /* needs comm data */
         case 0:           delvp = domain.delv_zeta(domain.lzetap(i)) ; break ;
         case ZETA_P_SYMM: delvp = domain.delv_zeta(i) ;         break ;
         case ZETA_P_FREE: delvp = Real_t(0.0) ;          break ;
         default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                   __FILE__, __LINE__);
            delvp = 0; /* ERROR - but quiets the compiler */
            break;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phizeta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm   < phizeta ) phizeta = delvm ;
      if ( delvp   < phizeta ) phizeta = delvp ;
      if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
      if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

      /* Remove length scale */

      if ( domain.vdov(i) > Real_t(0.) )  {
         qlin  = Real_t(0.) ;
         qquad = Real_t(0.) ;
      }
      else {
         Real_t delvxxi   = domain.delv_xi(i)   * domain.delx_xi(i)   ;
         Real_t delvxeta  = domain.delv_eta(i)  * domain.delx_eta(i)  ;
         Real_t delvxzeta = domain.delv_zeta(i) * domain.delx_zeta(i) ;

         if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
         if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
         if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

         Real_t rho = domain.elemMass(i) / (domain.volo(i) * vnew[i]) ;

         qlin = -qlc_monoq * rho *
            (  delvxxi   * (Real_t(1.) - phixi) +
               delvxeta  * (Real_t(1.) - phieta) +
               delvxzeta * (Real_t(1.) - phizeta)  ) ;

         qquad = qqc_monoq * rho *
            (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
               delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
               delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
      }

      domain.qq(i) = qquad ;
      domain.ql(i) = qlin  ;
   }
}


void ApplyMaterialPropertiesForElemsP1_darts(Domain& domain, Real_t vnew[],Index_t lw,Index_t hi)
{

    /* Expose all of the variables needed for material evaluation */
    Real_t eosvmin = domain.eosvmin() ;
    Real_t eosvmax = domain.eosvmax() ;

//#pragma omp parallel
    {
       // Bound the updated relative volumes with eosvmin/max
       if (eosvmin != Real_t(0.)) {
//#pragma omp for firstprivate(numElem)
          for(Index_t i=lw ; i<hi ; ++i) {
             if (vnew[i] < eosvmin)
                vnew[i] = eosvmin ;
          }
       }

       if (eosvmax != Real_t(0.)) {
//#pragma omp for nowait firstprivate(numElem)
          for(Index_t i=lw ; i<hi ; ++i) {
             if (vnew[i] > eosvmax)
                vnew[i] = eosvmax ;
          }
       }

       // This check may not make perfect sense in LULESH, but
       // it's representative of something in the full code -
       // just leave it in, please
//#pragma omp for nowait firstprivate(numElem)
       for (Index_t i=lw; i<hi; ++i) {
          Real_t vc = domain.v(i) ;
          if (eosvmin != Real_t(0.)) {
             if (vc < eosvmin)
                vc = eosvmin ;
          }
          if (eosvmax != Real_t(0.)) {
             if (vc > eosvmax)
                vc = eosvmax ;
          }
          if (vc <= 0.) {
			  std::cout<<"ApplyMaterialPropertiesForElemsP1_darts!"<<std::endl;
			  exit(VolumeError);
          }
       }
    }

}


void ApplyMaterialPropertiesForElemsP2_darts(Domain& domain, Real_t vnew[])
{
    for (Int_t r=0 ; r<domain.numReg() ; r++) {
       Index_t numElemReg = domain.regElemSize(r);
       Index_t *regElemList = domain.regElemlist(r);
       Int_t rep;
       //Determine load imbalance for this region
       //round down the number with lowest cost
       if(r < domain.numReg()/2)
	 rep = 1;
       //you don't get an expensive region unless you at least have 5 regions
       else if(r < (domain.numReg() - (domain.numReg()+15)/20))
         rep = 1 + domain.cost();
       //very expensive regions
       else
	 rep = 10 * (1+ domain.cost());
       EvalEOSForElems(domain, vnew, numElemReg, regElemList, rep);
    }

}


void EvalEOSForElemsP1_darts(Domain& domain, Real_t *vnewc,Real_t* p_new,Real_t *e_new, Real_t *q_new, Real_t *pbvc, Real_t *bvc,Int_t numElemReg, Index_t *regElemList, Int_t rep)
{
   Real_t  e_cut = domain.e_cut() ;
   Real_t  p_cut = domain.p_cut() ;
//   Real_t  ss4o3 = domain.ss4o3() ;
   Real_t  q_cut = domain.q_cut() ;

   Real_t eosvmax = domain.eosvmax() ;
   Real_t eosvmin = domain.eosvmin() ;
   Real_t pmin    = domain.pmin() ;
   Real_t emin    = domain.emin() ;
   Real_t rho0    = domain.refdens() ;

   // These temporaries will be of different size for 
   // each call (due to different sized region element
   // lists)
   Real_t *e_old = Allocate<Real_t>(numElemReg) ;
   Real_t *delvc = Allocate<Real_t>(numElemReg) ;
   Real_t *p_old = Allocate<Real_t>(numElemReg) ;
   Real_t *q_old = Allocate<Real_t>(numElemReg) ;
   Real_t *compression = Allocate<Real_t>(numElemReg) ;
   Real_t *compHalfStep = Allocate<Real_t>(numElemReg) ;
   Real_t *qq_old = Allocate<Real_t>(numElemReg) ;
   Real_t *ql_old = Allocate<Real_t>(numElemReg) ;
   Real_t *work = Allocate<Real_t>(numElemReg) ;
 
   //loop to add load imbalance based on region number 
   for(Int_t j = 0; j < rep; j++) {
      /* compress data, minimal set */
//#pragma omp parallel
      {
//#pragma omp for nowait firstprivate(numElemReg)
         for (Index_t i=0; i<numElemReg; ++i) {
            Index_t elem = regElemList[i];
            e_old[i] = domain.e(elem) ;
            delvc[i] = domain.delv(elem) ;
            p_old[i] = domain.p(elem) ;
            q_old[i] = domain.q(elem) ;
            qq_old[i] = domain.qq(elem) ;
            ql_old[i] = domain.ql(elem) ;
         }

//#pragma omp for firstprivate(numElemReg)
         for (Index_t i = 0; i < numElemReg ; ++i) {
            Index_t elem = regElemList[i];
            Real_t vchalf ;
            compression[i] = Real_t(1.) / vnewc[elem] - Real_t(1.);
            vchalf = vnewc[elem] - delvc[i] * Real_t(.5);
            compHalfStep[i] = Real_t(1.) / vchalf - Real_t(1.);
         }

      /* Check for v > eosvmax or v < eosvmin */
         if ( eosvmin != Real_t(0.) ) {
//#pragma omp for nowait firstprivate(numElemReg, eosvmin)
            for(Index_t i=0 ; i<numElemReg ; ++i) {
               Index_t elem = regElemList[i];
               if (vnewc[elem] <= eosvmin) { /* impossible due to calling func? */
                  compHalfStep[i] = compression[i] ;
               }
            }
         }
         if ( eosvmax != Real_t(0.) ) {
//#pragma omp for nowait firstprivate(numElemReg, eosvmax)
            for(Index_t i=0 ; i<numElemReg ; ++i) {
               Index_t elem = regElemList[i];
               if (vnewc[elem] >= eosvmax) { /* impossible due to calling func? */
                  p_old[i]        = Real_t(0.) ;
                  compression[i]  = Real_t(0.) ;
                  compHalfStep[i] = Real_t(0.) ;
               }
            }
         }

//#pragma omp for nowait firstprivate(numElemReg)
         for (Index_t i = 0 ; i < numElemReg ; ++i) {
            work[i] = Real_t(0.) ; 
         }
      }
      CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                         p_old, e_old,  q_old, compression, compHalfStep,
                         vnewc, work,  delvc, pmin,
                         p_cut, e_cut, q_cut, emin,
                         qq_old, ql_old, rho0, eosvmax,
                         numElemReg, regElemList);
   }

   Release(&work) ;
   Release(&ql_old) ;
   Release(&qq_old) ;
   Release(&compHalfStep) ;
   Release(&compression) ;
   Release(&q_old) ;
   Release(&p_old) ;
   Release(&delvc) ;
   Release(&e_old) ;
}


void EvalEOSForElemsP1TPP1_darts(Domain &domain,Real_t *vnewc,Index_t lw,Index_t hi,Index_t *regElemList,Real_t *e_old,Real_t *delvc,Real_t *p_old,Real_t *q_old,Real_t *compression,Real_t *compHalfStep,Real_t *qq_old,Real_t *ql_old ){


//#pragma omp for nowait firstprivate(numElemReg)
	for (Index_t i=lw; i<hi; ++i) {
		Index_t elem = regElemList[i];
		e_old[i] = domain.e(elem) ;
		delvc[i] = domain.delv(elem) ;
		p_old[i] = domain.p(elem) ;
		q_old[i] = domain.q(elem) ;
		qq_old[i] = domain.qq(elem) ;
		ql_old[i] = domain.ql(elem) ;
	}

//#pragma omp for firstprivate(numElemReg)
	for (Index_t i = lw; i < hi ; ++i) {
		Index_t elem = regElemList[i];
		Real_t vchalf ;
		compression[i] = Real_t(1.) / vnewc[elem] - Real_t(1.);
		vchalf = vnewc[elem] - delvc[i] * Real_t(.5);
		compHalfStep[i] = Real_t(1.) / vchalf - Real_t(1.);
	}

}

void EvalEOSForElemsP1TPP2_darts(Real_t *vnewc, Index_t lw,Index_t hi,Index_t *regElemList,Real_t eosvmin,Real_t eosvmax,Real_t *p_old,Real_t *compression,Real_t *compHalfStep,Real_t *work){

/* Check for v > eosvmax or v < eosvmin */
	if ( eosvmin != Real_t(0.) ) {
//#pragma omp for nowait firstprivate(numElemReg, eosvmin)
		for(Index_t i=lw ; i<hi ; ++i) {
			Index_t elem = regElemList[i];
			if (vnewc[elem] <= eosvmin) { /* impossible due to calling func? */
				compHalfStep[i] = compression[i] ;
			}
		}
    }
	if ( eosvmax != Real_t(0.) ) {
//#pragma omp for nowait firstprivate(numElemReg, eosvmax)
		for(Index_t i=lw ; i<hi ; ++i) {
			Index_t elem = regElemList[i];
			if (vnewc[elem] >= eosvmax) { /* impossible due to calling func? */
				p_old[i]        = Real_t(0.) ;
				compression[i]  = Real_t(0.) ;
				compHalfStep[i] = Real_t(0.) ;
			}
		}
	}

//#pragma omp for nowait firstprivate(numElemReg)
	for (Index_t i = lw ; i < hi ; ++i) {
		work[i] = Real_t(0.) ; 
	}

}

void EvalEOSForElemsP2_darts(Domain& domain, Index_t *regElemList,Real_t *p_new,Real_t *e_new,Real_t *q_new, Index_t lw,Index_t hi)
{

   for (Index_t i=lw; i<hi; ++i) {
      Index_t elem = regElemList[i];
      domain.p(elem) = p_new[i] ;
      domain.e(elem) = e_new[i] ;
      domain.q(elem) = q_new[i] ;
   }
}

void CalcEnergyForElemsP1_darts(Real_t *e_new, Index_t lw,Index_t hi,Real_t *e_old,Real_t *delvc,Real_t *p_old,Real_t *q_old,Real_t *work,Real_t emin){

	//#pragma omp parallel for firstprivate(length, emin)
	for (Index_t i = lw ; i < hi ; ++i) {
		e_new[i] = e_old[i] - Real_t(0.5) * delvc[i] * (p_old[i] + q_old[i])+ Real_t(0.5) * work[i];
	
		if (e_new[i]  < emin ) {
			e_new[i] = emin ;
		}
	}
}

void CalcEnergyForElemsP3_darts(Real_t *e_new,Real_t *q_new,Index_t lw,Index_t hi,Real_t *delvc,Real_t *p_old,Real_t *q_old,Real_t *pHalfStep,Real_t * compHalfStep,Real_t *bvc,Real_t *pbvc,Real_t *ql_old,Real_t *qq_old,Real_t rho0){

//#pragma omp parallel for firstprivate(length, rho0)
	for (Index_t i = lw ; i < hi ; ++i) {
		Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep[i]) ;
		if ( delvc[i] > Real_t(0.) ) {
			q_new[i] /* = qq_old[i] = ql_old[i] */ = Real_t(0.) ;
		}
		else {
			Real_t ssc = ( pbvc[i] * e_new[i]+ vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0 ;
	
			if ( ssc <= Real_t(.1111111e-36) ) {
				ssc = Real_t(.3333333e-18) ;
			} else {
				ssc = SQRT(ssc) ;
			}
	
			q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;
		}
	
		e_new[i] = e_new[i] + Real_t(0.5) * delvc[i]
				* (  Real_t(3.0)*(p_old[i]     + q_old[i])
				- Real_t(4.0)*(pHalfStep[i] + q_new[i])) ;
	}

}


void CalcEnergyForElemsP4_darts(Real_t *e_new,Index_t lw,Index_t hi,Real_t *work, Real_t e_cut,Real_t emin){

//#pragma omp parallel for firstprivate(length, emin, e_cut)
	for (Index_t i = lw ; i < hi ; ++i) {
	
		e_new[i] += Real_t(0.5) * work[i];
		
		if (FABS(e_new[i]) < e_cut) {
			e_new[i] = Real_t(0.)  ;
		}
		if (     e_new[i]  < emin ) {
			e_new[i] = emin ;
		}
	}

}

void CalcEnergyForElemsP6_darts(Real_t *e_new,Real_t *p_new,Real_t *q_new,Real_t *vnewc,Index_t *regElemList,Index_t lw,Index_t hi,Real_t *delvc,Real_t *p_old,Real_t *q_old,Real_t *pHalfStep, Real_t *bvc,Real_t *pbvc,Real_t *ql_old,Real_t *qq_old,Real_t e_cut,Real_t emin,Real_t rho0){

//#pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
   for (Index_t i = lw ; i < hi ; ++i){
      const Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
      Index_t elem = regElemList[i];
      Real_t q_tilde ;

      if (delvc[i] > Real_t(0.)) {
         q_tilde = Real_t(0.) ;
      }
      else {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= Real_t(.1111111e-36) ) {
            ssc = Real_t(.3333333e-18) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_tilde = (ssc*ql_old[i] + qq_old[i]) ;
      }

      e_new[i] = e_new[i] - (  Real_t(7.0)*(p_old[i]     + q_old[i])
                               - Real_t(8.0)*(pHalfStep[i] + q_new[i])
                               + (p_new[i] + q_tilde)) * delvc[i]*sixth ;

      if (FABS(e_new[i]) < e_cut) {
         e_new[i] = Real_t(0.)  ;
      }
      if (     e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }
}

void CalcEnergyForElemsP8_darts(Real_t *q_new,Real_t *e_new,Real_t *p_new,Real_t *vnewc,Index_t *regElemList,Index_t lw,Index_t hi,Real_t *bvc,Real_t *pbvc,Real_t *delvc,Real_t *ql_old,Real_t *qq_old,Real_t q_cut,Real_t rho0){

//#pragma omp parallel for firstprivate(length, rho0, q_cut)
   for (Index_t i = lw ; i < hi ; ++i){
      Index_t elem = regElemList[i];

      if ( delvc[i] <= Real_t(0.) ) {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[elem] * vnewc[elem] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= Real_t(.1111111e-36) ) {
            ssc = Real_t(.3333333e-18) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_new[i] = (ssc*ql_old[i] + qq_old[i]) ;

         if (FABS(q_new[i]) < q_cut) q_new[i] = Real_t(0.) ;
      }
   }
}


void CalcPressureForElemsP1_darts(Real_t* bvc,
                          Real_t* pbvc, 
                          Real_t* compression,
                          Index_t lw,Index_t hi)
{
//#pragma omp parallel for firstprivate(length)
   for (Index_t i =lw; i < hi ; ++i) {
      Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
      bvc[i] = c1s * (compression[i] + Real_t(1.));
      pbvc[i] = c1s;
   }

}

void CalcPressureForElemsP2_darts(Real_t* p_new, Real_t* bvc,
                          Real_t* e_old,
                          Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t lw,Index_t hi, Index_t *regElemList)
{

//#pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
   for (Index_t i = lw ; i < hi ; ++i){
      Index_t elem = regElemList[i];
      
      p_new[i] = bvc[i] * e_old[i] ;

      if    (FABS(p_new[i]) <  p_cut   )
         p_new[i] = Real_t(0.0) ;

      if    ( vnewc[elem] >= eosvmax ) /* impossible condition here? */
         p_new[i] = Real_t(0.0) ;

      if    (p_new[i]       <  pmin)
         p_new[i]   = pmin ;
   }
}

void CalcSoundSpeedForElems_darts(Domain &domain,
                            Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, 
                            Index_t lw,Index_t hi, Index_t *regElemList)
{
//#pragma omp parallel for firstprivate(rho0, ss4o3)
   for (Index_t i = lw; i < hi ; ++i) {
      Index_t elem = regElemList[i];
      Real_t ssTmp = (pbvc[i] * enewc[i] + vnewc[elem] * vnewc[elem] *
                 bvc[i] * pnewc[i]) / rho0;
      if (ssTmp <= Real_t(.1111111e-36)) {
         ssTmp = Real_t(.3333333e-18);
      }
      else {
         ssTmp = SQRT(ssTmp);
      }
      domain.ss(elem) = ssTmp ;
   }
}


void UpdateVolumesForElems_darts(Domain &domain, Real_t *vnew,Real_t v_cut,Index_t lw,Index_t hi)
{
//#pragma omp parallel for firstprivate(length, v_cut)
      for(Index_t i=lw ; i<hi ; ++i) {
         Real_t tmpV = vnew[i] ;

         if ( FABS(tmpV - Real_t(1.0)) < v_cut )
            tmpV = Real_t(1.0) ;

         domain.v(i) = tmpV ;
      }

   return ;
}



void CalcCourantConstraintForElemsP1_darts(Domain &domain,Index_t *regElemlist,Real_t qqc, Real_t& dtcourant,Index_t *courant_elem_per_thread,Real_t *dtcourant_per_thread,Index_t Id,Index_t lw,Index_t hi)
{
//#pragma omp parallel firstprivate(length, qqc)
   {
      Real_t   qqc2 = Real_t(64.0) * qqc * qqc ;
      Real_t   dtcourant_tmp = dtcourant;
      Index_t  courant_elem  = -1 ;

      Index_t thread_num = Id;

//#pragma omp for 
      for (Index_t i = lw ; i < hi ; ++i) {
         Index_t indx = regElemlist[i] ;
         Real_t dtf = domain.ss(indx) * domain.ss(indx) ;

         if ( domain.vdov(indx) < Real_t(0.) ) {
            dtf = dtf
                + qqc2 * domain.arealg(indx) * domain.arealg(indx)
                * domain.vdov(indx) * domain.vdov(indx) ;
         }

         dtf = SQRT(dtf) ;
         dtf = domain.arealg(indx) / dtf ;

         if (domain.vdov(indx) != Real_t(0.)) {
            if ( dtf < dtcourant_tmp ) {
               dtcourant_tmp = dtf ;
               courant_elem  = indx ;
            }
         }
      }

      dtcourant_per_thread[thread_num]    = dtcourant_tmp ;
      courant_elem_per_thread[thread_num] = courant_elem ;
   }

   return ;

}

void CalcCourantConstraintForElemsP2_darts(Real_t& dtcourant,Index_t *courant_elem_per_thread,Real_t *dtcourant_per_thread,Index_t numCores)
{

   for (Index_t i = 1; i < numCores; ++i) {
      if (dtcourant_per_thread[i] < dtcourant_per_thread[0] ) {
         dtcourant_per_thread[0]    = dtcourant_per_thread[i];
         courant_elem_per_thread[0] = courant_elem_per_thread[i];
      }
   }

   if (courant_elem_per_thread[0] != -1) {
      dtcourant = dtcourant_per_thread[0] ;
   }

   return ;

}


void CalcHydroConstraintForElemsP1_darts(Domain &domain,Index_t *regElemlist, Real_t dvovmax, Real_t& dthydro,Index_t *hydro_elem_per_thread, Real_t *dthydro_per_thread,Index_t Id,Index_t lw,Index_t hi)
{

//#pragma omp parallel firstprivate(length, dvovmax)
   {
      Real_t dthydro_tmp = dthydro ;
      Index_t hydro_elem = -1 ;

      Index_t thread_num = Id;

//#pragma omp for
      for (Index_t i = lw ; i < hi ; ++i) {
         Index_t indx = regElemlist[i] ;

         if (domain.vdov(indx) != Real_t(0.)) {
            Real_t dtdvov = dvovmax / (FABS(domain.vdov(indx))+Real_t(1.e-20)) ;

            if ( dthydro_tmp > dtdvov ) {
                  dthydro_tmp = dtdvov ;
                  hydro_elem = indx ;
            }
         }
      }

      dthydro_per_thread[thread_num]    = dthydro_tmp ;
      hydro_elem_per_thread[thread_num] = hydro_elem ;
   }

   if (hydro_elem_per_thread[0] != -1) {
      dthydro =  dthydro_per_thread[0] ;
   }

   return ;
}

void CalcHydroConstraintForElemsP2_darts(Real_t& dthydro,Index_t *hydro_elem_per_thread, Real_t *dthydro_per_thread,Index_t numCores)
{
   for (Index_t i = 1; i < numCores; ++i) {
      if(dthydro_per_thread[i] < dthydro_per_thread[0]) {
         dthydro_per_thread[0]    = dthydro_per_thread[i];
         hydro_elem_per_thread[0] =  hydro_elem_per_thread[i];
      }
   }

   if (hydro_elem_per_thread[0] != -1) {
      dthydro =  dthydro_per_thread[0] ;
   }

   return ;
}
