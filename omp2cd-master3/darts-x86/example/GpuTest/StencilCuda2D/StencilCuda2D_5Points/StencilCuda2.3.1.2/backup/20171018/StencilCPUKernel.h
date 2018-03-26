
#ifndef SPENCILCPUKERNEL_H
#define SPENCILCPUKERNEL_H

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstdio>

void computeInner_stencil2d_v2(double *dst,double *src,size_t n_rows,size_t n_cols);

void copyColumns_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride );

void copyRows_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride);


#endif
