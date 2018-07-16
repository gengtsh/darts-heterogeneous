
#ifndef SPENCILCPUKERNEL_H
#define SPENCILCPUKERNEL_H

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstdio>

void computeInner_stencil2d_v2(double *dst,double *src,size_t rpos2,size_t n_cols);

#endif
