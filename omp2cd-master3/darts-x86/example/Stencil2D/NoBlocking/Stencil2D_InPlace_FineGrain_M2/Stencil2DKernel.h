#ifndef SPENCIL2D_KERNEL_H
#define SPENCIL2D_KERNEL_H

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include "DARTS.h"

/**
* Naive 4pt stencil for 2D array
*/
void stencil2d_seq(double *dst,double *src,const uint64_t n_rows,const uint64_t n_cols,uint64_t n_tsteps);

uint32_t computeRowDecomposition(const uint64_t n_rows, const uint64_t n_cols);


void computeInner_stencil2d(uint64_t BlockPosition,double *InitialMatrix,double *NewMatrix, uint64_t BlockM, uint64_t BlockN,const uint64_t InitialN);

void computeInner_stencil2d_v2(double *dst,double *src,size_t rpos2,size_t n_cols);

static inline void swap_ptr(void** left, void** right) {
	    void* tmp = *left;
		    *left     = *right;
		    *right    = tmp;
}

#define SWAP_PTR(left,right) swap_ptr((void**)left,(void**)right)


#endif
