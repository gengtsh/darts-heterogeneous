
#include <stdint.h>
#include <stdlib.h>
#include "StencilCPUKernel.h"


void computeInner_stencil2d_v2(double *dst,double *src,size_t rpos2,size_t n_cols){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;

	for (size_t i = 0; i < rpos2; ++i) {
	    for (size_t j = 0; j < n_cols-2; ++j) {
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
	    }
	}
}

