
#include <stdint.h>
#include <stdlib.h>
#include "StencilCPUKernel.h"
#include <iostream>

void computeInner_stencil2d_v2(double *dst,double *src,size_t n_rows,size_t n_cols){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;

	for (size_t i = 1; i < n_rows; ++i) {
	    for (size_t j = 1; j < n_cols-1; ++j) {
	       // DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
		}
	}
}

