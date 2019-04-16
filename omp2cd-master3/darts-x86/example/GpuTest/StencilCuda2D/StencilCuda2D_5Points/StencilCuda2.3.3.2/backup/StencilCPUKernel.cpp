
#include <stdint.h>
#include <stdlib.h>
#include "StencilCPUKernel.h"
#include <iostream>
//#include <math.h>
#include <cmath>
//#define DEBUG_COPY

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

void copyColumns_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride){
	int n_cols_cp = std::ceil(1.0*n_cols/stride);

	typedef double (*Array2D)[n_cols];
	
	Array2D SRC = (Array2D) src;
	
	for (size_t i = 0;i<n_rows;++i){
		dst[i] = src[i*stride];
#ifdef DEBUG_COPY
		std::cout<<"shared_cols["<<i<<"] = "<< dst[i]<<std::endl;
#endif
	}
}


void copyRows_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride){
	int cpRow;
	int sz = std::ceil(1.0*n_rows/stride);
	for (int i = 0;i<sz;++i){
		cpRow = i*n_cols*stride;
		std::copy(src+cpRow, src+n_cols*2+cpRow, dst+i*n_cols*2);
	}

#ifdef DEBUG_COPY

		std::cout<<"shared_Rows:"<<std::endl;
	for(int i=0;i<sz;++i){	
		for(int j =0;j<2*n_cols;++j){
			std::cout<<dst[2*n_cols*i+j]<<",";
		}
		std::cout<<"\n";
	}
#endif

}


void computeBlock_stencil25(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_rows_ck,size_t n_cols_ck){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;

	for (size_t i = 1; i < n_rows_ck; ++i) {
        for (size_t j = 1; j < n_cols_ck; ++j) {
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
        }
    }
}
