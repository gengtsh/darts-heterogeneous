

#include <stdint.h>
#include <stdlib.h>
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

/**
* Naive 4pt stencil for 2D array
*/
void stencil2d_seq(double *dst,double *src,const uint64_t n_rows,const uint64_t n_cols,uint64_t n_tsteps){
    
	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
		for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
			}
		}
		SWAP_PTR(&DST,&SRC);
	}

}

uint32_t computeRowDecomposition(const uint64_t n_rows, const uint64_t n_cols){

	uint64_t nTotalSize = (n_rows-2)*(n_cols-2);//the total number of pixel in inner matrix (the whole matrix - left_edge - right_edge - upper_edge - down_edge)
	uint64_t nRowsCut = nTotalSize/TOTAL_TILE_SZ ;	
	uint64_t nThreads = (g_nCU+1)*g_nSU*N_THREADS;//(nCU+1)*nSU is available cores, N_THREADS are number of threads per cores

	uint32_t nRowsCut_final = (nRowsCut<1)?1:nThreads;
//	uint32_t nRowsCut_final = (nRowsCut<1)?1:((nRowsCut<nThreads)?nRowsCut:nThreads);
	return nRowsCut_final;
}



void computeInner_stencil2d(uint64_t BlockPosition,double* InitialMatrix,double *NewMatrix, uint64_t BlockM, uint64_t BlockN,const uint64_t InitialN){

	uint64_t i,j;
	uint64_t i_begin = BlockPosition/InitialN;
	uint64_t j_begin = BlockPosition%InitialN;
	uint64_t i_end = i_begin + BlockM-1;
	uint64_t j_end = j_begin + BlockN-1;
	
	for(i=i_begin;i<=i_end;i++)
		for(j=j_begin;j<=j_end;j++)
			NewMatrix[i*InitialN +j]=(InitialMatrix[(i-1)*InitialN +j]+InitialMatrix[(i+1)*InitialN +j]+InitialMatrix[i*InitialN +j-1]+InitialMatrix[i*InitialN +j+1])/4;	
	
	return;
}


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
