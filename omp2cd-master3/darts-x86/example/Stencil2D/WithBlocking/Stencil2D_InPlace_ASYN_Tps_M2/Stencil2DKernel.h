#ifndef SPENCIL2D_KERNEL_H
#define SPENCIL2D_KERNEL_H

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include "DARTS.h"
#include "common.h"

//#include "Stencil2DPartition.h"

/**
* Naive 4pt stencil for 2D array
*/
void stencil2d_seq(double *dst,double *src,const uint64_t n_rows,const uint64_t n_cols,uint64_t n_tsteps);

uint32_t computeRowDecomposition(const uint64_t n_rows, const uint64_t n_cols);


void computeInner_stencil2d(uint64_t BlockPosition,double *InitialMatrix,double *NewMatrix, uint64_t BlockM, uint64_t BlockN,const uint64_t InitialN);


void computeBlock_stencil25(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_rows_ck,size_t n_cols_ck);


void computeInner_stencil2d(uint64_t BlockPosition,double *Initial,double *share,uint64_t BlockM,uint64_t BlockN,const uint64_t InitialM, const uint64_t InitialN);
void copyLine_stencil2d(uint64_t BlockPosition,double *Initial,double *Copy,const uint64_t InitialN);


void computeInner_stencil2dbk(double *Initial,double *share,uint64_t nRowsChunk,uint64_t nColsChunk,const uint64_t nRows, const uint64_t nCols);


void computeInner_stencil2dv2(bool IsLeftCBk,bool IsLastRBk,double *Initial,double *share,uint64_t nRowsChunk,uint64_t nColsChunk,const uint64_t nRows, const uint64_t nCols,double *upper, double *lower, double *current,double *left,double *right);

void StencilSeq ( double* __restrict__ dst,    double* __restrict__ src,const size_t     n_rows, const size_t     n_cols,const size_t     n_tsteps );

void print_results(const double *results, const size_t  n_rows_st,const size_t n_rows_ed, const size_t  n_cols_st, const size_t n_cols_ed,const size_t n_rows, const size_t n_cols);

#endif
