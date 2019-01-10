
#ifndef SPENCILCPUKERNEL_H
#define SPENCILCPUKERNEL_H

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstdio>

void computeInner_stencil2d_v2(double *dst,double *src,size_t n_rows,size_t n_cols);

void copyColumns_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride );

void copyRows_stencil2d(double *dst,double *src,size_t n_rows,size_t n_cols,size_t stride);


void computeInner_stencil37(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_slices);


void computeBlock_stencil37(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_slices,size_t n_rows_ck,size_t n_cols_ck,size_t n_slices_ck);

void calcNextEP(int *arrncEdge, int *arrngEdge,int *arrncEdgeLeft, int *arrngEdgeLeft,int *cpos, int *gpos,int *arrnEdge, int *cvarEdge, int *gvarEdge,int *arrncEdgeMin,int *arrngEdgeMin,int cCnt,int gCnt, int nId,char *home);

void calcNextEEL(int * arrncgEdge, int *arrncgEdgeLeft,int *arrnEdge, int *varEdge,int *arrncgEdgeMin,int id, int homeCnt, int visitCnt);


#endif
