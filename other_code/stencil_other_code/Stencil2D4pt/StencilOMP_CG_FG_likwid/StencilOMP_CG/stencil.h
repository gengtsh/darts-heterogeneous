#ifndef STENCIL_H_GUARD
#define STENCIL_H_GUARD
#include "util.h"
#include <omp.h>

#define SWAP(type,left,right) do { \
    type tmp = left;               \
    left     = right;              \
    right    = tmp;                \
} while(0)

static inline void swap_ptr(void** left, void** right) {
    void* tmp = *left;
    *left     = *right;
    *right    = tmp;
}

#define SWAP_PTR(left,right) swap_ptr((void**)left,(void**)right)


//#define THREAD_TEST 1
#define GRID_TILECPU_X 48 
#define GRID_TILECPU_Y 48

typedef struct stencil_arg_s {
    double *dst, 
		   *src;
    size_t  n_rows, 
			n_cols,
			n_tsteps;
} stencil_arg_t;

typedef void (*stencil_funptr_t)(double*, double*, const size_t, const size_t, const size_t);
typedef struct stencil_s {
    stencil_arg_t *arg;
    stencil_funptr_t stencil;
} stencil_t;

#define STENCIL_COMPUTE(c,d) (c)((d)->dst,(d)->src,(d)->n_rows,(d)->n_cols,(d)->n_tsteps)
void* stencil_run(void* arg);

void stencil_init        ( stencil_t    *stencil, stencil_funptr_t code, 
                           const size_t  n_rows,  const size_t     n_cols, 
                           size_t        n_tsteps );
void stencil_arg_copy    ( stencil_arg_t* dest, stencil_arg_t* sourc );

void stencil2D4pt        ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );
void stencil2D4pt_omp    ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );


void init_loop ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols );

void stencil2D4pt_omp_v2 ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );


void stencil2D4pt_omp_simd ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );

void stencil2D4pt_simd2  ( double* restrict dst,    double* restrict src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );

void test_threads();


void stencil2D4ptSeq_wb ( double*__restrict__ dst,    double* __restrict__ src, const size_t     n_rows, const size_t     n_cols,const size_t     n_tsteps );

void computeBlock_stencil25(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_rows_ck,size_t n_cols_ck);

#endif // STENCIL_H_GUARD
