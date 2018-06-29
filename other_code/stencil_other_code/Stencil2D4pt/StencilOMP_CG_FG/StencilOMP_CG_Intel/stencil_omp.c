
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>

#include "stencil.h"
#include <pmmintrin.h>

void test_threads(){
#   pragma omp parallel
    {
        printf("test: threads # %d in %d threads\n", omp_get_thread_num(),omp_get_num_threads());
    }
}


/**
 *  Naïve 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt ( double* restrict dst,    double* restrict src, 
               const size_t     n_rows, const size_t     n_cols,
               const size_t     n_tsteps )
{
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

/**
 *  Naïve parallelization of a 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt_omp ( double* restrict dst,    double* restrict src, 
                   const size_t     n_rows, const size_t     n_cols,
                   const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];

    Array2D DST = (Array2D) dst,
            SRC = (Array2D) src;

    for (size_t ts = 0; ts < n_tsteps; ++ts) {
#       pragma omp parallel for default(none) shared(DST,SRC) firstprivate(n_rows,n_cols,n_tsteps) 
        for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
            }
        }
        SWAP_PTR(&DST,&SRC);
    }
}


void
init_loop ( double* dst, double *src, const size_t n_rows, const size_t n_cols )
{
 
    typedef double (*Array2D)[n_cols];

        Array2D DST = (Array2D) dst,
                SRC = (Array2D) src;

#pragma omp for nowait schedule(static) 
        for (size_t i = 0; i < n_rows; ++i) 
            for (size_t j = 0; j < n_cols; ++j){ 
                DST[i][j]=SRC[i][j];
            }
    
}

/**
 *  Less naïve parallelization of a 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt_omp_v2 ( double* restrict dst,    double* restrict src, 
                      const size_t     n_rows, const size_t     n_cols,
                      const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];

#       pragma omp parallel default(none) shared(src, dst) firstprivate(n_rows,n_cols,n_tsteps)
    {
        Array2D DST = (Array2D) dst,
                SRC = (Array2D) src;
        size_t n_ts = n_tsteps;

        while (n_ts-- > 0) {
//#           pragma omp for nowait
#           pragma omp for nowait schedule(static)
            for (size_t i = 1; i < n_rows-1; ++i) {
                for (size_t j = 1; j < n_cols-1; ++j) {
                    DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
                }
            }
            SWAP_PTR(&DST,&SRC);
//#           pragma omp single nowait
//            {
//              --n_ts;
//            } // barrier here => flush!
#           pragma omp barrier
//#           pragma omp flush
        }
    }
}

void
stencil2D4pt_omp_simd ( double* restrict dst,    double* restrict src, 
                      const size_t     n_rows, const size_t     n_cols,
                      const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];

#       pragma omp parallel default(none) shared(src, dst) firstprivate(n_rows,n_cols,n_tsteps)
    {
        Array2D DST = (Array2D) dst,
                SRC = (Array2D) src;
        size_t n_ts = n_tsteps;

        while (n_ts-- > 0) {
//#           pragma omp for nowait
#           pragma omp for nowait schedule(static)
            for (size_t i = 1; i < n_rows-1; ++i) {
#               pragma omp simd 
                for (size_t j = 1; j < n_cols-1; ++j) {
                    DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
                }
            }
            SWAP_PTR(&DST,&SRC);
//#           pragma omp single nowait
//            {
//              --n_ts;
//            } // barrier here => flush!
#           pragma omp barrier
//#           pragma omp flush
        }
    }
}


void
stencil2D4pt_simd2(double* restrict dst,    double* restrict src, 
                   const size_t     n_rows, const size_t     n_cols,
                   const size_t     n_tsteps)
{
    typedef double (*Array2D)[n_cols];

    Array2D DST = (Array2D) dst,
            SRC = (Array2D) src;

    for (size_t ts = 0; ts < n_tsteps; ++ts) {
        for (size_t i = 1; i < n_rows-1; ++i) {
            size_t j;
            for (j = 1; j < n_cols-1 - 1; j+=2) { // from 1 .. (bound - unroll factor)
// Iter: j+0    DST[ i ][j+0] = (SRC[i-1][j+0] + SRC[i+1][j+0] + SRC[i][j-1+0] + SRC[i][j+1+0]) / 4;
// Iter: j+1    DST[ i ][j+1] = (SRC[i-1][j+1] + SRC[i+1][j+1] + SRC[i][j-1+1] + SRC[i][j+1+1]) / 4;
// Example with i = 1, j = 1:
//              DST[1,1] = (SRC[0,1] + SRC[2,1] + SRC[1,0] + SRC[1,2]) / 4;
//              DST[1,2] = (SRC[0,2] + SRC[2,2] + SRC[1,1] + SRC[1,3]) / 4;

                /* 
                 * DATA LOADING 
                 * Note: vsh1 and vsh2 contain values required by either
                 * iteration (i,j) or iteration (i,j+1) 
                 */
                __m128d v00  = _mm_load_pd ( &SRC[i-1][j-1] ), // v00  = [ (0,0), (0,1) ]
                        v01  = _mm_load_pd ( &SRC[i+1][j-1] ), // v01  = [ (2,0), (2,1) ]
                        v10  = _mm_load_pd ( &SRC[i-1][j+1] ), // v10  = [ (0,2), (0,3) ]
                        v11  = _mm_load_pd ( &SRC[i+1][j+1] ), // v11  = [ (2,2), (2,3) ]
                        vsh1 = _mm_load_pd ( &SRC[ i ][j-1] ), // vsh1 = [ (1,0), (1,1) ]
                        vsh2 = _mm_load_pd ( &SRC[ i ][j+1] ); // vsh2 = [ (1,2), (1,3) ]

                /*
                 * DATA SHUFFLING
                 * Build v00, v01, v10, v11, so that the right values are in the
                 * right slot
                 */
                v00          = _mm_shuffle_pd ( v00,  v01,  _MM_SHUFFLE2(1,1) ); // v00  = [ (0,1), (2,1) ]
                v01          = _mm_shuffle_pd ( vsh1, vsh2, _MM_SHUFFLE2(0,0) ); // v01  = [ (1,0), (1,2) ]
                v10          = _mm_shuffle_pd ( v10,  v11,  _MM_SHUFFLE2(0,0) ); // v10  = [ (0,2), (2,2) ]
                v11          = _mm_shuffle_pd ( vsh1, vsh2, _MM_SHUFFLE2(1,1) ); // v11  = [ (1,1), (1,3) ]

                /*
                 * DATA PROCESSING 
                 * We perform a partial sum, THEN a horizontal sum directly in
                 * the result vector <res>. We can then divide each member of
                 * <res> by 4.
                 */

                __m128d res, div_by;
                div_by = _mm_set_pd  ( 4.0,    4.0 );
                v00    = _mm_add_pd  ( v00,    v01 ); // [ (v00.left + v01.left), (v00.right + v01.right) ]
                v10    = _mm_add_pd  ( v10,    v11 ); // [ (v10.left + v11.left), (v10.right + v11.right) ]
                res    = _mm_hadd_pd ( v00,    v10 ); // [ (v00.left + v00.right), (v10.left + v10.right) ]
                res    = _mm_div_pd  ( res, div_by );

                /* 
                 * DATA WRITE BACK
                 * We are writing to "odd" cells in the array (read: unaligned).
                 * So we need to use the 'u' version of the SSE intrinsics,
                 * otherwise we will crash the program
                 * Possible optimizations include unrolling by a bigger factor
                 * (e.g., 4), and thus possibly manage to use more aligned
                 * stores. Not likely to happen:
                 * 1 iteration = 3 "private" vector regs + 2 shared "vector" regs
                 */
                _mm_storeu_pd ( &DST[i][j], res );

            }

            for (; j < n_cols-1; ++j) { // epilog
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
            }
        }

        SWAP_PTR(&DST,&SRC);
    }
}

void*
stencil_run(void* arg)
{
    stencil_t* stencil = (stencil_t*)arg;
    STENCIL_COMPUTE(stencil->stencil,stencil->arg);
    return NULL;
}

