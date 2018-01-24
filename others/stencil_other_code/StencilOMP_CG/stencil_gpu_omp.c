#include "stencil.h"

#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_thread() 0
#endif


/**
 *  Naïve 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt ( double* restrict dst,    double* restrict src, 
               const size_t     n_rows, const size_t     n_cols,
               const size_t     n_tsteps )
{
    /* 
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
    */
    printf("Computing sequential stencil\n");
}

/**
 *  Naïve 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt_gpu_kernel ( double* restrict dst,    double* restrict src, 
                          const size_t     n_rows, const size_t     n_cols,
                          const size_t     n_tsteps )
{
    /*
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
    */
    printf("Sending data to board: Processing %lu rows\n", n_rows);
    printf("Launching GPU kernel\n");
}

void
stop_gpu(void)
{
    printf("Sync-ing all GPU threads\n");
    printf("Copying data back from the board\n");
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
    size_t n_ts = n_tsteps;

    printf("Checking whether we need the GPU...\n");

#   pragma omp parallel default(none) shared(src, dst, n_ts) firstprivate(n_rows,n_cols,n_tsteps)
    {
        Array2D DST = (Array2D) dst,
                SRC = (Array2D) src;
        char message[101];

        while (n_ts > 0) {
#           pragma omp master
            {
                size_t gpu_n_rows = (2 * n_rows) / 3;
                stencil2D4pt_gpu_kernel((double*)DST, (double*)SRC, gpu_n_rows, n_cols, 1UL); // n_tsteps is ignored
            }
#           pragma omp for nowait
            for (size_t i = (2 * n_rows) / 3; i < n_rows-1; ++i) {
                snprintf(message, sizeof(message), "[%d]\tcomputing a stencil row %lu", omp_get_thread_num(), i);
                puts(message);
//              for (size_t j = 1; j < n_cols-1; ++j) {
//                  DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 4;
//              }
            }
            SWAP_PTR(&DST,&SRC);
            // Keep thread ID #0 as the "GPU controller"
#           pragma omp master
            {
                stop_gpu();
            } 
            // Anyone can decrement the time step counter
#           pragma omp single
            {
                --n_ts;
            } // barrier -> flush(n_ts) !
        }
    }
}
