#include "stencil.h"

#include <stdio.h>

#include <pmmintrin.h>

#include <omp.h>

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
#           pragma omp for nowait
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


//stencil fine grain
void stencil2D4pt_omp_v3 ( double* restrict dst,    double* restrict src, 
                      const size_t     n_rows, const size_t     n_cols,
                      const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];
	size_t nThrds = omp_get_num_threads();
	size_t nThrdsMinus1=nThrds-1;
	size_t chunk= (n_rows-2)/nThrds;
	

	int **dCalc;
	int **dSwap;

	dCalc = smalloc(sizeof(int*)*nThrds);
	dSwap = smalloc(sizeof(int*)*nThrds);
	for (size_t i=0;i<nThrds;++i){
		dCalc[i]=smalloc(sizeof(int)*(n_tsteps+1));
		dSwap[i]=smalloc(sizeof(int)*(n_tsteps+1));
	}

	Array2D *dstPtr  = smalloc(sizeof(Array2D)*nThrds);
	Array2D *srcPtr  = smalloc(sizeof(Array2D)*nThrds);

	#pragma omp parallel shared(src,dst,dCalc,dSwap,dstPtr,srcPtr) firstprivate(n_rows,n_cols,n_tsteps,chunk,nThrds,nThrdsMinus1)
	{
		#pragma omp single
		{

			#pragma omp taskgroup
			{
				#pragma omp task
				{
					
					for(size_t i=0;i<nThrds;++i){
						size_t pos1=1+chunk*i*n_cols+n_cols;
						
						dstPtr[i] = (Array2D) (dst+pos1);
						srcPtr[i] = (Array2D) (src+pos1);

					}
				}
				#pragma omp taskwait
				
				
				for(int ts=n_tsteps;ts>0;--ts){
					//calc 0 chunk
					#pragma omp task depend(out:dCalc[0][ts]) depend(in:dSwap[0][ts+1]) depend(in:dSwap[1][ts+1]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
					{
						
						for (size_t i = 0; i < chunk; ++i) {
            			    for (size_t j = 0; j < n_cols-2; ++j) {
            			        dstPtr[0][i][j] = (srcPtr[0][i-1][j] + srcPtr[0][i+1][j] + srcPtr[0][i][j-1] + srcPtr[0][i][j+1]) / 4;
            			    }
            			}

				
					//	printf("compute,ts=%d \n",ts);

					//	printf("src0:\n");
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][-1][-1],srcPtr[0][-1][0],srcPtr[0][-1][1], srcPtr[0][-1][2]);
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][0][-1],srcPtr[0][0][0],srcPtr[0][0][1], srcPtr[0][0][2]);
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][1][-1],srcPtr[0][1][0],srcPtr[0][1][1], srcPtr[0][1][2]);

					//	printf("dst0:\n");
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][-1][-1],dstPtr[0][-1][0],dstPtr[0][-1][1], dstPtr[0][-1][2]);
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][0][-1],dstPtr[0][0][0],dstPtr[0][0][1], dstPtr[0][0][2]);
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][1][-1],dstPtr[0][1][0],dstPtr[0][1][1], dstPtr[0][1][2]);
					
					//	printf("src0:\n");
					//	printf("%f,%f,%f,%f, \n",src[0],src[1],src[2], src[3]);
					//	printf("%f,%f,%f,%f, \n",src[0+n_cols],src[1+n_cols],src[2+n_cols], src[3+n_cols]);
					//	printf("%f,%f,%f,%f, \n",src[0+n_cols*2],src[1+n_cols*2],src[2+n_cols*2], src[3+n_cols*2]);
					//	
					//	printf("dst0:\n");
					//	printf("%f,%f,%f,%f, \n",dst[0],dst[1],dst[2], dst[3]);
					//	printf("%f,%f,%f,%f, \n",dst[0+n_cols],dst[1+n_cols],dst[2+n_cols], dst[3+n_cols]);
					//	printf("%f,%f,%f,%f, \n",dst[0+n_cols*2],dst[1+n_cols*2],dst[2+n_cols*2], dst[3+n_cols*2]);
					
					}

					//calc last chunk
					#pragma omp task depend(out:dCalc[nThrdsMinus1][ts]) depend(in:dSwap[nThrdsMinus1-1][ts+1]) depend(in:dSwap[nThrdsMinus1][ts+1]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
					{
						size_t pos2=n_rows-2-chunk*nThrdsMinus1;
						
						for (size_t i = 0; i < pos2; ++i) {
            			    for (size_t j = 0; j < n_cols-2; ++j) {
            			        dstPtr[nThrdsMinus1][i][j] = (srcPtr[nThrdsMinus1][i-1][j] + srcPtr[nThrdsMinus1][i+1][j] + srcPtr[nThrdsMinus1][i][j-1] + srcPtr[nThrdsMinus1][i][j+1]) / 4;
            			    }
            			}
					}

					//calc others chunks
					for(size_t k=1;k<nThrdsMinus1;++k){
						#pragma omp task depend(out:dCalc[k][ts]) depend(in:dSwap[k-1][ts+1]) depend(in:dSwap[k][ts+1]) depend(in:dSwap[k+1][ts+1]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
						{
							for (size_t i = 0; i < chunk; ++i) {
            				    for (size_t j = 0; j < n_cols-2; ++j) {
									dstPtr[k][i][j] = (srcPtr[k][i-1][j] + srcPtr[k][i+1][j] + srcPtr[k][i][j-1] + srcPtr[k][i][j+1]) / 4;
								}
            				}
						}
					}

					//swap 0 chunk
					#pragma omp task depend(out:dSwap[0][ts]) depend(in:dCalc[0][ts]) depend(in:dCalc[1][ts]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
					{

						
						SWAP_PTR(&(dstPtr[0]),&(srcPtr[0]));

					//	printf("SWAP,ts=%d\n",ts);
				
					//	printf("src0:\n");
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][-1][-1],srcPtr[0][-1][0],srcPtr[0][-1][1], srcPtr[0][-1][2]);
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][0][-1],srcPtr[0][0][0],srcPtr[0][0][1], srcPtr[0][0][2]);
					//	printf("%f,%f,%f,%f, \n",srcPtr[0][1][-1],srcPtr[0][1][0],srcPtr[0][1][1], srcPtr[0][1][2]);
					//	
					//	printf("dst0:\n");
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][-1][-1],dstPtr[0][-1][0],dstPtr[0][-1][1], dstPtr[0][-1][2]);
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][0][-1],dstPtr[0][0][0],dstPtr[0][0][1], dstPtr[0][0][2]);
					//	printf("%f,%f,%f,%f, \n",dstPtr[0][1][-1],dstPtr[0][1][0],dstPtr[0][1][1], dstPtr[0][1][2]);
					}


					//swap last chunk
					#pragma omp task depend(out:dSwap[nThrdsMinus1][ts]) depend(in:dCalc[nThrdsMinus1-1][ts]) depend(in:dCalc[nThrdsMinus1][ts]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
					{
						SWAP_PTR(&(dstPtr[nThrdsMinus1]),&(dstPtr[nThrdsMinus1]));

					}
					//swap other chunks
					for(size_t k=1;k<nThrdsMinus1;++k)
					{
						#pragma omp task depend(out:dSwap[k][ts]) depend(in:dCalc[k-1][ts]) depend(in:dCalc[k][ts]) depend(in:dCalc[k+1][ts]) shared(dst,src,dstPtr,srcPtr) firstprivate(n_cols)
						{

							SWAP_PTR(&(dstPtr[k]),&srcPtr[k]);
						}
					}
					

				}
			}
		}
	}
	
	for (size_t i=0;i<nThrds;++i){
		sfree(dCalc[i]);
		sfree(dSwap[i]);
	}
	
	sfree(dCalc);
	sfree(dSwap);
	sfree(dstPtr);
	sfree(srcPtr);


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

