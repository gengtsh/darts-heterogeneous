//SAMPLE OF A 2D STENCIL CODE
//Made by: Jose Monsalve

#include <stdio.h>
#include <stdlib.h>
#include "conf.h"
#define NUM_ITERATIONS 10

void SWAP(void ** dst, void ** src);
void cpu_stencil2D_4pt(double * o_dst, double * o_src, const int N);

int main()
{
        double * dst = malloc( sizeof(double) * SIZE_MATRIX * SIZE_MATRIX) ;//[SIZE_MATRIX][SIZE_MATRIX];
        double * src = malloc( sizeof(double) * SIZE_MATRIX * SIZE_MATRIX) ;//[SIZE_MATRIX][SIZE_MATRIX];
        for (int i = 0 ; i < SIZE_MATRIX; i++)
        {
                for (int j = 0 ; j < SIZE_MATRIX ; j++)
		{
                        src[i*SIZE_MATRIX+j]=1;
			dst[i*SIZE_MATRIX+j]=0;
		}
        }
        cpu_stencil2D_4pt(dst,src,SIZE_MATRIX);
#if PRINT_RESULTS == 1 
        for (int i = 0 ; i < SIZE_MATRIX; i++)
        {
                for (int j = 0 ; j < SIZE_MATRIX ; j++)
                        printf(" %f ",src[i*SIZE_MATRIX + j]);
                printf("\n");
        }
        return 0;
#endif
}

void cpu_stencil2D_4pt(double * o_dst, double * o_src, const int N)
{
        int time_step = NUM_ITERATIONS;

        double (*dst)[N] = (double (*)[N]) o_dst; //destination matrix, where the new computes are going to be stores
        double (*src)[N] = (double (*)[N]) o_src; //source matrix, where the samples are taken from the previous iteration

        while (time_step > 0) 
        {
                for (size_t i = 1; i < N-1; ++i) 
                {
                        for (size_t j = 1; j < N-1; ++j) 
                        {
                                dst[i][j] = ( src[i][j-1] + src[i][j+1] + src[i-1][j] + src[i+1][j] ) /4;
			} 
                }
                SWAP((void **) &dst,(void **) &src);
                time_step--;
        }
}

void SWAP(void ** dst, void ** src)
{
    void * temp = *dst;
    *dst = *src;
    *src = temp;
} 
