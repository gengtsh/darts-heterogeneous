//#ifdef __cplusplus
//extern "C" {
//#endif

#include <cuda.h>
#include <stdio.h>
extern "C"{
#include "Hello_cuda.h"
}

const int N=32;
const int M=1000;
const int blocksize =32;

__global__
void hello(int *a){
	int tx=threadIdx.x;
	printf ("tx: %d \n",tx);
	a[tx]=tx*tx;	
}

extern "C"
void Hello_cuda()
{
	size_t sz =N*sizeof(int);
	int *a= (int *)malloc(sz);
	for (int i=0;i<N;++i){
			a[i]=2;
	}
	
	int *d_addr;
	cudaMalloc( (void **)&d_addr,sz );
	cudaMemcpy( d_addr, a, sz, cudaMemcpyHostToDevice );

	dim3 dimBlock( blocksize, 1 );
	dim3 dimGrid( 1, 1 );
	hello<<<dimGrid, dimBlock>>>(d_addr);
	cudaMemcpy( a, d_addr, sz, cudaMemcpyDeviceToHost );
	cudaFree( d_addr );

	for (int i=0;i<N;++i){
			if(a[i]!= i*i){
				fprintf(stderr,"a[%d] = %d != %d !\n",i,a[i],i*i)	;
			}
	}
	free(a);

}



//#ifdef __cplusplus
//}
//#endif
