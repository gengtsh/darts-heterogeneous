extern "C" {
#include <cuda.h>
#include "conf.h"
#include "stencil.h"
}
#include <stdio.h>
#define ROTATE_DOWN(val,MAX) ((val-1==-1)?MAX-1:val-1)
#define ROTATE_UP(val,MAX) ((val+1)%MAX)
/**
  * GPU Device kernel for the for 2D stencil
  * First attempt during hackaton
  * M = Rows, N = Cols INCLUDING HALOS
  * In this version now we replace the size of the shared memory to be just 3 rows (actually 1+HALO*2) rows 
  */

//__global__ void gpu_stencil2D_4pt_hack4(double * dst, double * src, int M, int N)
//{
////	printf("kernel begin!\n");
//	//Declaring the shared memory array for source
//	extern	__shared__ double shared_mem[] ;
//
//	//indexes
//	int i, j, k,curRow;
//                           //Cols   *  numRows/Tile * tileIndex  
//	int base_global_row = ( N ) * ( GRID_TILE_Y * blockIdx.y ); 
//	int base_global_col = blockDim.x*blockIdx.x;
//	int base_global_idx = base_global_row + base_global_col ;
//	int center = 1,north = 0,south = 2; //indexes for the current location in the shared memory
//	int t = threadIdx.x;
//	
//	//copy the shared memory to fill the pipeline
//	bool rowLeft = (blockIdx.y==(gridDim.y-1))&&(M%GRID_TILE_Y<3)&&(M%GRID_TILE_Y>0);
//	int  numRowLeft =(rowLeft)?(3-M%GRID_TILE_Y):0;
//	bool noColsLeft = (base_global_col +t )<N;
//	bool noColsLeft2= (base_global_col+t+2)<N;
//	for (i = 0 ; i < 1+HALO*2-numRowLeft ; i ++ ){
//		k = base_global_idx+i*N+t;
//		j = i*(blockDim.x+2) + t;
//		shared_mem [j] = (noColsLeft)?src[k]:0.0;
//		if((t<2) &&(noColsLeft)){
//			shared_mem[j+blockDim.x]=src[k+blockDim.x];
//		}
//	}
//		
//	__syncthreads();
//
//	int tt = (((blockIdx.y+1)*GRID_TILE_Y)>M)?(M%GRID_TILE_Y): GRID_TILE_Y;
//	int ss = (((M%GRID_TILE_Y)==1)&&(blockIdx.y ==(gridDim.y-2)))?-1:0;
//
//	int lastRow = ((blockIdx.y == (gridDim.y-1))?-1:1)+tt +ss ;
////	printf("lastRow:%d \n",lastRow );
//	//Pipelined copy one row and process it
//	for ( curRow = HALO; curRow < lastRow; curRow+=1 )
//	{
//		//Stencil computation
//							//top             + bottom              + left                + right
//		j = threadIdx.x+HALO;
//		if(noColsLeft2){
//			dst[base_global_idx + curRow*N + j] =(shared_mem[north*(blockDim.x+2)+j] + shared_mem[south*(blockDim.x+2)+j] + shared_mem[center*(blockDim.x+2)+j-1] + shared_mem[center*(blockDim.x+2)+j+1] )/5.5;
//		}
//		__syncthreads();
//		//We are copying from src to shared memory.
//		k=base_global_idx+(curRow+2)*N+threadIdx.x;
//		if(k<M*N){
//			shared_mem [north*(blockDim.x+2)+threadIdx.x] =(noColsLeft)? src[k]:0.0;
//		}
//		if((t<2)&&(noColsLeft)&&(k<M*N)){
//			shared_mem[north*(blockDim.x+2)+threadIdx.x+blockDim.x]=src[k+blockDim.x];
//		}	
//		center = ROTATE_UP(center,3);
//		south  = ROTATE_UP(south,3);
//		north  = ROTATE_UP(north,3);
//		__syncthreads();
//	}
//
////	printf("kernel finish!\n");
//}



__global__ void gpu_stencil2D_4pt_hack4(double * dst, double * src, int M, int N)
{
//	printf("kernel begin!\n");
	//Declaring the shared memory array for source
	extern	__shared__ double shared_mem[] ;

	//indexes
	int i,j, k,curRow;
                           //Cols   *  numRows/Tile * tileIndex  
	int base_global_row = ( GRID_TILE_Y * blockIdx.y ); 
	int base_global_col = blockDim.x*blockIdx.x;
	int base_global_idx = N*base_global_row + base_global_col ;
	int center = 1,north = 0,south = 2; //indexes for the current location in the shared memory
	int t = threadIdx.x;
	
	//copy the shared memory to fill the pipeline
	bool legalCol = (base_global_col +t )<N;
	bool legalCol2= (base_global_col+t+2)<N;
	bool legalColn= (base_global_col+t+blockDim.x)<N;
	for (i = 0 ; i < 1+HALO*2 ; i ++ ){
		k = base_global_idx+i*N+t;
		j = i*(blockDim.x+2) + t;
		bool legalRow = (base_global_row+i)<M;
		shared_mem [j] =legalRow?( legalCol?src[k]:0.0):0.0;
		if((t<2)&&legalColn&&legalRow){
			shared_mem[j+blockDim.x]=src[k+blockDim.x];
		}
	}
		
	__syncthreads();

	//Pipelined copy one row and process it
	for ( curRow = HALO; curRow < GRID_TILE_Y+1; curRow+=1 )
	{
		//Stencil computation
		//top + bottom + left + right

		j = threadIdx.x+HALO;
		bool legalRow1 =( base_global_row+curRow+1)<M;
		if((legalCol2)&&(legalRow1)){
			dst[base_global_idx + curRow*N + j] =(shared_mem[north*(blockDim.x+2)+j] + shared_mem[south*(blockDim.x+2)+j] + shared_mem[center*(blockDim.x+2)+j-1] + shared_mem[center*(blockDim.x+2)+j+1] )/5.5;
		}
		__syncthreads();
		//We are copying from src to shared memory.
		int nextRow2 = base_global_row+curRow+2;
		bool legalRow2 = nextRow2<M;
		k = base_global_col+nextRow2*N+t;

		shared_mem [north*(blockDim.x+2)+t] =(legalRow2&&legalCol)?src[k]:0.0;

		if((t<2)&&legalColn&&legalRow2){
			shared_mem[north*(blockDim.x+2)+t+blockDim.x]=src[k+blockDim.x];
		}	
		
		center = ROTATE_UP(center,3);
		south  = ROTATE_UP(south,3);
		north  = ROTATE_UP(north,3);
		__syncthreads();
	}

//	printf("kernel finish!\n");
}


__global__ void gpu_stencil2D_4pt_hack2(double * dst, double * src, int M, int N)
{
//	printf("kernel begin!\n");
	//Declaring the shared memory array for source
	__shared__ double shared_mem[ 1 + HALO*2 ] [ GRID_TILE_X + HALO*2]; //1 is the row I am modifying
	//double * shSrc = shared_mem;

	//indexes
	int i, j, curRow;
                           //Cols   *  numRows/Tile * tileIndex  
	int base_global_idx = ( N ) * ( GRID_TILE_Y * blockIdx.y ) + GRID_TILE_X*blockIdx.x;
	
	int center = 1,north = 0,south = 2; //indexes for the current location in the shared memory

	//copy the shared memory to fill the pipeline
	for (i = 0 ; i < 1+HALO*2 ; i ++ )
		for (j = threadIdx.x ; j < GRID_TILE_X+2*HALO ; j+=blockDim.x)
		{
			shared_mem [i][j] = src[base_global_idx + i*N + j];
		}
	__syncthreads();
	//Pipelined copy one row and process it
	for ( curRow = HALO; curRow < GRID_TILE_Y; curRow+=1 )
	{
		//Stencil computation
		for (j = threadIdx.x + HALO ; j < GRID_TILE_X+HALO ; j+=blockDim.x)
		{
							//top             + bottom              + left                + right
			dst[base_global_idx + curRow*N + j] = (shared_mem[north][j] + shared_mem[south][j] + shared_mem[center][j-1] + shared_mem[center][j+1] )/5.5;
		}
		
		__syncthreads();
		//We are copying from dst to shared memory.
		for (j = threadIdx.x ; j < GRID_TILE_X+2*HALO ; j+=blockDim.x)
		{
			shared_mem [north][j] = src[base_global_idx + (curRow+2)*N + j];
		}
	
		center = ROTATE_UP(center,3);
		south  = ROTATE_UP(south,3);
		north = ROTATE_UP(north,3);
		__syncthreads();
	}

	//Dranning the pipeline
	for (j = threadIdx.x + HALO ; j < GRID_TILE_X+HALO ; j+=blockDim.x)
	{
							//top             + bottom              + left                + right
		dst[base_global_idx + curRow*N + j] = (shared_mem[north][j] + shared_mem[south][j] + shared_mem[center][j-1] + shared_mem[center][j+1] )/5.5;
	}
	__syncthreads();

//	printf("kernel finish!\n");
}



///**
//  * GPU Device kernel for the for 2D stencil
//  * First attempt during hackaton
//  * M = Rows, N = Cols INCLUDING HALOS
//  */
//__global__ void gpu_stencil2D_4pt_hack1(double * dst, double * src, int M, int N)
//{
//
//	//Declaring the shared memory array for source
//	__shared__ double shared_mem[GRID_TILE_Y + HALO*2 ] [ GRID_TILE_X + HALO*2];
//	//double * shSrc = shared_mem;
//
//	//indexes
//	int i, j;
//
//                           //Cols   *  numRows/Tile * tileIndex  
//	int base_global_idx = ( N ) * ( GRID_TILE_Y * blockIdx.y ) + GRID_TILE_X*blockIdx.x;
//
//	//We are copying from dst to shared memory.
//	for (i = 0 ; i < GRID_TILE_Y+2*HALO ; i ++ )
//		for (j = threadIdx.x ; j < GRID_TILE_X+2*HALO ; j+=blockDim.x)
//		{
//			shared_mem [i][j] = src[base_global_idx + i*N + j];
//		}
//
//	__syncthreads();
//
//	//Stencil computation
//	for (i = HALO ; i < GRID_TILE_Y+HALO ; i ++ )
//		for (j = threadIdx.x + HALO ; j < GRID_TILE_X+HALO ; j+=blockDim.x)
//		{
//			                                //top             + bottom              + left                + right
//			dst[base_global_idx + i*N + j] = (shared_mem[i-1][j] + shared_mem[i+1][j] + shared_mem[i][j-1] + shared_mem[i][j+1] )/5.5;
//		}
//
//	__syncthreads();
//}




/**
  * GPU Device kernel for the for 2D stencil
  * M = Rows, N = Cols
  */
__global__ void gpu_stencil2D_4pt(double * dst, double * src, int M, int N)
{
	//Declaring the shared memory array for source
	extern __shared__ double shared_mem[];
	double * shSrc = shared_mem;

	//indexes
	int i, j;

	//neighbor's values 
	double north, south, east, west;



	//SharedMem Collumns Dimension
	int smColDim = HALO*2+blockDim.y*TILE_SIZE;
	int smRowDim = HALO*2+blockDim.x*TILE_SIZE;

	//Copying to shared memory

	//Inner part
	for ( i = 0 ; i < TILE_SIZE ; i++ )
	{
		for ( j = 0 ; j < TILE_SIZE ; j++ )
		{
			int globalIndex=HALO*N+blockIdx.x*blockDim.x*TILE_SIZE*N+threadIdx.x*TILE_SIZE*N+i*N+blockIdx.y*blockDim.y*TILE_SIZE+threadIdx.y*TILE_SIZE+j+HALO;
			int shMemIndex=HALO*smColDim+threadIdx.x*smColDim*TILE_SIZE+i*smColDim+HALO+threadIdx.y*TILE_SIZE+j;
			shSrc[shMemIndex]=src[globalIndex];
		}
	}

	//Halos

	if (threadIdx.x == 0 && threadIdx.y == 0 ) 
	{

		int indexTopHalo, indexBottomHalo, indexLeftHalo, indexRightHalo;
		//For Bottom and top row
		for ( i = 0 ; i < HALO ; i++ )
		{
			for ( j = 0 ; j < smColDim ; j++ )
			{
				indexTopHalo = (blockIdx.x*blockDim.x*TILE_SIZE+i)*N + (blockIdx.y*blockDim.y*TILE_SIZE) + j;
				indexBottomHalo = (HALO + (blockIdx.x+1)*blockDim.x*TILE_SIZE)*N + (blockIdx.y*blockDim.y*TILE_SIZE)+j;
				shSrc[i*smColDim+j] = src[indexTopHalo];
				shSrc[(HALO+blockDim.x*TILE_SIZE+i)*smColDim + j] = src[indexBottomHalo];
			}
		}
		
		//For right and left Columns
		for ( i = 0 ; i < HALO ; i++ )
		{
			for ( j = 0 ; j < smRowDim-HALO*2; j ++ )
			{
				indexLeftHalo = (HALO+blockIdx.x*blockDim.x*TILE_SIZE+j)*N + (blockIdx.y*blockDim.y*TILE_SIZE)+i;
				indexRightHalo = (HALO+blockIdx.x*blockDim.x*TILE_SIZE+j)*N + ((blockIdx.y+1)*blockDim.y*TILE_SIZE)+HALO+i;
				shSrc[(HALO+j)*smColDim+i] = src[indexLeftHalo];
				shSrc[(HALO+j+1)*smColDim-HALO+i] = src[indexRightHalo];
			}
		}
	}

	__syncthreads();



	for ( i = 0 ; i < TILE_SIZE ; i++ )
	{
		for ( j = 0 ; j < TILE_SIZE ; j++ )
		{
			int globalIndex=HALO*N+blockIdx.x*blockDim.x*TILE_SIZE*N+threadIdx.x*TILE_SIZE*N+i*N+blockIdx.y*blockDim.y*TILE_SIZE+threadIdx.y*TILE_SIZE+j+HALO;
			int shMemIndex=HALO*smColDim+threadIdx.x*smColDim*TILE_SIZE+i*smColDim+HALO+threadIdx.y*TILE_SIZE+j;


			//Getting the neighbohrs
			north = shSrc[shMemIndex-smColDim];
			south = shSrc[shMemIndex+smColDim];
			east  = shSrc[shMemIndex+1];
			west  = shSrc[shMemIndex-1];
			//Real Stencil operation
			dst[globalIndex] = ( north + south + east + west )/5.5;
//			dst[globalIndex] = ( north + south + east + west )/4;
		}
	}

	__syncthreads();
}




/**
 *  Naïve 4pt stencil code for 2D arrays. 
 */
void
stencil2D4pt ( double* __restrict__ dst,    double* __restrict__ src, 
               const size_t     n_rows, const size_t     n_cols,
               const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];
    volatile Array2D DST = (Array2D) dst,
            SRC = (Array2D) src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
        for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1])/5.5;
            }
        }
        SWAP_PTR(&DST,&SRC);
    }
}

extern "C"
void
stencil2D4pt_gpu( double * __restrict__ dst, double* __restrict__ src,
		  const size_t M, const size_t N, 
		  const size_t NUM_ITERATIONS)//M Rows by N Columns
{
		
	double size = sizeof(double) * M * N;

	//device memory allocation
	double * d_dst, * d_src;
	cudaMalloc( (void **) &d_dst, size);
	cudaMalloc( (void **) &d_src, size);
	
	//dimmensions for indexes
	// TODO the -2 is to remove the borders
	dim3 dimBlock(MAX_BLOCK_DIM,MAX_BLOCK_DIM);
	int gridx = (N-2)/(MAX_BLOCK_DIM*TILE_SIZE) + (((N-2)%(MAX_BLOCK_DIM*TILE_SIZE) == 0)? 0:1 ) ;
	int gridy = (M-2)/(MAX_BLOCK_DIM*TILE_SIZE) + (((M-2)%(MAX_BLOCK_DIM*TILE_SIZE) == 0)? 0:1 ) ;
	dim3 dimGrid(gridx,gridy);

	//Shared memory size = inside + halo
	int shMemSize=MAX_BLOCK_DIM*TILE_SIZE*MAX_BLOCK_DIM*TILE_SIZE*sizeof(double)+(HALO*MAX_BLOCK_DIM*TILE_SIZE+HALO*HALO)*4*sizeof(double);
	
	//Hackaton dimensions
	dim3 dimGrid_hack1((N-HALO*2)/GRID_TILE_X,(M-HALO*2)/GRID_TILE_Y);

	//Copying the device memory
	cudaMemcpy(d_src, src, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_dst, dst, size, cudaMemcpyHostToDevice);

	//printf("CUDA Stencil Code running... cycles = %d. dim = %d by %d \n",NUM_ITERATIONS,M,N);
	
    int time_step = NUM_ITERATIONS;

    while (time_step-- > 0) 
    {
    	//gpu_stencil2D_4pt<<<dimGrid,dimBlock,shMemSize>>>(d_dst,d_src,M,N);
		//gpu_stencil2D_4pt_hack1<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N); //JOSE Hackathon!
		//printf("before: d_src[10] = %ld",d_src[10]);

		gpu_stencil2D_4pt_hack2<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N);
		//Inline swapping.
		
		//printf("after: d_src[10] = %ld",d_src[10]);
		double * temp;
		if ( NUM_ITERATIONS%2 ==0 || time_step !=0)
		{
			temp=d_src;
			d_src=d_dst;
			d_dst=temp;
		}
	}
	
	
	//Copying memory back from device to DRAM
	//cudaMemcpy(src, d_src, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(dst, d_dst, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(src, d_src, size, cudaMemcpyDeviceToHost);
	
	//Free device memory
	cudaFree(d_src); cudaFree(d_dst);
}

//void*
//stencil_run(void* arg)
//{
//    stencil_t* stencil = (stencil_t*)arg;
//    STENCIL_COMPUTE(stencil->stencil,stencil->arg);
//    return NULL;
//}


void gpu_kernel4(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * d_src, int M, int N){
		int sharedMemSize = sizeof(double)*(1+HALO*2)*(dimBlock.x+2);
#ifdef CUDA_DARTS_DEBUG
		printf("sharedMemSize: %dKB, total sharedMemSize: %dKB\n",sharedMemSize/1024, sharedMemSize*dimGrid.x*dimGrid.y/1024);
#endif
		gpu_stencil2D_4pt_hack4<<<dimGrid,dimBlock,sharedMemSize>>>(d_dst,d_src,M,N);
}

void gpu_kernel1(dim3 dimGrid_hack1,double * d_dst, double * d_src, int M, int N){
		gpu_stencil2D_4pt_hack2<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N);
}

void gpu_kernel3(cudaStream_t &stream,dim3 dimGrid_hack1,double * d_dst, double * d_src, int M, int N){
		gpu_stencil2D_4pt_hack2<<<dimGrid_hack1,NUM_THREADS,0,stream>>>(d_dst,d_src,M,N);

}

void gpu_kernel2(dim3 dimGrid_hack1,double *dst, double *src, double size, size_t ts, double * d_dst, double * d_src, int M, int N){
	double * tmp;
	while (--ts!=0){
		printf("ts:%ld \n", ts);
		gpu_stencil2D_4pt_hack2<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N);
		tmp = d_src;
		d_src = d_dst;
		d_dst=tmp;
	}
}
