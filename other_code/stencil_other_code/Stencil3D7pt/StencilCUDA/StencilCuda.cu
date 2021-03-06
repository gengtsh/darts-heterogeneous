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
#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("kernel begin!\n");
	}
#endif
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
#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("kernel finish!\n");
	}
#endif
}


__global__ void gpu_stencil2D_4pt_hack5_cp_rows(double * dst, double * shared_cols, double *shared_rows,int tile_y,int M, int N){


#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("copy rows begin!\n");
	}
#endif

	int base_global_row = (tile_y  * blockIdx.y ); 
	int base_global_col = blockDim.x*blockIdx.x;
	int base_global_idx = N*base_global_row + base_global_col ;
	int nextRow = base_global_row+1;
	bool legalNextRow = (nextRow<M)?1:0;
	int t = threadIdx.x;
	bool legalCurCol = (base_global_col + t)<N;
	int idx = (base_global_row/tile_y)*2*N + t+base_global_col;
	int idx_nextrow = idx + N;
	if(legalCurCol){
		shared_rows[idx] = dst[base_global_idx + t];
	}
	if(legalNextRow&&legalCurCol){
		shared_rows[idx_nextrow] = dst[base_global_idx + N+t];
	}
	__syncthreads();


#ifdef CUDA_CUDA_DEBUG
//	if(threadIdx.x==0){
//		printf("blockIdx.x = %d,blockIdx.y = %d\n",blockIdx.x,blockIdx.y);
//	}
//	if(blockIdx.y==1 && threadIdx.x==0){
//		printf("addr: %d\n",idx_nextrow);
//	}
	if(blockIdx.y==0 && blockIdx.x==2 && (t==0 || t==1)){	
		printf("addr:%d, val = %f\n", idx_nextrow,shared_rows[idx_nextrow]);
	}
#endif

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("copy rows finish!\n");
	}
#endif
}

__global__ void gpu_stencil2D_4pt_hack5_cp_cols(double * dst, double * shared_cols, double *shared_rows,int tile_x,int tile_y, int M, int N){

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.y==0)){
		printf("copy cols begin!\n");
	}
#endif

	int base_global_row = tile_y  * blockIdx.y; 
	int base_global_col = tile_x  * blockIdx.x;
	int base_global_idx = N*base_global_row + base_global_col ;
	int nextCol = base_global_col+1;
	bool legalNextCol = (nextCol<N);
	int t = threadIdx.y;
	int idx = 2*M*blockIdx.x + t + base_global_row;
	int idx_nextCol = idx + M ;
	bool legalCurRow = (base_global_row + t)<M;
	if(legalCurRow){
		shared_cols[idx] = dst[base_global_idx + t*N];
	}
	if(legalNextCol && legalCurRow){
		shared_cols[idx_nextCol] = dst[base_global_idx + t*N+1];
	}
	__syncthreads();


#ifdef CUDA_CUDA_DEBUG
//	if(threadIdx.y==0){
//		printf("blockDimy = %d\n",blockDim.y);
//	}
	if(blockIdx.x==1 && t<5){
		printf("addr: %d ,%f,\n",idx_nextCol,shared_cols[idx_nextCol]);
	}
#endif

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.y==0)){
		printf("copy cols finish!\n");
	}
#endif
}

__global__ void gpu_stencil2D_4pt_hack5(double * dst, double * shared_cols, double *shared_rows,int tile_y,int M, int N)
{
#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("kernel begin!\n");
	}
#endif
	//Declaring the shared memory array for source
	extern	__shared__ double shared_mem[] ;

	//indexes
	int i,j, k,curRow;
                           //Cols   *  numRows/Tile * tileIndex  
	int base_global_row = ( tile_y * blockIdx.y ); 
	int base_global_col = blockDim.x*blockIdx.x;
	int base_global_idx = N*base_global_row + base_global_col ;
	int center = 1,north = 0,south = 2; //indexes for the current location in the shared memory
	int t = threadIdx.x;

	//copy the shared memory to fill the pipeline
	bool legalCol = (base_global_col + t )<N;
	bool legalCol1 = (base_global_col + t +1)<N;
	bool legalCol2= (base_global_col+t+2)<N;
	bool legalColn= (base_global_col+t+blockDim.x)<N;

	shared_mem [t] = shared_rows[base_global_col + t + blockIdx.y * N*2];
	if(t==2 || t==3){
		shared_mem [blockDim.x + t-2] = shared_rows[base_global_col+blockIdx.y*N*2+blockDim.x  + t-2];
	}
#ifdef CUDA_CUDA_DEBUG
//	if(blockIdx.y==0 && blockIdx.x==1 && (t==2||t==3)){
//		printf("addr: %d,val: %f\n",blockDim.x+t-2,shared_mem[blockDim.x+t-2]);
//	}
#endif
	for (i = 1 ; i < 1+HALO*2 ; i ++ ){
		k = base_global_idx+i*N+t;
		j = i*(blockDim.x+2) + t;
		bool legalRow = (base_global_row+i)<M;
		shared_mem [j+1] =legalRow?( legalCol1?dst[k+1]:0.0):0.0;
	
		if((t==1)&&legalColn&&legalRow){
			shared_mem[j+blockDim.x]=(blockIdx.x == (gridDim.x-1))?dst[k+blockDim.x]:shared_cols[blockIdx.x*2*M+3*M+i+base_global_row];
		}
		if(t==0){
			shared_mem[j] = shared_cols[blockIdx.x*2*M+base_global_row+i];
		}
	}
		
	__syncthreads();


#ifdef CUDA_CUDA_DEBUG
	if(blockIdx.y==0 && blockIdx.x==1 && (t==1||t==0)){
		printf("addr: %d,val: %f\n",blockDim.x+(blockDim.x+2)+t,shared_mem[blockDim.x + (blockDim.x+2)+t]);

		printf("addr: %d,val: %f\n",2*(blockDim.x+2)+blockDim.x+t,shared_mem[2*(blockDim.x+2) + blockDim.x+t]);
	}
#endif

	//Pipelined copy one row and process it
	for ( curRow = HALO; curRow < tile_y; curRow+=1 )
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

		shared_mem [north*(blockDim.x+2)+t+1] =(legalRow2&&legalCol1)?dst[k+1]:0.0;

		if((t==1)&&legalColn&&legalRow2){
			shared_mem[north*(blockDim.x+2)+t+blockDim.x]=(blockIdx.x == (gridDim.x-1))?dst[k+blockDim.x]:shared_cols[blockIdx.x*2*M+3*M+nextRow2];
		}	
		if((t==0)&&legalRow2){
			shared_mem[north*(blockDim.x+2)+t] = shared_cols[nextRow2+blockIdx.x*2*M];
		}
		
		center = ROTATE_UP(center,3);
		south  = ROTATE_UP(south,3);
		north  = ROTATE_UP(north,3);
		__syncthreads();
	}
	int lastRow1 = base_global_row+curRow+1;
	bool legalLastRow1 = (lastRow1)<M;
	
	if(legalLastRow1){
		
		shared_mem[south*(blockDim.x+2)+t] = shared_rows[base_global_col + t + blockIdx.y * N*2+N*3] ;

		if(t==2 || t==3){
			shared_mem [south*(blockDim.x+2) +blockDim.x + t-2] = shared_rows[base_global_col+blockIdx.y*N*2+3*N+blockDim.x  + t-2];
		}
#ifdef CUDA_CUDA_DEBUG
//		if((blockIdx.x==1)&&((t==2)||(t==3))&&(blockIdx.y==0)){ 
//			printf("addr: %d,val: %f\n",base_global_col+blockIdx.y*N*2+3*N+blockDim.x  + t-2 ,shared_mem[south*(blockDim.x+2)+blockDim.x+t-2]);
//		}
#endif
		__syncthreads();
	}
	if((legalCol2)&& legalLastRow1){
			dst[base_global_idx + curRow*N + j] =(shared_mem[north*(blockDim.x+2)+j] + shared_mem[south*(blockDim.x+2)+j] + shared_mem[center*(blockDim.x+2)+j-1] + shared_mem[center*(blockDim.x+2)+j+1] )/5.5;
	
	}


#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(threadIdx.x==0)){
		printf("kernel finish!\n");
	}
#endif
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

//extern "C"
//void
//stencil2D4pt_gpu( double * __restrict__ dst, double* __restrict__ src,
//		  const size_t M, const size_t N, 
//		  const size_t NUM_ITERATIONS)//M Rows by N Columns
//{
//		
//	double size = sizeof(double) * M * N;
//
//	//device memory allocation
//	double * d_dst, * d_src;
//	cudaMalloc( (void **) &d_dst, size);
//	cudaMalloc( (void **) &d_src, size);
//	
//	//dimmensions for indexes
//	// TODO the -2 is to remove the borders
//	dim3 dimBlock(MAX_BLOCK_DIM,MAX_BLOCK_DIM);
//	int gridx = (N-2)/(MAX_BLOCK_DIM*TILE_SIZE) + (((N-2)%(MAX_BLOCK_DIM*TILE_SIZE) == 0)? 0:1 ) ;
//	int gridy = (M-2)/(MAX_BLOCK_DIM*TILE_SIZE) + (((M-2)%(MAX_BLOCK_DIM*TILE_SIZE) == 0)? 0:1 ) ;
//	dim3 dimGrid(gridx,gridy);
//
//	//Shared memory size = inside + halo
//	int shMemSize=MAX_BLOCK_DIM*TILE_SIZE*MAX_BLOCK_DIM*TILE_SIZE*sizeof(double)+(HALO*MAX_BLOCK_DIM*TILE_SIZE+HALO*HALO)*4*sizeof(double);
//	
//	//Hackaton dimensions
//	dim3 dimGrid_hack1((N-HALO*2)/GRID_TILE_X,(M-HALO*2)/GRID_TILE_Y);
//
//	//Copying the device memory
//	cudaMemcpy(d_src, src, size, cudaMemcpyHostToDevice);
//	cudaMemcpy(d_dst, dst, size, cudaMemcpyHostToDevice);
//
//	//printf("CUDA Stencil Code running... cycles = %d. dim = %d by %d \n",NUM_ITERATIONS,M,N);
//	
//    int time_step = NUM_ITERATIONS;
//
//    while (time_step-- > 0) 
//    {
//    	//gpu_stencil2D_4pt<<<dimGrid,dimBlock,shMemSize>>>(d_dst,d_src,M,N);
//		//gpu_stencil2D_4pt_hack1<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N); //JOSE Hackathon!
//		//printf("before: d_src[10] = %ld",d_src[10]);
//
//		gpu_stencil2D_4pt_hack2<<<dimGrid_hack1,NUM_THREADS>>>(d_dst,d_src,M,N);
//		//Inline swapping.
//		
//		//printf("after: d_src[10] = %ld",d_src[10]);
//		double * temp;
//		if ( NUM_ITERATIONS%2 ==0 || time_step !=0)
//		{
//			temp=d_src;
//			d_src=d_dst;
//			d_dst=temp;
//		}
//	}
//	
//	
//	//Copying memory back from device to DRAM
//	//cudaMemcpy(src, d_src, size, cudaMemcpyDeviceToHost);
//	cudaMemcpy(dst, d_dst, size, cudaMemcpyDeviceToHost);
//	cudaMemcpy(src, d_src, size, cudaMemcpyDeviceToHost);
//	
//	//Free device memory
//	cudaFree(d_src); cudaFree(d_dst);
//}




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
		printf("sharedMemSize: %d B, total sharedMemSize: %d B\n",sharedMemSize, sharedMemSize*dimGrid.x*dimGrid.y);
#endif
		gpu_stencil2D_4pt_hack4<<<dimGrid,dimBlock,sharedMemSize>>>(d_dst,d_src,M,N);
#ifdef CUDA_DARTS_DEBUG
		printf("gpu kernel return to host, but kernel haven't finished!\n");
#endif

}
void gpu_kernel5(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N){
		int sharedMemSize = sizeof(double)*(1+HALO*2)*(dimBlock.x+2);
#ifdef CUDA_DARTS_DEBUG
		printf("sharedMemSize: %d B, total sharedMemSize: %d B\n",sharedMemSize, sharedMemSize*dimGrid.x*dimGrid.y);
#endif
		gpu_stencil2D_4pt_hack5<<<dimGrid,dimBlock,sharedMemSize>>>(d_dst,sharedCols,sharedRows,tile_y,M,N);
#ifdef CUDA_DARTS_DEBUG
		printf("gpu kernel return to host, but kernel haven't finished!\n");
#endif

}


void gpu_kernel5_cp_rows(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N){

		gpu_stencil2D_4pt_hack5_cp_rows<<<dimGrid,dimBlock>>>(d_dst,sharedCols,sharedRows,tile_y,M,N);
}

void gpu_kernel5_cp_cols(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_x,int tile_y,int M, int N){

		gpu_stencil2D_4pt_hack5_cp_cols<<<dimGrid,dimBlock>>>(d_dst,sharedCols,sharedRows,tile_x,tile_y,M,N);
}


void gpu_kernel5_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N){
		int sharedMemSize = sizeof(double)*(1+HALO*2)*(dimBlock.x+2);
#ifdef CUDA_DARTS_DEBUG
		printf("Kernel5 stream: sharedMemSize: %d B, total sharedMemSize: %d B\n",sharedMemSize, sharedMemSize*dimGrid.x*dimGrid.y);
#endif

        gpu_stencil2D_4pt_hack5<<<dimGrid,dimBlock,sharedMemSize,stream>>>(d_dst,sharedCols,sharedRows,tile_y,M,N);
#ifdef CUDA_DARTS_DEBUG
		printf("gpu kernel return to host, but kernel haven't finished!\n");
#endif


}


void gpu_kernel5_stream_cp_rows(cudaStream_t &stream ,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N){

#ifdef CUDA_DARTS_DEBUG
		printf("Kernel5 stream: copy Rows. \n");
#endif
        gpu_stencil2D_4pt_hack5_cp_rows<<<dimGrid,dimBlock,0,stream>>>(d_dst,sharedCols,sharedRows,tile_y,M,N);

}

void gpu_kernel5_stream_cp_cols(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_x,int tile_y,int M, int N){

#ifdef CUDA_DARTS_DEBUG
		printf("Kernel5 stream: copy cols. \n");
#endif
        gpu_stencil2D_4pt_hack5_cp_cols<<<dimGrid,dimBlock,0,stream>>>(d_dst,sharedCols,sharedRows,tile_x,tile_y,M,N);

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
bool checkGpu(cudaStream_t *stream, size_t n){
    for (size_t i=0;i<n;++i){
  
#ifdef CUDA_DARTS_DEBUG
		printf("checkGpu: %d \n",i);
#endif
        if (cudaSuccess != cudaStreamQuery(stream[i]))
            return false;
    }
    return true;
}


extern "C"
void
stencil2D4pt_gpu( double * __restrict__ h_dst, double* __restrict__ h_src, const size_t nRows, const size_t nCols, const size_t timestep)//M Rows by N Columns
{


	double *d_dst ;
	double *d_sharedCols ;
	double *d_sharedRows ;
	double d_size = sizeof(double) * nRows * nCols;
	int64_t d_size_sharedCols ;
	int64_t d_size_sharedRows ;
	
	size_t gpuMemMax = 0;
	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
    gpuMemMax =(2*GB)> gpu_mem_valid_t?gpu_mem_avail_t: 2*GB;

    int tile_y = GRID_TILE_Y;
    int tile_x = NUM_THREADS;
   
	d_size_sharedCols = sizeof(double)*nRows* (std::ceil(1.0*nCols/NUM_THREADS)) *2;
	d_size_sharedRows = sizeof(double)*nCols*(std::ceil(1.0*nRows/tile_y))*2;
	double req_size = sizeof(double)* nRows*nCols +  d_size_sharedCols + d_size_sharedRows  ;

	int nRowsGpu;
	int nRowsGpuMax;
    
	int nGPU = 1;
	int gpuPos = 0;

	cudaError err1,err2,err3,err4,err5;
	
	if(req_size<gpuMemMax){
		nGPU = 1;
		nRowsGpu = nRows;

		int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
		int blockDimy = 1;
		int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
		int gridDimy = std::ceil(1.0*nRowsGpu/tile_y); //GRID_TILE_Y=10, it needs to change.
		
		dim3 dimGrid(gridDimx,gridDimy);
		dim3 dimBlock(blockDimx,blockDimy);
		d_size_sharedCols = sizeof(double)*nRowsGpu*gridDimx*2;
		d_size_sharedRows = sizeof(double)*nCols*gridDimy*2;

		err1 = cudaMalloc( (void **) &d_dst, d_size);
		err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
		err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);

#ifdef CUDA_ERROR_CHECKING
		if(err1!=cudaSuccess){
			printf("GpuKernelWithAllTimeSteps: cuda malloc1: %s \n",cudaGetErrorString(err1));
			exit(-1);
		}
		if(err2!=cudaSuccess){
			printf("GpuKernelWithAllTimeSteps: cuda malloc2: %s \n",cudaGetErrorString(err2));
			exit(-1);
		}

		if(err3!=cudaSuccess){
			printf("GpuKernelWithAllTimeSteps: cuda malloc3: %s \n ",cudaGetErrorString(err3));
			exit(-1);
		}

#endif
		size_t pos1 = gpuPos*nCols;	

		err4 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);

#ifdef CUDA_ERROR_CHECKING
		if(err4!=cudaSuccess){
		
			printf("GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst: %s \n ",cudaGetErrorString(err4));
			exit(-1);
		}
#endif

		int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
		int blockDimy_rows = 1;
		int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
		int gridDimy_rows = std::ceil(1.0*nRowsGpu/tile_y);
	
		int blockDimx_cols = 1 ;
		int blockDimy_cols = (nRowsGpu>NUM_THREADS)?NUM_THREADS:nRows;
		int gridDimx_cols = gridDimx;
		int gridDimy_cols = std::ceil(1.0*nRowsGpu/blockDimy_cols);


		dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
		dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
	
		dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
		dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);

		size_t ts = timestep; 
	
		while(ts-- >0){
			gpu_kernel5_cp_rows(dimGrid_rows,dimBlock_rows,d_dst, d_sharedCols, d_sharedRows, tile_y,nRowsGpu, nCols);
			gpu_kernel5_cp_cols(dimGrid_cols,dimBlock_cols,d_dst, d_sharedCols, d_sharedRows, tile_x,tile_x,nRowsGpu, nCols);
			gpu_kernel5(dimGrid,dimBlock,d_dst,d_sharedCols,d_sharedRows,tile_y,nRowsGpu,nCols);
		}
		
		err5 = cudaDeviceSynchronize();

#ifdef CUDA_ERROR_CHECKING
		if(err5!=cudaSuccess){
			printf("GpuKernelWithAllTimeSteps: cuda deviceSynchronize:  %s \n ",cudaGetErrorString(err5));
			exit(-1);
		}
#endif
	
#ifdef VERIFICATION
        if(timestep%2==0){
			SWAP_PTR(&h_dst ,&h_src);
        }
#endif
	
        err1=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
		if(err1!=cudaSuccess){
			printf(" GpuKernelWithAllTimeSteps: cuda memcpyDeviceToHost:  %s \n ",cudaGetErrorString(err1));
			
			exit(-1);
		}
#endif
	
	    err1 = cudaFree(d_dst);
		err2 = cudaFree(d_sharedRows);
	    err3 = cudaFree(d_sharedCols);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_dst :  %s \n ",cudaGetErrorString(err1));
			exit(-1);
	    }
	
	    if(err2!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedRows :  %s \n ",cudaGetErrorString(err2));
			exit(-1);
	    }
	
	    if(err3!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedCols:  %s \n ",cudaGetErrorString(err3));
			exit(-1);
	    }
#endif
	
	
	}else{
		nGPU = std::ceil(req_size/gpuMemMax); 
		nRowsGpu = nRows;
		int nStream = 4 ;
		cudaStream_t *stream ;
		stream = new cudaStream_t[nStream];
		for(int i=0;i<nStream;++i){
			cudaStreamCreate(&stream[i]);
		}

		int vnStream = nStream*nGPU;
		int nTile_y = nRows/(tile_y * vnStream);
		
		int chunk = nTile_y*tile_y;
		int chunk2= chunk+2;
		
		int nRowsGpuBlock = nStream*chunk2 + nRows-nGPU*nStream*chunk;
		int64_t nRowsGpuStream;
		int64_t d_size_stream;
		
		int blockDimx =( (nCols-2)>NUM_THREADS)?NUM_THREADS:(nCols-2);
		int blockDimy = 1;
		int gridDimx = std::ceil(1.0*(nCols-2)/blockDimx);
		int gridDimy = std::ceil(1.0*nRowsGpuBlock/tile_y); 
		
		dim3 dimBlock(blockDimx,blockDimy);
		//dim3 dimGrid(gridDimx,gridDimy);
		
		d_size = sizeof(double)*nRowsGpuBlock*nCols;
		d_size_sharedCols = sizeof(double) * nRowsGpuBlock*gridDimx*2 ;
		d_size_sharedRows = sizeof(double) * nCols* gridDimy*2;

		err1 = cudaMalloc( (void **) &d_dst, d_size);
		err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
		err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);


#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	        
			printf("GpuKernelPureGpuWithStreams: cuda malloc d_dst:  %s \n ",cudaGetErrorString(err1));
			exit(-1);
	    }
	    if(err2!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda mallock d_sharedRows :  %s \n ",cudaGetErrorString(err2));
			exit(-1);
	    }
	    
	    if(err3!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda mallock d_sharedCols :  %s \n ",cudaGetErrorString(err3));
			exit(-1);
	    }
#endif
    
	    int blockDimx_rows =( nCols>NUM_THREADS)?NUM_THREADS:nCols;
	    int blockDimy_rows = 1;
	    int gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	    int gridDimy_rows; 
	    
	    int blockDimx_cols = 1 ;
	    int blockDimy_cols;
	    int gridDimx_cols = gridDimx;  
	    int gridDimy_cols; 
	    
	    dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows);
	    //dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
	
	    int h_pos;
	    int d_pos;
	    int pos0 = gpuPos*nCols;
		
		size_t ts = timestep;

	    while(ts-- >0){
	        for (size_t i = 0; i<nGPU; ++i){
	            for (size_t j =0; j<nStream;++j){
	                int ps = i*nStream+j;
	                nRowsGpuStream = ((i==(nGPU-1))&&(j==(nStream-1)))? (nRows-ps*chunk) :chunk2;
	                h_pos = pos0+ps*chunk*nCols;
	                d_pos = j*chunk2*nCols; 
	
#ifdef CUDA_ERROR_CHECKING
	                err3 = cudaGetLastError();
	                if(cudaSuccess != err3){
						printf("GpuKernelPureGpuWithStreams multiple streams: kernel5 stream error :  %s \n ",cudaGetErrorString(err3));
						exit(-1);
	                }
#endif
	
	                d_size_stream = sizeof(double)*nCols*nRowsGpuStream;
	                err1 = cudaMemcpyAsync(d_dst+d_pos, h_src+h_pos, d_size_stream, cudaMemcpyHostToDevice,stream[j]);
#ifdef CUDA_ERROR_CHECKING
	                if(err1!=cudaSuccess){
						printf("GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device :  %s \n ",cudaGetErrorString(err1));
						
						exit(-1);
	                }
#endif
	            
	                gridDimy_rows = std::ceil(1.0*nRowsGpuStream/tile_y);
	                dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows);
	                gpu_kernel5_stream_cp_rows(stream[j],dimGrid_rows,dimBlock_rows,d_dst+d_pos , d_sharedCols, d_sharedRows+j*nTile_y*2*nCols, tile_y,nRowsGpuStream, nCols);
	            
#ifdef CUDA_ERROR_CHECKING
	                err3 = cudaGetLastError();
	                if(cudaSuccess != err3){
						printf("GpuKernelWithStream multiple streams: kernel5 cuda cp rows :  %s \n ",cudaGetErrorString(err3));
						exit(-1);
	                }
#endif
	            
	                blockDimy_cols = (nRowsGpuStream>NUM_THREADS)?NUM_THREADS:nRowsGpuStream;
	                gridDimy_cols = std::ceil(1.0*nRowsGpuStream/blockDimy_cols);
	        	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols);
	                dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols);
	                int addrCol = j*chunk2*2*gridDimx_cols; 
	                gpu_kernel5_stream_cp_cols(stream[j],dimGrid_cols,dimBlock_cols,d_dst+d_pos, d_sharedCols+addrCol, d_sharedRows, tile_x,tile_x,nRowsGpuStream, nCols);
	
#ifdef CUDA_ERROR_CHECKING
	                err3 = cudaGetLastError();
	                if(cudaSuccess != err3){
						printf("GpuKernelWithStream multiple streams: kernel5 cuda cp cols :  %s \n ",cudaGetErrorString(err3));
						exit(-1);
	                }
#endif
	            
		            int gridDimy_stream = std::ceil(1.0*nRowsGpuStream/tile_y);
		            dim3 dimGrid_stream(gridDimx,gridDimy_stream);
	                gpu_kernel5_stream(stream[j] ,dimGrid_stream,dimBlock,d_dst+d_pos,d_sharedCols+addrCol,d_sharedRows+j*nTile_y*2*nCols,tile_y,nRowsGpuStream,nCols);
	            
#ifdef CUDA_ERROR_CHECKING
	                err3 = cudaGetLastError();
	                if(cudaSuccess != err3){
						printf("GpuKernelWithStream multiple streams: kernel5 cuda computation :  %s \n ",cudaGetErrorString(err3));
						exit(-1);
	                }
#endif
	            
		            err3=cudaMemcpyAsync(h_dst+h_pos+nCols, d_dst+d_pos+nCols,d_size_stream-(nCols)*2*sizeof(double), cudaMemcpyDeviceToHost,stream[j]);
	
#ifdef CUDA_ERROR_CHECKING
	                err3 = cudaGetLastError();
	                if(cudaSuccess != err3){
						printf("GpuKernelWithStream multiple streams: kernel5 Asyn Memory copy from device to host :  %s \n ",cudaGetErrorString(err3));
						exit(-1);
	                }
#endif
	            }
	        }
	   
		    err4 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
		    if(err4!=cudaSuccess){
				printf("GpuKernelPureGpuWithStreams: cuda deviceSynchronize: %s \n ",cudaGetErrorString(err4));
				exit(-1);
		    }
#endif
		    SWAP_PTR(&h_dst ,&h_src);
	    }
	
		err4 = cudaDeviceSynchronize();
	    err1 = cudaFree(d_dst);
		err2 = cudaFree(d_sharedRows);
	    err3 = cudaFree(d_sharedCols);

#ifdef CUDA_ERROR_CHECKING
		if(err4!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda deviceSynchronize :  %s \n ",cudaGetErrorString(err4));
			exit(-1);
		}
	
	    if(err1!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_dst :  %s \n ",cudaGetErrorString(err1));
			exit(-1);
	    }
	
	    if(err2!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedRows :  %s \n ",cudaGetErrorString(err2));
			exit(-1);
	    }
	
	    if(err3!=cudaSuccess){
			printf("GpuKernelPureGpuWithStreams: cuda memcpy free d_sharedCols:  %s \n ",cudaGetErrorString(err3));
			exit(-1);
	    }
#endif
		if(timestep%2){
			SWAP_PTR(&h_dst ,&h_src);
		}

		for(int i=0;i<nStream;++i){
			cudaStreamDestroy(stream[i]);
		}
		delete [] stream;

	}

}

__global__ void gpu_stencil37_hack1_cp_slices(double * dst, double * shared_rows, double *shared_cols,double *shared_slices,int n_rows, int n_cols,int n_slices,int tile_x,int tile_y, int tile_z){

#ifdef CUDA_DARTS_DEBUG
    if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)){
		printf("copy slices begin!\n");
        printf("gridDim.x=%d,gridDim.y=%d,gridDim.z=%d\n",gridDim.x,gridDim.y,gridDim.z);
        printf("blockDim.x=%d,blockDim.y=%d,blockDim.z=%d\n",blockDim.x,blockDim.y,blockDim.z);
        printf("tile_x=%d,tile_y=%d,tile_z=%d\n",tile_x,tile_y,tile_z);
	}
#endif
    int base_global_slice = tile_z * blockIdx.z;
	int base_global_row   = tile_y * blockIdx.y;
	int base_global_col   = blockDim.x * blockIdx.x;

	int area = n_rows*n_cols;
    int base_global_idx = base_global_slice*area + base_global_row * n_cols + base_global_col;
    
    int nextSlice = base_global_slice+1;
    bool legalNextSlice = (nextSlice<n_slices);
	int tx = threadIdx.x;
	bool legalCurCol = (base_global_col + tx)<n_cols;
    
    for(int ty=0;ty<tile_y;++ty){ 
        bool legalCurRow = (base_global_row + ty)<n_rows;
        int idx = blockIdx.z*area*2 + (base_global_row+ty)*n_cols + base_global_col+tx ;
        int idx_dst = base_global_idx + ty*n_cols+tx;
    	if(legalCurCol&&legalCurRow){
    		shared_slices[idx] = dst[idx_dst];
    	}
    	if(legalNextSlice&&legalCurCol&&legalCurRow){
    		shared_slices[idx+area] = dst[idx_dst+area];
    	}

    }
    __syncthreads();

#ifdef CUDA_CUDA_DEBUG
	if(blockIdx.z ==0 && blockIdx.y==0 && blockIdx.x==1 ){
	//	printf("shared_slices: addr:%d, val = %f\n",n_cols*n_rows + threadIdx.x,shared_slices[n_cols*n_rows+threadIdx.x]);
	    if(threadIdx.x==0||threadIdx.x==1||threadIdx.x==2){
            int addr = n_cols*n_rows + blockDim.x*blockIdx.x+threadIdx.x;
            int addr1 = n_cols*n_rows + blockDim.x*blockIdx.x+threadIdx.x+n_cols;
            int addr2 = n_cols*n_rows + blockDim.x*blockIdx.x+threadIdx.x+n_cols*2;
	    	printf("blockIdx.x=%d, blockIdx.y=%d, blockIdx.z=%d,shared_slices: addr= %d, val= %f\n",blockIdx.x, blockIdx.y, blockIdx.z, addr,shared_slices[addr]);
	    	printf("blockIdx.x=%d, blockIdx.y=%d, blockIdx.z=%d,shared_slices: addr= %d, val= %f\n",blockIdx.x, blockIdx.y, blockIdx.z, addr1,shared_slices[addr1]);
	    	printf("blockIdx.x=%d, blockIdx.y=%d, blockIdx.z=%d,shared_slices: addr= %d, val= %f\n",blockIdx.x, blockIdx.y, blockIdx.z, addr2,shared_slices[addr2]);
        }
    }
#endif

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)){
		printf("copy slices end!\n");
	}
#endif
}


__global__ void gpu_stencil37_hack1_cp_rows(double * dst, double * shared_rows, double *shared_cols,double *shared_slices,int n_rows, int n_cols,int n_slices,int tile_x,int tile_y, int tile_z){

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)){
		printf("copy rows begin\n");
        printf("gridDim.x=%d,gridDim.y=%d,gridDim.z=%d\n",gridDim.x,gridDim.y,gridDim.z);
        printf("blockDim.x=%d,blockDim.y=%d,blockDim.z=%d\n",blockDim.x,blockDim.y,blockDim.z);
        printf("tile_x=%d,tile_y=%d,tile_z=%d\n",tile_x,tile_y,tile_z);
	}
#endif
    int base_global_slice = tile_z * blockIdx.z;
	int base_global_row   = tile_y  * blockIdx.y;
	int base_global_col   = blockDim.x*blockIdx.x;

	int dst_area = n_rows*n_cols;
    int s_area = gridDim.y*n_cols*2; 
    
    int base_global_idx = base_global_slice*dst_area + base_global_row * n_cols + base_global_col;
    
    int nextRow = base_global_row+1;
	bool legalNextRow = nextRow<n_rows;

    int tx = threadIdx.x;
	bool legalCurCol = (base_global_col + tx)<n_cols;
    
    for(int tz=0;tz<tile_z;++tz){ 
        bool legalCurSlice = (base_global_slice + tz)<n_slices;
        int idx_dst =base_global_idx + tz*dst_area+ tx  ;
        int idx = (base_global_slice+tz)*s_area + blockIdx.y*n_cols*2+blockIdx.x*blockDim.x+ tx  ;
        if(legalCurCol && legalCurSlice){
    		shared_rows[idx] = dst[idx_dst];
    	}
        if(legalCurCol && legalCurSlice && legalNextRow){
    		shared_rows[idx+n_cols] = dst[idx_dst+n_cols];
    	}


    }
    __syncthreads();

#ifdef CUDA_CUDA_DEBUG
	if(blockIdx.y==0 && blockIdx.x==0 &&blockIdx.z==0 ){
        if((threadIdx.x==0 || threadIdx.x==1 || threadIdx.x==2 ) && threadIdx.y==0){
            
            int addr0 = base_global_idx+0*dst_area+threadIdx.x;
            int addr  = base_global_slice+blockIdx.x*blockDim.x + threadIdx.x;
            int addr1 = s_area*(base_global_slice+1)+n_cols+blockIdx.x*blockDim.x+ threadIdx.x;
            int addr2 = s_area*(base_global_slice+2)+n_cols+blockIdx.x*blockDim.x+ threadIdx.x;
		    printf("blockIdx.x=%d, blockIdx.y=%d,blockIdx.z=%d,dst      : z:%d, addr:%d, val = %f\n",blockIdx.x, blockIdx.y,blockIdx.z,0,addr0,dst[addr0]);
		    printf("blockIdx.x=%d, blockIdx.y=%d,blockIdx.z=%d,shared_rows: z:%d, addr:%d, val = %f\n",blockIdx.x, blockIdx.y,blockIdx.z,0,addr,shared_rows[addr]);
		    printf("blockIdx.x=%d, blockIdx.y=%d,blockIdx.z=%d,shared_rows: z:%d, addr:%d, val = %f\n",blockIdx.x, blockIdx.y,blockIdx.z,1,addr1,shared_rows[addr1]);
		    printf("blockIdx.x=%d, blockIdx.y=%d,blockIdx.z=%d,shared_rows: z:%d, addr:%d, val = %f\n",blockIdx.x, blockIdx.y,blockIdx.z,2,addr2,shared_rows[addr2]);
        }
        if(threadIdx.x==0 && threadIdx.y==0){
            int addr =  2*s_area+n_cols+256;
            int addr1 = 2*dst_area+n_cols+256;
            printf("shared_rows: addr:%d, val:%f\n", addr, shared_rows[addr]);  
            printf("dst        : addr:%d, val:%f\n", addr1, dst[addr1]);  
        }
	}
#endif

#ifdef CUDA_DARTS_DEBUG
	
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)){
		printf("copy rows end!\n");
	}
#endif
}


__global__ void gpu_stencil37_hack1_cp_cols(double * dst, double * shared_rows, double *shared_cols,double *shared_slices,int n_rows, int n_cols,int n_slices,int tile_x,int tile_y, int tile_z){

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.y==0)&& threadIdx.x==0 && threadIdx.z==0){
		printf("copy cols begin\n");
        printf("gridDim.x=%d,gridDim.y=%d,gridDim.z=%d\n",gridDim.x,gridDim.y,gridDim.z);
        printf("blockDim.x=%d,blockDim.y=%d,blockDim.z=%d\n",blockDim.x,blockDim.y,blockDim.z);
        printf("tile_x=%d,tile_y=%d,tile_z=%d\n",tile_x,tile_y,tile_z);
	}
#endif
    int base_global_slice = tile_z * blockIdx.z;
	int base_global_row   = blockDim.y * blockIdx.y;
	int base_global_col   = tile_x * blockIdx.x;

	int area_dst = n_rows*n_cols;
    int area_shared = gridDim.x*n_rows*2; 
    
#ifdef CUDA_CUDA_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.y==0&&threadIdx.x==0&&threadIdx.z==0)){
        printf("area_shared=%d\n",area_shared);
	}
#endif
    int base_global_idx = base_global_slice*area_dst + base_global_row * n_cols + base_global_col;
    
    int nextCol= base_global_col+1;
	bool legalNextCol = (nextCol<n_cols)?1:0;
    
    int ty = threadIdx.y;
	bool legalCurRow = (base_global_row + ty)<n_rows;
    
    for(int tz=0;tz<tile_z;++tz){ 
        bool legalCurSlice = (base_global_slice + tz)<n_slices;
        int idx_dst =base_global_idx + tz*area_dst + ty*n_cols ;
        int idx = (base_global_slice+tz)*area_shared + blockIdx.x*2*n_rows+blockIdx.y*blockDim.y+ty; 

        if(legalCurRow && legalCurSlice){
    		shared_cols[idx] = dst[idx_dst];
    	}
        if(legalCurRow && legalCurSlice && legalNextCol){
    		shared_cols[idx + n_rows] = dst[idx_dst + 1];
        }

        __syncthreads();
    }
    __syncthreads();

#ifdef CUDA_CUDA_DEBUG
	if(blockIdx.z ==0 && blockIdx.y==0 && blockIdx.x==0 && (threadIdx.x==0)){
//		printf("shared_cols: addr:%d, val = %f\n", threadIdx.y,shared_cols[threadIdx.y]);
	}
#endif

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.y==0 && threadIdx.x==0 && threadIdx.z==0)){
		printf("copy cols end!\n");
	}
#endif
}


void gpu_kernel37_cp_slices(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy slices begin!\n");
//        printf("dimBlock.x: %d, dimBlock.y: %d,dimBlock.z: %d\n",dimBlock.x,dimBlock.y,dimBlock.z);
//        printf("dimGrid.x: %d, dimGrid.y: %d,dimGrid.z: %d\n",dimGrid.x,dimGrid.y,dimGrid.z);
//#endif
		gpu_stencil37_hack1_cp_slices<<<dimGrid,dimBlock>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy slices finish!\n");
//#endif
}

void gpu_kernel37_cp_rows(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows,double * sharedCols, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy rows begin!\n");
//#endif
		gpu_stencil37_hack1_cp_rows<<<dimGrid,dimBlock>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy rows finish!\n");
//#endif
}

void gpu_kernel37_cp_cols(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy cols begin!\n");
//#endif
		gpu_stencil37_hack1_cp_cols<<<dimGrid,dimBlock>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);
//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy cols finish!\n");
//#endif
}

void gpu_kernel37_cp_slices_stream(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy slices begin!\n");
//        printf("dimBlock.x: %d, dimBlock.y: %d,dimBlock.z: %d\n",dimBlock.x,dimBlock.y,dimBlock.z);
//        printf("dimGrid.x: %d, dimGrid.y: %d,dimGrid.z: %d\n",dimGrid.x,dimGrid.y,dimGrid.z);
//#endif
		gpu_stencil37_hack1_cp_slices<<<dimGrid,dimBlock,0,stream>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy slices finish!\n");
//#endif
}

void gpu_kernel37_cp_rows_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows,double * sharedCols, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy rows begin!\n");
//#endif
		gpu_stencil37_hack1_cp_rows<<<dimGrid,dimBlock,0,stream>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy rows finish!\n");
//#endif
}

void gpu_kernel37_cp_cols_stream(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z){

//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy cols begin!\n");
//#endif
		gpu_stencil37_hack1_cp_cols<<<dimGrid,dimBlock,0,stream>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);
//#ifdef CUDA_DARTS_DEBUG
//		printf("gpu_kernel37 copy cols finish!\n");
//#endif
}

__global__ void gpu_stencil37_hack2(double * dst, double * shared_rows, double * shared_cols, double * shared_slices,int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z){

#ifdef CUDA_DARTS_DEBUG
	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)&&(threadIdx.y==0)){
		printf("3D kernel begin!\n");
//        printf("blockIdx.x = %d,blockIdx.y = %d, blockIdx.z = %d\n", blockIdx.x,blockIdx.y,blockIdx.z);
//        printf("threadIdx.x = %d,threadIdx.y = %d, threadIdx.z = %d \n", threadIdx.x,threadIdx.y,threadIdx.z);
        printf("gridDim.x=%d,gridDim.y=%d,gridDim.z=%d\n",gridDim.x,gridDim.y,gridDim.z);
        printf("blockDim.x=%d,blockDim.y=%d,blockDim.z=%d\n",blockDim.x,blockDim.y,blockDim.z);
        printf("tile_x=%d,tile_y=%d,tile_z=%d\n",tile_x,tile_y,tile_z);
	}
#endif
    
    int base_global_slice = tile_z * blockIdx.z;
	int base_global_row   = tile_y * blockIdx.y;
	int base_global_col   = tile_x * blockIdx.x;
    int global_area = n_rows*n_cols;
   
    int base_global_idx = base_global_slice * global_area + base_global_row * n_cols + base_global_col;
    int num_rows = ((base_global_row + tile_y+2)<n_rows)?(tile_y+2):(n_rows-base_global_row);
    int num_cols = ((base_global_col + tile_x+2)<n_cols)?(tile_x+2):(n_cols-base_global_col);
    int num_slices = ((base_global_slice + tile_z+2)<n_slices)?(tile_z+2):(n_slices-base_global_slice);

    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int tx1 = threadIdx.x + 1;
    int ty1 = threadIdx.y + 1;
	bool legalCol    = (base_global_col + tx   )<n_cols;
	bool legalCol1   = (base_global_col + tx + 1)<n_cols;
	bool legalCol2   = (base_global_col + tx + 2)<n_cols;
    bool legalColN   = (base_global_col + blockDim.x     )<n_cols;
    bool legalColN1  = (base_global_col + blockDim.x + 1 )<n_cols;
    bool legalColN2  = (base_global_col + blockDim.x + 2 )<n_cols;
    bool legalColNX  = (base_global_col + blockDim.x + tx)<n_cols;
    
    bool legalSlice1  = (base_global_slice + 1       )<n_slices;
    bool legalSlice2  = (base_global_slice + 2       )<n_slices;
    bool legalSliceN  = (base_global_slice + tile_z  )<n_slices;
    bool legalSliceN1 = (base_global_slice + tile_z+1)<n_slices;
   
    bool legalRow   = (base_global_row + ty    )<n_rows;
    bool legalRow1  = (base_global_row + ty + 1)<n_rows;
    bool legalRow2  = (base_global_row + ty + 2)<n_rows;
    bool legalRowN  = (base_global_row + blockDim.y     )<n_rows;
    bool legalRowN1 = (base_global_row + blockDim.y + 1 )<n_rows;
    bool legalRowN2 = (base_global_row + blockDim.y + 2 )<n_rows;
    bool legalRowNY = (base_global_row + blockDim.y + ty)<n_rows;
    
    //Declaring the shared memory array for source
	extern	__shared__ double shared_mem[] ;
    
    
    //====================================copy first 3 slices to shared_mem[]=================================//
    int s_stride_x = blockDim.x + 2;
    int s_stride_y = blockDim.y + 2;
    int shared_area = s_stride_x*s_stride_y; 
    int sslices_area = global_area;
    int slices_idx = (blockIdx.z)*global_area*2 +blockIdx.y*tile_y*n_cols + blockIdx.x*tile_x;
    //----copy first two slices from shared_slices----//
    //--copy x: 1~blockDim.x , y: 1~blockDimy
    shared_mem[(ty+1)*s_stride_x+tx+1] = (legalRow1&&legalCol1)?shared_slices[slices_idx+(ty+1)*n_cols+tx+1]:0;
    shared_mem[(ty+1)*s_stride_x+tx+1+shared_area] = (legalRow1&&legalCol1&&legalSlice1)?shared_slices[slices_idx+(ty+1)*n_cols+tx+1+sslices_area]:0;

    //--copy y=0,y=blockDim.y+1, x=1~blockDim.x --//
    if(ty ==0){
        shared_mem[ty*s_stride_x+ tx+1] = (legalRow&&legalCol1)?shared_slices[slices_idx+ty*n_cols+tx+1]:0;
        shared_mem[ty*s_stride_x+ tx+1+shared_area] = (legalRow&&legalCol1&&legalSlice1)?shared_slices[slices_idx+ty*n_cols+tx+1+sslices_area]:0;
    }
    if(ty==1){
        shared_mem[(blockDim.y+ty)*s_stride_x+ tx+1] = (legalRowNY&&legalCol1)?shared_slices[slices_idx+(ty+blockDim.y)*n_cols+tx+1]:0;
        shared_mem[(blockDim.y+ty)*s_stride_x+ tx+1+shared_area] = (legalRowNY&&legalCol1&&legalSlice1)?shared_slices[slices_idx+(ty+blockDim.y)*n_cols+tx+1+sslices_area]:0;
    }
    //--copy x= 0, x=blockDim.x+1, y=1~blockDim.x--//
    if(tx==0){
        shared_mem[(ty+1)*s_stride_x + tx] = (legalRow1)?shared_slices[slices_idx+(ty+1)*n_cols+tx]:0;
        shared_mem[(ty+1)*s_stride_x + tx+shared_area] = (legalRow1&&legalSlice1)?shared_slices[slices_idx+(ty+1)*n_cols+tx+sslices_area]:0;
    }
    if(tx==1){
        shared_mem[(ty+1)*s_stride_x + blockDim.x + tx] = (legalRow1&&legalColNX)?shared_slices[slices_idx+(ty+1)*n_cols+tx+blockDim.x]:0;
        shared_mem[(ty+1)*s_stride_x + blockDim.x + tx+shared_area] = (legalRow1&&legalColNX&&legalSlice1)?shared_slices[slices_idx+(ty+1)*n_cols+tx+blockDim.x+sslices_area]:0;
    }
    //----copy third plane from shared_rows, shared_cols, dst----//
    int srows_area = gridDim.y*n_cols*2;
    int scols_area = gridDim.x*n_rows*2;
    int dst_idx= base_global_idx + global_area*2;
    int s_idx  = shared_area*2;
    int srows_idx = srows_area*(base_global_slice+2) + blockIdx.y*n_cols*2 +blockIdx.x*tile_x;
    int scols_idx = scols_area*(base_global_slice+2) + blockIdx.x*n_rows*2 +blockIdx.y*tile_y;
    if(legalSlice2){
        //--copy x=1~blockDim.x, y=1~blockDim.y from dst to shared_mem--//
        if(ty>0){
            shared_mem[s_idx+(ty+1)*s_stride_x+tx+1] = (legalCol1&&legalRow1)?dst[dst_idx+(ty+1)*n_cols+tx+1]:0; 
        }
        //--copy y=0, y=blockDim.y+1 , x=1~blockDim.x from shared_rows to shared_mem--/
        if(ty==0){ //y=0, x=1~blockDim.x
            shared_mem[s_idx+tx+1] = (legalCol1)?shared_rows[srows_idx+tx+1]:0;
            shared_mem[s_idx+(ty+1)*s_stride_x+tx+1] = (legalCol1&&legalRow1)?shared_rows[srows_idx+(ty+1)*n_cols+tx+1]:0;
        }
        if(ty==1){//y=blockDim.y+1, x=1~blockDim.x
            shared_mem[s_idx+(ty+blockDim.y)*s_stride_x+tx+1] = (legalCol1&&legalRowN2)?shared_rows[srows_idx+(ty+2)*n_cols+tx+1]:((legalCol1&&legalRowNY)?dst[dst_idx+(ty+blockDim.y)*n_cols+tx+1]:0);
        }
        //--copy x=0, x=blockDim.x+1, y= 0~blockDim.y+1 from shared_cols to shared_mem--/
        //--[0,0],[0,blockDim.y+1],[blockDim.x+1,0],[blockDim.x+1, blockDim.y+1] never be used
        if(tx==0){  //x=0, y=1~blockDim.y   
            shared_mem[s_idx+(ty+1)*s_stride_x] = (legalRow1)?shared_cols[scols_idx+ty+1]:0;
        }
        if(tx==1){  //x=blockDim.x+1,y=1~blockDim.y
            shared_mem[s_idx+(ty+1)*s_stride_x+tx+blockDim.x] = (legalColN2&&legalRow1)?shared_cols[scols_idx+(tx+2)*n_rows+ty+1]:((legalColNX&&legalRow1)?dst[dst_idx+(ty+1)*n_cols+tx+blockDim.x]:0);
        }
    }
    
	__syncthreads();

#ifdef CUDA_CUDA_DEBUG
	if(blockIdx.z==0 && blockIdx.x==0 && blockIdx.y==0 ){
        if(threadIdx.y==0 || threadIdx.y==1 || threadIdx.y ==2){
            int s_s=2;
//            printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,0,threadIdx.y,threadIdx.x,shared_mem[0*shared_area+threadIdx.y*s_stride_x+threadIdx.x]);
//            printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,1,blockIdx.x,blockIdx.y, threadIdx.y,threadIdx.x,shared_mem[shared_area+threadIdx.y*s_stride_x+threadIdx.x]);
//            printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,2,blockIdx.x,blockIdx.y, threadIdx.y,threadIdx.x,shared_mem[2*shared_area+threadIdx.y*s_stride_x+threadIdx.x]);
        }
        if(threadIdx.x==0&&threadIdx.y==0){
            int addr= 2*shared_area+s_stride_x+1;
            int addr1=2*srows_area+n_cols+1;
            printf("shared_mem   : addr: %d,val: %f\n",addr,shared_mem[addr]);
            printf("shared_rows  : addr: %d,val: %f\n",addr1,shared_rows[addr1]);
        }
    }

	__syncthreads();
#endif
    //====================================copy first 3 slices to shared_mem[]=================================//

    //==============================compute plus copy 1 slices to shared_mem[]===============================//

    int center = 1;
    int north  = 0;
    int south  = 2;
    int curSlice;
    int lenSlice = (legalSliceN)? (tile_z):(n_slices-base_global_slice-1);
    for (curSlice = HALO; curSlice < lenSlice ; curSlice+=1){
        //----compute slice----//
        if(legalCol2 && legalRow2){
            dst[base_global_idx + curSlice*global_area + ty1*n_cols + tx1] = 
                 (   shared_mem[center*shared_area+ty1*s_stride_x+tx]  + shared_mem[center*shared_area+ty1*s_stride_x+tx+2]
                 +   shared_mem[center*shared_area+ty*s_stride_x+tx1]  + shared_mem[center*shared_area+(ty+2)*s_stride_x+tx1]
                 +   shared_mem[north*shared_area +ty1*s_stride_x+tx1] + shared_mem[south*shared_area+ty1*s_stride_x+tx1]
                 +   shared_mem[center*shared_area+ty1*s_stride_x+tx1] )/7.5 ; 
        }
		__syncthreads();

#ifdef CUDA_CUDA_DEBUG
    	if(blockIdx.z==0 && blockIdx.x==0 && blockIdx.y==0 ){
            if((threadIdx.y==0)&&(threadIdx.x==0)){
                printf("dst addr: %d\n", base_global_idx + curSlice*global_area + ty1*n_cols + tx1);
                //printf("curSlice: %d, lenSlice: %d\n", curSlice,lenSlice);
                printf("curSlice: %d, lenSlice: %d,gridDim.x=%d,gridDim.y=%d,gridDim.z=%d\n", curSlice,lenSlice,gridDim.x,gridDim.y,gridDim.z);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d, val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,center,ty1,tx,center*shared_area+ty1*s_stride_x+tx,shared_mem[center*shared_area+ty1*s_stride_x+tx]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,center,ty1,tx+2,center*shared_area+ty1*s_stride_x+tx+2,shared_mem[center*shared_area+ty1*s_stride_x+tx+2]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,center,ty,tx1,center*shared_area+ty*s_stride_x+tx1,shared_mem[center*shared_area+ty*s_stride_x+tx1]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,center,ty+2,tx1,center*shared_area+(ty+2)*s_stride_x+tx1,shared_mem[center*shared_area+(ty+2)*s_stride_x+tx1]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,north,ty1,tx1,north*shared_area+(ty1)*s_stride_x+tx1,shared_mem[north*shared_area+(ty1)*s_stride_x+tx1]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,south,ty1,tx1,south*shared_area+(ty1)*s_stride_x+tx1,shared_mem[south*shared_area+(ty1)*s_stride_x+tx1]);
                printf("blockIdx.x=%d,blockIdx.y=%d,blockIdx.z=%d, shared_mem   : addr: z:%d,y:%d,x=%d,addr: %d,val: %f\n",blockIdx.x,blockIdx.y,blockIdx.z,center,ty1,tx1,center*shared_area+(ty1)*s_stride_x+tx1,shared_mem[center*shared_area+(ty1)*s_stride_x+tx1]);
            }
//            if(threadIdx.x==0 && threadIdx.y==0){
//                int addr = 10*n_rows*n_cols + n_cols+1;
//                int addr1 = 10*srows_area + 2*n_cols+1;
//                printf("dst addr: %d,val: %f\n", addr, dst[addr] );
//                printf("shared_row addr: %d,val: %f\n", addr1, shared_rows[addr1] );
//            }
//                
//            if(legalCol1 && legalRow1){
//                printf("addr: %d\n",base_global_idx+curSlice*global_area+ty1*n_cols+tx1 );
//            }
        }
#endif
        //----copy next slice to shared_mem[]----//
        int ssSlice = curSlice+2;
        int g_ssSlice = base_global_slice+ssSlice;
        srows_idx = srows_area*g_ssSlice + blockIdx.y*n_cols*2 +blockIdx.x*tile_x;
        scols_idx = scols_area*g_ssSlice + blockIdx.x*n_rows*2 +blockIdx.y*tile_y;
        bool legalSliceSS = g_ssSlice<n_slices;
        dst_idx= base_global_idx + global_area*ssSlice;
        s_idx  = shared_area*north;
        if(legalSliceSS){
            //--copy x=1~blockDim.x, y=1~blockDim.y from dst to shared_mem--//
            if(ty>0){
                shared_mem[s_idx+(ty+1)*s_stride_x+tx+1] = (legalCol1&&legalRow1)?dst[dst_idx+(ty+1)*n_cols+tx+1]:0; 
            }
            //--copy y=0, y=blockDim.y+1 , x=1~blockDim.x from shared_rows to shared_mem--/
            if(ty==0){ //y=0, x=1~blockDim.x
                shared_mem[s_idx+tx+1] = (legalCol1)?shared_rows[srows_idx+tx+1]:0;
                shared_mem[s_idx+(ty+1)*s_stride_x+tx+1] = (legalCol1&&legalRow1)?shared_rows[srows_idx+(ty+1)*n_cols+tx+1]:0;
            }
            if(ty==1){//y=blockDim.y+1, x=1~blockDim.x
                shared_mem[s_idx+(ty+blockDim.y)*s_stride_x+tx+1] = (legalCol1&&legalRowN2)?shared_rows[srows_idx+(ty+2)*n_cols+tx+1]:((legalCol1&&legalRowNY)?dst[dst_idx+(ty+blockDim.y)*n_cols+tx+1]:0);
            }
            //--copy x=0, x=blockDim.x+1, y= 0~blockDim.y+1 from shared_cols to shared_mem--/
            //--[0,0],[0,blockDim.y+1],[blockDim.x+1,0],[blockDim.x+1, blockDim.y+1] never be used
            if(tx==0){  //x=0, y=1~blockDim.y   
                shared_mem[s_idx+(ty+1)*s_stride_x] = (legalRow1)?shared_cols[scols_idx+ty+1]:0;
            }
            if(tx==1){  //x=blockDim.x+1,y=1~blockDim.y
                shared_mem[s_idx+(ty+1)*s_stride_x+tx+blockDim.x] = (legalColN2&&legalRow1)?shared_cols[scols_idx+(tx+2)*n_rows+ty+1]:((legalColNX&&legalRow1)?dst[dst_idx+(ty+1)*n_cols+tx+blockDim.x]:0);
            }
        }
        
        center = ROTATE_UP(center,3);
		south  = ROTATE_UP(south,3);
		north  = ROTATE_UP(north,3);
		__syncthreads();
    }

    //==============================compute plus copy 1 slices to shared_mem[]===============================//

    //=========================copy plus compute last slice in one grid to shared_mem[]==========================//
    int g_nnSlice = base_global_slice+ tile_z + 1;
    bool legalSliceNN = g_nnSlice < n_slices;
    if(legalSliceNN){
        //----copy----//
        //--copy x: 1~blockDim.x , y: 1~blockDimy
        slices_idx = (blockIdx.z+1)*global_area*2 +global_area + blockIdx.y*tile_y*n_cols + blockIdx.x*tile_x;
        
        s_idx = south*shared_area;

        //--copy x: 1~blockDim.x , y: 1~blockDimy
        shared_mem[s_idx+(ty+1)*s_stride_x+tx+1] = (legalRow1&&legalCol1)?shared_slices[slices_idx+(ty+1)*n_cols+tx+1]:0;
        
        //--copy y=0,y=blockDim.y+1, x=1~blockDim.x --//
        if(ty ==0){
            shared_mem[s_idx+ty*s_stride_x+ tx+1] = (legalRow&&legalCol1)?shared_slices[slices_idx+ty*n_cols+tx+1]:0;
        }
        if(ty==1){
            shared_mem[s_idx+(blockDim.y+ty)*s_stride_x+ tx+1] = (legalRowNY&&legalCol1)?shared_slices[slices_idx+(ty+blockDim.y)*n_cols+tx+1]:0;
        }
        //--copy x= 0, x=blockDim.x+1, y=1~blockDim.x--//
        if(tx==0){
            shared_mem[s_idx+(ty+1)*s_stride_x + tx] = (legalRow1)?shared_slices[slices_idx+(ty+1)*n_cols+tx]:0;
        }
        if(tx==1){
            shared_mem[s_idx+(ty+1)*s_stride_x + blockDim.x + tx] = (legalRow1&&legalColNX)?shared_slices[slices_idx+(ty+1)*n_cols+tx+blockDim.x]:0;
        }
        __syncthreads();
        //----compute----//
        if(legalCol2 && legalRow2){
            dst[base_global_idx + curSlice*global_area + ty1*n_cols + tx1] = 
                 (   shared_mem[center*shared_area+ty1*s_stride_x+tx]  + shared_mem[center*shared_area+ty1*s_stride_x+tx+2]
                 +   shared_mem[center*shared_area+ty*s_stride_x+tx1]  + shared_mem[center*shared_area+(ty+2)*s_stride_x+tx1]
                 +   shared_mem[north*shared_area +ty1*s_stride_x+tx1] + shared_mem[south*shared_area+ty1*s_stride_x+tx1]
                 +   shared_mem[center*shared_area+ty1*s_stride_x+tx1] )/7.5 ; 
        }
    }

	__syncthreads();
    //=========================copy plus compute last slice in one grid to shared_mem[]==========================//
#ifdef CUDA_DARTS_DEBUG

	if((blockIdx.x==0)&&(blockIdx.y==0)&&(blockIdx.z==0)&&(threadIdx.x==0)&&(threadIdx.y==0)&&(threadIdx.z==0)){
		printf("3D kernel finish!\n");
	}
#endif

}


void gpu_kernel37(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices,int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z){
	int sharedMemSize = sizeof(double)*(1+HALO*2)*((tile_x+2)*(tile_y+2));
#ifdef CUDA_DARTS_DEBUG
		printf("sharedMemSize: %d B, total sharedMemSize: %d B\n",sharedMemSize, sharedMemSize*dimGrid.x*dimGrid.y*dimGrid.z);
		printf("gpu_kernel37: dimGrid.x= %d dimGrid.y= %d, dimGrid.z= %d\n",dimGrid.x,dimGrid.y,dimGrid.z);
#endif
		gpu_stencil37_hack2<<<dimGrid,dimBlock,sharedMemSize>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);
#ifdef CUDA_DARTS_DEBUG
		printf("gpu kernel37 return to host, but kernel haven't finished!\n");
#endif
}


void gpu_kernel37_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices,int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z){
	int sharedMemSize = sizeof(double)*(1+HALO*2)*((tile_x+2)*(tile_y+2));
#ifdef CUDA_DARTS_DEBUG
		printf("sharedMemSize: %d B, total sharedMemSize: %d B\n",sharedMemSize, sharedMemSize*dimGrid.x*dimGrid.y*dimGrid.z);
		printf("gpu_kernel37: dimGrid.x= %d dimGrid.y= %d, dimGrid.z= %d\n",dimGrid.x,dimGrid.y,dimGrid.z);
#endif
		gpu_stencil37_hack2<<<dimGrid,dimBlock,sharedMemSize,stream>>>(d_dst,sharedRows,sharedCols,sharedSlices,n_rows,n_cols,n_slices,tile_x,tile_y,tile_z);
#ifdef CUDA_DARTS_DEBUG
		printf("gpu kernel37 return to host, but kernel haven't finished!\n");
#endif
}


extern "C"
void
stencil3D7pt_gpu( double * __restrict__ h_dst, double* __restrict__ h_src, const size_t nRows, const size_t nCols, const size_t nSlices, const size_t timestep)
{

	double *d_dst = NULL;
    double *d_sharedRows;
    double *d_sharedCols;
    double *d_sharedSlices;
    size_t ts = timestep;

    int nGPU =0; 
    int gpuPos=0;
    
    int tWL=nSlices;
	int gpuWL=0;
    int gpuWLMin=3;
    int gpuWLMax=0;
    double req_size=0;

	int tile_x = 0;
	int tile_y = 0;
    int tile_z = 0;
    
    int blockDimx = 0;
    int blockDimy = 0;
    int blockDimz = 0;

    int gridDimx = 0;
    int gridDimy = 0;
    int gridDimz = 0;

    double d_size = 0;
    int64_t d_size_sharedCols = 0;
    int64_t d_size_sharedRows = 0;
    int64_t d_size_sharedSlices = 0;

    
    size_t gpuMemMax = 0;
    bool invokeStreams = false;
    int nStream = 4 ;
	cudaStream_t *stream ;

    int numThreads=0 ;

    int blockDimx_slices = 0;  
    int blockDimy_slices = 0; 
    int blockDimz_slices = 0; 
	int gridDimx_slices  = 0; 
	int gridDimy_slices  = 0; 
    int gridDimz_slices  = 0; 

    int blockDimx_rows   = 0; 
	int blockDimy_rows   = 0; 
    int blockDimz_rows   = 0; 
	int gridDimx_rows    = 0; 
	int gridDimy_rows    = 0; 
    int gridDimz_rows    = 0; 

	int blockDimx_cols   = 0; 
	int blockDimy_cols   = 0; 
	int blockDimz_cols   = 0; 
    int gridDimx_cols    = 0; 
	int gridDimy_cols    = 0; 
    int gridDimz_cols    = 0; 



	size_t gpu_mem_total_t = 0;
	size_t gpu_mem_avail_t = 0;
	size_t gpu_mem_valid_t = 0;
	cudaMemGetInfo(&gpu_mem_avail_t,&gpu_mem_total_t);
	gpu_mem_valid_t = gpu_mem_avail_t - XMB;
    
    gpuMemMax =(2*GB)> gpu_mem_valid_t?gpu_mem_avail_t: 2*GB;

	tile_x = (nCols>GRID_TILE37_X)?GRID_TILE37_X:nCols; //16
	tile_y = (nRows>GRID_TILE37_Y)?GRID_TILE37_Y:nRows; //16
    tile_z = (nSlices>GRID_TILE37_Z)?GRID_TILE37_Z:nSlices; //100 tile_z +2 < NUM_THREAD

    blockDimx = (nCols-2)> tile_x? tile_x:(nCols-2);
    blockDimy = (nRows-2)> tile_y? tile_y:(nRows-2);
    blockDimz = 1;

	gridDimx = std::ceil(1.0*(nCols)/blockDimx);
	gridDimy = std::ceil(1.0*(nRows)/blockDimy);
    gridDimz = std::ceil(1.0*(nSlices)/tile_z);
    d_size = sizeof(double)*nRows*nCols*nSlices;  
    d_size_sharedCols = sizeof(double)*nRows*nSlices*gridDimx*2;
    d_size_sharedRows = sizeof(double)*nCols*nSlices*gridDimy*2;
    d_size_sharedSlices = sizeof(double)*nRows*nCols*gridDimz*2;
    req_size = d_size + d_size_sharedCols + d_size_sharedRows + d_size_sharedSlices;

    gpuWLMax = std::floor(1.0*gpuMemMax/(sizeof(double)*(nRows*nCols+nRows*gridDimx*2+nCols*gridDimy*2+1.0*nRows*nCols*(1/tile_z))));
    
    cudaError err1,err2,err3,err4;
    
    numThreads= tile_x*tile_y;
    if(req_size < gpuMemMax){
        nGPU = 1;
        gpuWL = tWL;
        gpuPos = 0;


        blockDimx_slices = (nCols>numThreads)?numThreads:nCols; 
        blockDimy_slices = 1;
        blockDimz_slices = 1;
	    gridDimx_slices = std::ceil(1.0*nCols/blockDimx_slices);
	    gridDimy_slices = gridDimy;
        gridDimz_slices = gridDimz;

        blockDimx_rows = (nCols>numThreads)?numThreads:nCols;
	    blockDimy_rows = 1;
        blockDimz_rows = 1;
	    gridDimx_rows = std::ceil(1.0*nCols/blockDimx_rows);
	    gridDimy_rows = gridDimy;
        gridDimz_rows = gridDimz;

	    blockDimx_cols = 1 ;
	    blockDimy_cols = (nRows>numThreads)?numThreads:nRows;
	    blockDimz_cols = 1 ;
        gridDimx_cols = gridDimx;
	    gridDimy_cols = std::ceil(1.0*nRows/numThreads);
        gridDimz_cols = gridDimz;

	    err1 = cudaMalloc( (void **) &d_dst, d_size);
	    err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols);
	    err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows);
	    err4 = cudaMalloc( (void **) &d_sharedSlices, d_size_sharedSlices);
#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc1: %s \n ",cudaGetErrorString(err1));
	    	exit(-1);
	    }

	    if(err2!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc2: %s \n ",cudaGetErrorString(err2));
	    	exit(-1);
	    }

	    if(err3!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc3: %s \n ",cudaGetErrorString(err3));
	    	exit(-1);
	    }

	    if(err4!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc4: %s \n ",cudaGetErrorString(err4));
	    	exit(-1);
	    }

#endif
        int pos1 = gpuPos*nRows*nCols;
	    err1 = cudaMemcpy(d_dst, h_src+pos1, d_size, cudaMemcpyHostToDevice);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
		    printf("GpuKernelWithAllTimeSteps: cuda memcpyHostToDevice d_dst:%d \n ",cudaGetErrorString(err1));
		    exit(-1);
	    }
#endif


        dim3 dimGrid(gridDimx,gridDimy,gridDimz);
	    dim3 dimBlock(blockDimx,blockDimy,blockDimz);

	    dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
	    dim3 dimBlock_slices(blockDimx_slices,blockDimy_slices,blockDimz_slices);
        dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
	    dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows,blockDimz_rows);
	    dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);
	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols,blockDimz_cols);

	    while(ts-- >0){
            gpu_kernel37_cp_slices(dimGrid_slices,dimBlock_slices,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
            gpu_kernel37_cp_rows(dimGrid_rows,dimBlock_rows,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
            gpu_kernel37_cp_cols(dimGrid_cols,dimBlock_cols,d_dst,d_sharedRows, d_sharedCols, d_sharedSlices, nRows, nCols,nSlices,tile_x,tile_y,tile_z);
            gpu_kernel37(dimGrid,dimBlock,d_dst,d_sharedRows,d_sharedCols,d_sharedSlices,nRows,nCols,nSlices,tile_x,tile_y,tile_z);
	    }
	    
	    err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
		    printf("GpuKernelWithAllTimeSteps: DeviceSynchronize :%d \n ",cudaGetErrorString(err1));
		    exit(-1);
	    }
#endif

#ifdef VERIFICATION
        //if(ts%2==0){
    	//    SWAP_PTR(&h_dst ,&h_src);
        //}
#endif

        err1=cudaMemcpy(h_dst+pos1, d_dst,d_size, cudaMemcpyDeviceToHost);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
		    printf("GpuKernelWithAllTimeSteps: MemcpyDeviceToHost :%d \n ",cudaGetErrorString(err1));
		    exit(-1);
	    }
#endif

	    err1 = cudaFree(d_dst);
	    err2 = cudaFree(d_sharedCols);
	    err3 = cudaFree(d_sharedRows);
        err4 = cudaFree(d_sharedSlices);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free1: %s \n ",cudaGetErrorString(err1));
	    	exit(-1);
	    }

	    if(err2!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free2: %s \n ",cudaGetErrorString(err2));
	    	exit(-1);
	    }

	    if(err3!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free3: %s \n ",cudaGetErrorString(err3));
	    	exit(-1);
	    }

	    if(err4!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free4: %s \n ",cudaGetErrorString(err4));
	    	exit(-1);
	    }

#endif


    }else{
#ifdef CUDA_CUDA_DEBUG
		printf("cuda with streams\n");

#endif
        nGPU = std::ceil(req_size/gpuMemMax);
        gpuWL = tWL;
        gpuPos=0;
        invokeStreams = true;
        stream = new cudaStream_t[nStream];
        for(int i=0;i<nStream;++i){
            cudaStreamCreate(&stream[i]);
        }

        int vnStream = nStream*nGPU;

        int nTile_z = std::ceil(1.0*gpuWL/(tile_z*vnStream));
        int chunk = nTile_z*tile_z;
        int nSlicesChunk = nTile_z*tile_z + 2;
        int nSlicesChunkInit = nSlicesChunk;
        
        gridDimz = nTile_z;             //  chunk/tile_z
        
        int gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
	    d_size = sizeof(double)*nRows*nCols*nSlicesChunk;
		d_size_sharedRows   = sizeof(double)*nCols*gridDimy*2*nSlicesChunk;
        d_size_sharedCols   = sizeof(double)*nRows*gridDimx*2*nSlicesChunk;
        d_size_sharedSlices = sizeof(double)*nRows*nCols*2*gridDimz2;
       

	    double d_sizeInit = d_size;

	    err1 = cudaMalloc( (void **) &d_dst, d_size*nStream);
	    err2 = cudaMalloc( (void **) &d_sharedCols, d_size_sharedCols*nStream);
	    err3 = cudaMalloc( (void **) &d_sharedRows, d_size_sharedRows*nStream);
	    err4 = cudaMalloc( (void **) &d_sharedSlices, d_size_sharedSlices*nStream);
#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc1: %s \n ",cudaGetErrorString(err1));
	    	exit(-1);
	    }

	    if(err2!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc2: %s \n ",cudaGetErrorString(err2));
	    	exit(-1);
	    }

	    if(err3!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc3: %s \n ",cudaGetErrorString(err3));
	    	exit(-1);
	    }

	    if(err4!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: cuda malloc4: %s \n ",cudaGetErrorString(err4));
	    	exit(-1);
	    }

#endif
        
		blockDimx_slices = (nCols>numThreads)?numThreads:nCols; 
		blockDimy_slices = 1;
		blockDimz_slices = 1;
		gridDimx_slices = std::ceil(1.0*nCols/blockDimx_slices);
		gridDimy_slices = gridDimy;
		gridDimz_slices = gridDimz2;
		
		blockDimx_rows = (nCols>numThreads)?numThreads:nCols;
		blockDimy_rows = 1;
		blockDimz_rows = 1;
		gridDimx_rows = std::ceil(1.0*nCols/blockDimy_rows);
		gridDimy_rows = gridDimy;
		gridDimz_rows = gridDimz2;
		
		blockDimx_cols = 1 ;
		blockDimy_cols = (nRows>numThreads)?numThreads:nRows;
		blockDimz_cols = 1 ;
		gridDimx_cols = gridDimx;
		gridDimy_cols = std::ceil(1.0*nRows/blockDimy_cols);
		gridDimz_cols = gridDimz2;
        
	    dim3 dimBlock(blockDimx,blockDimy,blockDimz);
        
	    dim3 dimBlock_slices(blockDimx_slices,blockDimy_slices,blockDimz_slices);
        
	    dim3 dimBlock_rows(blockDimx_rows,blockDimy_rows,blockDimz_rows);
       
	    dim3 dimBlock_cols(blockDimx_cols,blockDimy_cols,blockDimz_cols);

        int pos1 = gpuPos*nRows*nCols;

#ifdef CUDA_CUDA_DEBUG
			printf("nGPU:%d\n",nGPU);
			printf("gridDimx=%d, gridDimy=%d, gridDimz=%d\n", gridDimx, gridDimy,gridDimz);
#endif
		while(ts-- >0){

            nSlicesChunk = nSlicesChunkInit;
            d_size = d_sizeInit; 
            gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
            gridDimz_slices = gridDimz2;
            gridDimz_rows = gridDimz2;
            gridDimz_cols = gridDimz2;
            gridDimz = gridDimz2;
			
			for (size_t i = 0; i<nGPU; ++i){
	            for (size_t j =0; j<nStream;++j){
    
                size_t ps = i*nStream+j;
                
                int64_t h_pos = ps*chunk*nRows*nCols;
                int64_t d_pos = j*nSlicesChunk*nRows*nCols; 
                int64_t s_srows_pos=j*nCols*gridDimy*2*nSlicesChunk;
                int64_t s_scols_pos=j*nRows*gridDimx*2*nSlicesChunk;
                int64_t s_sslices_pos=j*nRows*nCols*2*nTile_z;
                int64_t nSlicesLeft = gpuWL-ps*chunk;
                
#ifdef CUDA_CUDA_DEBUG
			printf("i:%d,j=%d,gpuWL=%d,chunk=%d,nSlicesLef=%d\n",i,j,gpuWL,chunk,nSlicesLeft);
#endif
                if(nSlicesLeft>1){
                
                    if(nSlicesLeft<=chunk+2){
                        d_size = sizeof(double)*nSlicesLeft*nRows*nCols;
                        nSlicesChunk = nSlicesLeft;
                        gridDimz2 = std::ceil(1.0*nSlicesChunk/tile_z);
                        gridDimz_slices = gridDimz2;
                        gridDimz_rows = gridDimz2;
                        gridDimz_cols = gridDimz2;
                        gridDimz = gridDimz2;
                    }
	                dim3 dimGrid_slices(gridDimx_slices,gridDimy_slices,gridDimz_slices);
                    dim3 dimGrid_rows(gridDimx_rows,gridDimy_rows,gridDimz_rows);
                    dim3 dimGrid_cols(gridDimx_cols,gridDimy_cols,gridDimz_cols);

                    dim3 dimGrid (gridDimx,gridDimy,gridDimz);
                
                    err1 = cudaMemcpyAsync(d_dst+d_pos, h_dst+h_pos, d_size, cudaMemcpyHostToDevice,stream[j]);
#ifdef CUDA_ERROR_CHECKING
                    if(err1!=cudaSuccess){
                        printf("GpuKernelWithStream multiple streams: cuda MemcpyAsync from host to device:%s \n ",cudaGetErrorString(err1));
                        exit(-1);
                    }
#endif
#ifdef CUDA_CUDA_DEBUG
					printf("kernel i=%d,j=%d!\n",i,j);
					printf("dimGrid_slices.x=%d,dimGrid_slices.y=%d,dimgrid_slices.z=%d\n",dimGrid_slices.x,dimGrid_slices.y,dimGrid_slices.z);
#endif
                    gpu_kernel37_cp_slices_stream(stream[j],dimGrid_slices,dimBlock_slices,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_cp_rows_stream(stream[j],dimGrid_rows,dimBlock_rows,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_cp_cols_stream(stream[j],dimGrid_cols,dimBlock_cols,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);
                    gpu_kernel37_stream(stream[j],dimGrid,dimBlock,d_dst+d_pos,d_sharedRows+s_srows_pos, d_sharedCols+s_scols_pos, d_sharedSlices+s_sslices_pos, nRows, nCols,nSlicesChunk,tile_x,tile_y,tile_z);

#ifdef CUDA_ERROR_CHECKING
                    err3 = cudaGetLastError();
                    if(cudaSuccess != err3){
                        printf("GpuKernelWithStream multiple streams: kernel5 cuda computation:%s \n ",cudaGetErrorString(err3));
                        exit(-1);
                    }
#endif
            
	                err3=cudaMemcpyAsync(h_dst+h_pos+nRows*nCols, d_dst+d_pos+nRows*nCols,d_size-sizeof(double)*nRows*nCols*2, cudaMemcpyDeviceToHost,stream[j]);

#ifdef CUDA_ERROR_CHECKING
                    err3 = cudaGetLastError();
                    if(cudaSuccess != err3){
                        printf("GpuKernelWithStream multiple streams: kernel5 Asyn Memory copy from device to host: %s \n",cudaGetErrorString(err3));
                        exit(-1);
                    }
#endif

                }
                
				}
            }
        }

	        err1 = cudaDeviceSynchronize();
#ifdef CUDA_ERROR_CHECKING
	        if(err1!=cudaSuccess){
		        printf("GpuKernelWithAllTimeSteps: DeviceSynchronize :%d \n ",cudaGetErrorString(err1));
		        exit(-1);
	        }
#endif



	    err1 = cudaFree(d_dst);
	    err2 = cudaFree(d_sharedCols);
	    err3 = cudaFree(d_sharedRows);
        err4 = cudaFree(d_sharedSlices);

#ifdef CUDA_ERROR_CHECKING
	    if(err1!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free1: %s \n ",cudaGetErrorString(err1));
	    	exit(-1);
	    }

	    if(err2!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free2: %s \n ",cudaGetErrorString(err2));
	    	exit(-1);
	    }

	    if(err3!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free3: %s \n ",cudaGetErrorString(err3));
	    	exit(-1);
	    }

	    if(err4!=cudaSuccess){
	    	printf("GpuKernelWithAllTimeSteps: free4: %s \n ",cudaGetErrorString(err4));
	    	exit(-1);
	    }

#endif


		for(int i=0;i<nStream;++i){
			cudaStreamDestroy(stream[i]);
		}
		delete [] stream;

    }

}
