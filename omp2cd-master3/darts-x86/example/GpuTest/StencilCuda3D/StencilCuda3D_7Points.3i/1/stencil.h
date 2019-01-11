#ifndef STENCIL_H_GUARD
#define STENCIL_H_GUARD

#define VERIFICATION 
//#define CUDA_DARTS_DEBUG
#define VERIFICATION_PRINT 
#define CUDA_ERROR_CHECKING
//#define CUDA_CUDA_DEBUG

#ifdef __cplusplus

#include <cuda.h>
#include <cuda_runtime.h>
#include<device_launch_parameters.h>
//#include<conio.h>
extern "C"
{
#endif


//#define KB 1024
//#define MB 1024*1024
//#define XMB 20*MB 
//#define GB 1024*MB

#define kb(x) (size_t (x)<<10)
#define mb(x) (size_t (x)<<20)
#define gb(x) (size_t (x)<<30)
#define KB kb(1)
#define MB mb(1)
#define XMB mb(20)
#define GB gb(1)

//#include "util.h"
#define SWAP(type,left,right) do { \
    type tmp = left;               \
    left     = right;              \
    right    = tmp;                \
} while(0)

inline void swap_ptr(void** left, void** right) {

	void* tmp = *left;
    *left     = *right;
    *right    = tmp;
}

#define SWAP_PTR(left,right) swap_ptr((void**)left,(void**)right)

//typedef struct stencil_arg_s {
//    double *dst, 
//		   *src;
//    size_t  n_rows, 
//			n_cols,
//			n_tsteps;
//} stencil_arg_t;
//
//typedef void (*stencil_funptr_t)(double*, double*, const size_t, const size_t, const size_t);
//typedef struct stencil_s {
//    stencil_arg_t *arg;
//    stencil_funptr_t stencil;
//} stencil_t;
//
//#define STENCIL_COMPUTE(c,d) (c)((d)->dst,(d)->src,(d)->n_rows,(d)->n_cols,(d)->n_tsteps)
//void* stencil_run(void* arg);
//
//void stencil_init        ( stencil_t    *stencil, stencil_funptr_t code, 
//                           const size_t  n_rows,  const size_t     n_cols, 
//                           size_t        n_tsteps );
//void stencil_arg_copy    ( stencil_arg_t* dest, stencil_arg_t* sourc );


void gpuAssert(cudaError_t code,const char *file,const int line,bool abort);   
#define gpuErrchk(ans){gpuAssert((ans),__FILE__,__LINE__,abort);}

//struct st_array3d_double{
//   int width;
//   int height;
//   int depth;
////   double src[height][width][depth];
//}


//Struct AA{
//    int a;
//};



void test_copy3d(dim3 dimGrid,dim3 dimBlock, double *dst,int tile_x,int tile_y,int tile_z );

void stencil2D4pt        ( double* __restrict__ dst,    double* __restrict__ src, 
                           const size_t     n_rows, const size_t     n_cols,
                           const size_t     n_tsteps );

void stencil2D4pt_gpu	 ( double * __restrict__ dst, double* __restrict__ src,
			   const size_t M, const size_t N, 
			   const size_t NUM_ITERATIONS);//M Rows by N Columns


void gpu_kernel1(dim3 dimGrid_hack1,double * d_dst, double * d_src, int M, int N);

void gpu_kernel3(cudaStream_t &stream,dim3 dimGrid_hack1,double * d_dst, double * d_src, int M, int N);

void gpu_kernel2(dim3 dimGrid_hack1,double *dst, double *src, double size, size_t ts, double * d_dst, double * d_src, int M, int N);


void gpu_kernel4(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * d_src, int M, int N);
void gpu_kernel5(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols,double *sharedRows,int tile_y, int M, int N);

void gpu_kernel5_cp_rows(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N);

void gpu_kernel5_cp_cols(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_x,int tile_y,int M, int N);

void gpu_kernel5_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N);

void gpu_kernel5_stream_cp_rows(cudaStream_t &stream ,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_y,int M, int N);

void gpu_kernel5_stream_cp_cols(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, int tile_x,int tile_y,int M, int N);

bool checkGpu(cudaStream_t *stream, size_t n);


void gpu_kernel37_cp_slices(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37_cp_rows(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);
void gpu_kernel37_cp_cols(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37(dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices,int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z);

void gpu_kernel37_cp_slices_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37_cp_rows_stream(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);
void gpu_kernel37_cp_cols_stream(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37_stream( cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices,int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z);

void gpu_kernel37_cp_slices_stream_p(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices, int d_xpitch,int d_ypitch,int d_zpitch,int s_xpitch,int s_ypitch,int s_zpitch,  int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37_cp_rows_stream_p(cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices,int d_xpitch,int d_zpitch,int s_xpitch,int s_ypitch,int s_zpitch, int d_ypitch, int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);
void gpu_kernel37_cp_cols_stream_p(cudaStream_t &stream, dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedCols, double * sharedRows, double * sharedSlices,int d_xpitch,int d_ypitch, int d_zpitch,int s_xpitch,int s_ypitch,int s_zpitch,int n_rows, int n_cols, int n_slices,int tile_x,int tile_y, int tile_z);

void gpu_kernel37_stream_p( cudaStream_t &stream,dim3 dimGrid,dim3 dimBlock,double * d_dst, double * sharedRows, double * sharedCols, double * sharedSlices,int d_xpitch, int d_ypitch,int d_zpitch,int s_xpitch,int s_ypitch,int s_zpitch, int n_rows,int n_cols, int n_slices,int tile_x, int tile_y, int tile_z);

#ifdef __cplusplus
}
#endif
#endif // STENCIL_H_GUARD
