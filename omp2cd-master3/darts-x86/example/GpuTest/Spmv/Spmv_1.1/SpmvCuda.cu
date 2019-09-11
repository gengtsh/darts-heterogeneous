#include "cudacommon.h"
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include "OptionParser.h"
#include "ResultDatabase.h"
#include "Spmv.h"
#include "util.h"
#include "conf.h"

using namespace std;


texture<float, 1> vecTex;  // vector textures
texture<int2, 1>  vecTexD;

// Texture Readers (used so kernels can be templated)
struct texReaderSP {
   __device__ __forceinline__ float operator()(const int idx) const
   {
       return tex1Dfetch(vecTex, idx);
   }
};

struct texReaderDP {
   __device__ __forceinline__ double operator()(const int idx) const
   {
       int2 v = tex1Dfetch(vecTexD, idx);
#if (__CUDA_ARCH__ < 130)
       // Devices before arch 130 don't support DP, and having the
       // __hiloint2double() intrinsic will cause compilation to fail.
       // This return statement added as a workaround -- it will compile,
       // but since the arch doesn't support DP, it will never be called
       return 0;
#else
       return __hiloint2double(v.y, v.x);
#endif
   }
};

template <typename floatType>
void memcpyHostToDevice(floatType *dst, floatType *src, int size ){
    CUDA_SAFE_CALL(cudaMemcpy(dst, src, size * sizeof(floatType),cudaMemcpyHostToDevice));
};

template <typename floatType>
void memcpyDeviceTexture(const void* devPtr, size_t size ){

    if (sizeof(floatType) == sizeof(float))
    {
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
        CUDA_SAFE_CALL(cudaBindTexture(0, vecTex, devPtr, channelDesc,size * sizeof(float)));
    }else {
        cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
        CUDA_SAFE_CALL(cudaBindTexture(0, vecTexD, devPtr, channelDesc,size * sizeof(int2)));
    }
};

template void memcpyHostToDevice<double>(double *dst, double *src, int size );
template void memcpyHostToDevice<int>(int *dst, int *src, int size );

template void memcpyDeviceTexture<double>(const void* devPtr, size_t size );
template void memcpyDeviceTexture<int>(const void* devPtr, size_t size );


// Forward declarations for kernels
template <typename fpType, typename texReader>
__global__ void
spmv_csr_scalar_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out);


template <typename fpType, typename texReader>
__global__ void
spmv_csr_scalar_section_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out,
                       const int secStart);

template <typename fpType, typename texReader>
__global__ void
spmv_csr_vector_kernel(const fpType * __restrict__ val,
             	       const int    * __restrict__ cols,
		               const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out);

template <typename fpType, typename texReader>
__global__ void
spmv_csr_vector_section_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out,
                       const int secStart);


template <typename fpType, typename texReader>
__global__ void
spmv_ellpackr_kernel(const fpType * __restrict__ val,
		             const int    * __restrict__ cols,
		             const int    * __restrict__ rowLengths,
                     const int dim, fpType * __restrict__ out);

template <typename fpType>
__global__ void
zero(fpType * __restrict__ a, const int size);


template <typename floatType>
void csrTestScalar(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

    int deviceStart = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+deviceStart;
    int secStart = h_rowDelimiters[0] ;
    floatType *h_val = csrHost->getVal()+secStart;
    int *h_cols = csrHost->getCols()+secStart;
    floatType *h_vec = csrHost->getVec()+deviceStart;
    floatType *h_out = csrHost->getOut()+deviceStart;

    int numRows = csrDevice->getNumRows();
    //int numNonZeroes = csrDevice->getNumNonZeroes();
    int numNonZeroes = h_rowDelimiters[numRows]-secStart;

    //std::cout<<"secStart: "<<secStart<<std::endl;

    floatType *d_val = csrDevice->getVal();
    int *d_cols = csrDevice->getCols();
    int *d_rowDelimiters = csrDevice->getRowDelimiters();
    floatType *d_vec = csrDevice->getVec();
    floatType *d_out = csrDevice->getOut();

#ifdef CUDA_RECORD
    // Setup events for timing
    cudaEvent_t start, stop;
    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));
   
    // Transfer data to device
    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
#endif
    CUDA_SAFE_CALL(cudaMemcpy(d_val, h_val,   numNonZeroes * sizeof(floatType),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_cols, h_cols, numNonZeroes * sizeof(int),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_rowDelimiters, h_rowDelimiters,(numRows+1) * sizeof(int), cudaMemcpyHostToDevice));
    
#ifdef CUDA_RECORD
    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));


    float iTransferTime, oTransferTime;
    CUDA_SAFE_CALL(cudaEventElapsedTime(&iTransferTime, start, stop));
    iTransferTime *= 1.e-3;
#endif
    // Bind texture for position
    string suffix;
    if (sizeof(floatType) == sizeof(float)){
        suffix = "-SP";
    }else {
        suffix = "-DP";
    }

    // Setup thread configuration
    int nBlocksScalar = (int) ceil((floatType) numRows / BLOCK_SIZE);
    int nBlocksVector = (int) ceil(numRows /(floatType)(BLOCK_SIZE / WARP_SIZE));
    int passes = op->getOptionInt("passes");
    int iters  = op->getOptionInt("iterations");
    

#ifdef CUDA_RECORD
    // Results description info
    char atts[TEMP_BUFFER_SIZE];
    sprintf(atts, "%d_elements_%d_rows",numNonZeroes, numRows);
    string prefix = "";
    double gflop = 2 * (double) numNonZeroes / 1e9;
#endif

#ifdef DARTS_DEBUG
    cout << "CSR Scalar Kernel\n";
#endif
    //cout<<"passes is : " <<passes<<", iters is "<< iters<<std::endl;

    //for (int k=0; k<passes; k++)
    //{
        // Run Scalar Kernel
    
#ifdef CUDA_RECORD
        CUDA_SAFE_CALL(cudaEventRecord(start, 0));
#endif
        //for (int j = 0; j < iters; j++)
        //{
            if(suffix == "-DP"){
                spmv_csr_scalar_section_kernel<floatType, texReaderDP><<<nBlocksScalar, BLOCK_SIZE>>>
            (d_val, d_cols, d_rowDelimiters, numRows, d_out,secStart);
            }else{
                spmv_csr_scalar_section_kernel<floatType, texReaderSP><<<nBlocksScalar, BLOCK_SIZE>>>
            (d_val, d_cols, d_rowDelimiters, numRows, d_out,secStart);
            }
        //}
       
#ifdef CUDA_RECORD
        CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
        CUDA_SAFE_CALL(cudaEventSynchronize(stop));
        float scalarKernelTime;
        CUDA_SAFE_CALL(cudaEventElapsedTime(&scalarKernelTime, start, stop));
        // Transfer data back to host
        CUDA_SAFE_CALL(cudaEventRecord(start, 0));
#endif
        CUDA_SAFE_CALL(cudaMemcpy(h_out, d_out, numRows * sizeof(floatType),cudaMemcpyDeviceToHost));
        
#ifdef CUDA_RECORD
        CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
        CUDA_SAFE_CALL(cudaEventSynchronize(stop));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&oTransferTime, start, stop));
#endif
        cudaThreadSynchronize();
        
#ifdef CUDA_RECORD
        oTransferTime *= 1.e-3;

        scalarKernelTime = (scalarKernelTime / (float)iters) * 1.e-3;
        double totalTransfer = iTransferTime + oTransferTime;
        string startPoint = std::to_string(csrDevice->getStartPoint());
        string testName = prefix+"CSR-Scalar"+suffix+"-startPoint-"+startPoint;
    
        resultDB->AddResult(testName, atts, "Gflop/s",gflop/(scalarKernelTime));
        resultDB->AddResult(testName, atts, "Gflop/s",gflop / (scalarKernelTime+totalTransfer));
        //resultDB->AddResult(testName+"_PCIe", atts, "Gflop/s",gflop / (scalarKernelTime+totalTransfer));
#endif
    //}
}


template <typename floatType>
void csrTestVector(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

    int deviceStart = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+deviceStart;
    int secStart = h_rowDelimiters[0] ;
    floatType *h_val = csrHost->getVal()+secStart;
    int *h_cols = csrHost->getCols()+secStart;
    floatType *h_vec = csrHost->getVec()+deviceStart;
    floatType *h_out = csrHost->getOut()+deviceStart;

    int numRows = csrDevice->getNumRows();
    int numNonZeroes = csrDevice->getNumNonZeroes();

    //std::cout<<"secStart: "<<secStart<<std::endl;

    floatType *d_val = csrDevice->getVal();
    int *d_cols = csrDevice->getCols();
    int *d_rowDelimiters = csrDevice->getRowDelimiters();
    floatType *d_vec = csrDevice->getVec();
    floatType *d_out = csrDevice->getOut();

    // Setup events for timing
    cudaEvent_t start, stop;
    CUDA_SAFE_CALL(cudaEventCreate(&start));
    CUDA_SAFE_CALL(cudaEventCreate(&stop));
    
    // Transfer data to device
    CUDA_SAFE_CALL(cudaEventRecord(start, 0));
    CUDA_SAFE_CALL(cudaMemcpy(d_val, h_val,   numNonZeroes * sizeof(floatType),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_cols, h_cols, numNonZeroes * sizeof(int),cudaMemcpyHostToDevice));
   // CUDA_SAFE_CALL(cudaMemcpy(d_vec, h_vec, numRows * sizeof(floatType),cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_rowDelimiters, h_rowDelimiters,(numRows+1) * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
    CUDA_SAFE_CALL(cudaEventSynchronize(stop));


    float iTransferTime, oTransferTime;
    CUDA_SAFE_CALL(cudaEventElapsedTime(&iTransferTime, start, stop));
    iTransferTime *= 1.e-3;
    // Bind texture for position
    string suffix;
    if (sizeof(floatType) == sizeof(float)){
     //   cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
     //   CUDA_SAFE_CALL(cudaBindTexture(0, vecTex, d_vec, channelDesc,numRows * sizeof(float)));
        suffix = "-SP";

    }
    else {
    //    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<int2>();
    //    CUDA_SAFE_CALL(cudaBindTexture(0, vecTexD, d_vec, channelDesc,numRows * sizeof(int2)));
        suffix = "-DP";
    }

    // Setup thread configuration
    int nBlocksScalar = (int) ceil((floatType) numRows / BLOCK_SIZE);
    int nBlocksVector = (int) ceil(numRows /(floatType)(BLOCK_SIZE / WARP_SIZE));
    int passes = op->getOptionInt("passes");
    int iters  = op->getOptionInt("iterations");
    
    // Results description info
    char atts[TEMP_BUFFER_SIZE];
    sprintf(atts, "%d_elements_%d_rows", numNonZeroes, numRows);
    string prefix = "";
    double gflop = 2 * (double) numNonZeroes / 1e9;
    cout << "CSR vector Kernel\n";

    //cout<<"passes is : " <<passes<<", iters is "<< iters<<std::endl;

    //for (int k=0; k<passes; k++)
    //{
        // Run Scalar Kernel
        CUDA_SAFE_CALL(cudaEventRecord(start, 0));
        //for (int j = 0; j < iters; j++)
        //{
            if(suffix == "-DP"){
                spmv_csr_vector_section_kernel<floatType, texReaderDP><<<nBlocksVector, BLOCK_SIZE>>>
            (d_val, d_cols, d_rowDelimiters, numRows, d_out,secStart);
            }else{
                spmv_csr_vector_section_kernel<floatType, texReaderSP><<<nBlocksScalar, BLOCK_SIZE>>>
            (d_val, d_cols, d_rowDelimiters, numRows, d_out,secStart);
            }
        //}
        CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
        CUDA_SAFE_CALL(cudaEventSynchronize(stop));
        float vectorKernelTime;
        CUDA_SAFE_CALL(cudaEventElapsedTime(&vectorKernelTime, start, stop));
        // Transfer data back to host
        CUDA_SAFE_CALL(cudaEventRecord(start, 0));
        CUDA_SAFE_CALL(cudaMemcpy(h_out, d_out, numRows * sizeof(floatType),cudaMemcpyDeviceToHost));
        CUDA_SAFE_CALL(cudaEventRecord(stop, 0));
        CUDA_SAFE_CALL(cudaEventSynchronize(stop));
        CUDA_SAFE_CALL(cudaEventElapsedTime(&oTransferTime, start, stop));
        cudaThreadSynchronize();
        
        vectorKernelTime = (vectorKernelTime / (float)iters) * 1.e-3;
        string testName = prefix+"CSR-Vector"+suffix;
        double totalTransfer = iTransferTime + oTransferTime;
        
        resultDB->AddResult(testName, atts, "Gflop/s",gflop/(vectorKernelTime));
        resultDB->AddResult(testName+"_PCIe", atts, "Gflop/s",gflop / (vectorKernelTime+totalTransfer));
    
    //}
}


// ****************************************************************************
// Function: spmv_csr_scalar_kernel
//
// Purpose:
//   Computes sparse matrix - vector multiplication on the GPU using
//   the CSR data storage format, using a thread per row of the sparse
//   matrix; based on Bell (SC09) and Baskaran (IBM Tech Report)
//
// Arguments:
//   val: array holding the non-zero values for the matrix
//   cols: array of column indices for each element of the sparse matrix
//   rowDelimiters: array of size dim+1 holding indices to rows of the matrix
//                  last element is the index one past the last
//                  element of the matrix
//   dim: number of rows in the matrix
//   out: output - result from the spmv calculation
//
// Returns:  nothing
//           out indirectly through a pointer
//
// Programmer: Lukasz Wesolowski
// Creation: June 28, 2010
//
// Modifications:
//
// ****************************************************************************
template <typename fpType, typename texReader>
__global__ void
spmv_csr_scalar_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out)
{
    int myRow = blockIdx.x * blockDim.x + threadIdx.x;
    texReader vecTexReader;
    if (myRow < dim)
    {
        fpType t = 0.0f;
        int start = rowDelimiters[myRow];
        int end = rowDelimiters[myRow+1];
        for (int j = start; j < end; j++)
        {
            int col = cols[j];
            t += val[j] * vecTexReader(col);
#ifdef DARTS_DEBUG
            if(threadIdx.x <20&&blockIdx.x ==0){
                printf("val[%d]=%lf, vecTexReader(%d)=%lf\n",j,val[j],col,vecTexReader(col));
            }
#endif
        }
        out[myRow] = t;
    }
}



template <typename fpType, typename texReader>
__global__ void
spmv_csr_scalar_section_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out,
                       const int secStart)
{
    int myRow = blockIdx.x * blockDim.x + threadIdx.x;
    texReader vecTexReader;
    if (myRow < dim)
    {
        fpType t = 0.0f;
        int start = rowDelimiters[myRow]-secStart;
        int end = rowDelimiters[myRow+1]-secStart;
        for (int j = start; j < end; j++)
        {
            int col = cols[j];
            t += val[j] * vecTexReader(col);

#ifdef DARTS_DEBUG
           // if(threadIdx.x <20&&blockIdx.x ==0){
           //     printf("val[%d]=%lf, vecTexReader(%d)=%lf\n",j,val[j],col,vecTexReader(col));
           // }
#endif
        }
        out[myRow] = t;
    }
}


// ****************************************************************************
// Function: spmv_csr_vector_kernel
//
// Purpose:
//   Computes sparse matrix - vector multiplication on the GPU using
//   the CSR data storage format, using a warp per row of the sparse
//   matrix; based on Bell (SC09) and Baskaran (IBM Tech Report)
//
// Arguments:
//   val: array holding the non-zero values for the matrix
//   cols: array of column indices for each element of the sparse matrix
//   rowDelimiters: array of size dim+1 holding indices to rows of the matrix
//                  last element is the index one past the last
//                  element of the matrix
//   dim: number of rows in the matrix
//   out: output - result from the spmv calculation
//
// Returns:  nothing
//           out indirectly through a pointer
//
// Programmer: Lukasz Wesolowski
// Creation: June 28, 2010
//
// Modifications:
//
// ****************************************************************************
template <typename fpType, typename texReader>
__global__ void
spmv_csr_vector_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out)
{
    // Thread ID in block
    int t = threadIdx.x;
    // Thread ID within warp
    int id = t & (warpSize-1);
    int warpsPerBlock = blockDim.x / warpSize;
    // One row per warp
    int myRow = (blockIdx.x * warpsPerBlock) + (t / warpSize);
    // Texture reader for the dense vector
    texReader vecTexReader;

    __shared__ volatile fpType partialSums[BLOCK_SIZE];

    if (myRow < dim)
    {
        int warpStart = rowDelimiters[myRow];
        int warpEnd = rowDelimiters[myRow+1];
        fpType mySum = 0;
        for (int j = warpStart + id; j < warpEnd; j += warpSize)
        {
            int col = cols[j];
            mySum += val[j] * vecTexReader(col);
        }
        partialSums[t] = mySum;

        // Reduce partial sums
        if (id < 16) partialSums[t] += partialSums[t+16];
        if (id <  8) partialSums[t] += partialSums[t+ 8];
        if (id <  4) partialSums[t] += partialSums[t+ 4];
        if (id <  2) partialSums[t] += partialSums[t+ 2];
        if (id <  1) partialSums[t] += partialSums[t+ 1];

        // Write result
        if (id == 0)
        {
            out[myRow] = partialSums[t];
        }
    }
}



template <typename fpType, typename texReader>
__global__ void
spmv_csr_vector_section_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out,
                       const int secStart)
{
    // Thread ID in block
    int t = threadIdx.x;
    // Thread ID within warp
    int id = t & (warpSize-1);
    int warpsPerBlock = blockDim.x / warpSize;
    // One row per warp
    int myRow = (blockIdx.x * warpsPerBlock) + (t / warpSize);
    // Texture reader for the dense vector
    texReader vecTexReader;

    __shared__ volatile fpType partialSums[BLOCK_SIZE];

    if (myRow < dim)
    {
        int warpStart = rowDelimiters[myRow]-secStart;
        int warpEnd = rowDelimiters[myRow+1]-secStart;
        fpType mySum = 0;
        for (int j = warpStart + id; j < warpEnd; j += warpSize)
        {
            int col = cols[j];
            mySum += val[j] * vecTexReader(col);
        }
        partialSums[t] = mySum;

        // Reduce partial sums
        if (id < 16) partialSums[t] += partialSums[t+16];
        if (id <  8) partialSums[t] += partialSums[t+ 8];
        if (id <  4) partialSums[t] += partialSums[t+ 4];
        if (id <  2) partialSums[t] += partialSums[t+ 2];
        if (id <  1) partialSums[t] += partialSums[t+ 1];

        // Write result
        if (id == 0)
        {
            out[myRow] = partialSums[t];
        }
    }
}



// ****************************************************************************
// Function: spmv_ellpackr_kernel
//
// Purpose:
//   Computes sparse matrix - vector multiplication on the GPU using
//   the ELLPACK-R data storage format; based on Vazquez et al (Univ. of
//   Almeria Tech Report 2009)
//
// Arguments:
//   val: array holding the non-zero values for the matrix in column
//   major format and padded with zeros up to the length of longest row
//   cols: array of column indices for each element of the sparse matrix
//   rowLengths: array storing the length of each row of the sparse matrix
//   dim: number of rows in the matrix
//   out: output - result from the spmv calculation
//
// Returns:  nothing directly
//           out indirectly through a pointer
//
// Programmer: Lukasz Wesolowski
// Creation: June 29, 2010
//
// Modifications:
//
// ****************************************************************************
template <typename fpType, typename texReader>
__global__ void
spmv_ellpackr_kernel(const fpType * __restrict__ val,
                     const int    * __restrict__ cols,
                     const int    * __restrict__ rowLengths,
                     const int dim, fpType * __restrict__ out)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    texReader vecTexReader;

    if (t < dim)
    {
        fpType result = 0.0f;
        int max = rowLengths[t];
        for (int i = 0; i < max; i++)
        {
            int ind = i*dim+t;
            result += val[ind] * vecTexReader(cols[ind]);
        }
        out[t] = result;
    }
}

template <typename fpType>
__global__ void
zero(fpType * __restrict__ a, const int size)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < size) a[t] = 0;
}



template void csrTestScalar<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void csrTestScalar<int>(ResultDatabase* resultDB, OptionParser* op, CSRMM<int> *csrHost, CSRMM<int> *csrDevice );

template void csrTestVector<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void csrTestVector<int>(ResultDatabase* resultDB, OptionParser* op, CSRMM<int> *csrHost, CSRMM<int> *csrDevice );
