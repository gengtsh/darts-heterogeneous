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
#include "cusparse.h"

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


__global__ void
minusVal(int * __restrict__ a, const int val, const int n)
{
    //printf("blockIdx.x = %d, blockDim.x = %d, threadIdx.x = %d \n",blockIdx.x,blockDim.x,threadIdx.x);
    //printf("threadIdx.x = %d\n",threadIdx.x);
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) a[i] -= val;

}

template <typename floatType>
__global__ void
printVal(floatType * __restrict__ a, const int n)
{
    //printf("blockIdx.x = %d, blockDim.x = %d, threadIdx.x = %d \n",blockIdx.x,blockDim.x,threadIdx.x);
    //printf("threadIdx.x = %d\n",threadIdx.x);
     
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(blockIdx.x == 0 && threadIdx.x< n){
        printf("device val[%d] = %lf\n", i,(double)a[i]);
    }
}

// Forward declarations for kernels
template <typename fpType, typename texReader>
__global__ void
shoc_spmv_csr_scalar_kernel(const fpType * __restrict__ val,
                       const int    * __restrict__ cols,
                       const int    * __restrict__ rowDelimiters,
                       const int dim, fpType * __restrict__ out);


template <typename fpType, typename texReader>
__global__ void
shoc_spmv_csr_scalar_section_kernel(const fpType * __restrict__ val,
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
void SHOC_csrTestScalar(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

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
                shoc_spmv_csr_scalar_section_kernel<floatType, texReaderDP><<<nBlocksScalar, BLOCK_SIZE>>>
            (d_val, d_cols, d_rowDelimiters, numRows, d_out,secStart);
            }else{
                shoc_spmv_csr_scalar_section_kernel<floatType, texReaderSP><<<nBlocksScalar, BLOCK_SIZE>>>
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
void SHOC_csrStreamTestScalar(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

    int deviceStart = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+deviceStart;
    int secStart = h_rowDelimiters[0] ;
    floatType *h_val = csrHost->getVal()+secStart;
    int *h_cols = csrHost->getCols()+secStart;
    //floatType *h_vec = csrHost->getVec()+deviceStart;
    //floatType *h_vec = csrHost->getVec()+secStart;
    floatType *h_out = csrHost->getOut()+deviceStart;

    int numRows = csrDevice->getNumRows();
    //int numNonZeroes = csrDevice->getNumNonZeroes();
    int numNonZeroes = h_rowDelimiters[numRows]-secStart;
#ifdef DARTS_DEBUG 
    std::cout<<"deviceStart: "<<deviceStart<<",secStart: "<<secStart<<",numRows: "<<numRows<<std::endl;
#endif
    floatType *d_val = csrDevice->getVal();
    int *d_cols = csrDevice->getCols();
    int *d_rowDelimiters = csrDevice->getRowDelimiters();
    floatType *d_vec = csrDevice->getVec();
    floatType *d_out = csrDevice->getOut();


    // Bind texture for position
    string suffix;
    if (sizeof(floatType) == sizeof(float)){
        suffix = "-SP";
    }else {
        suffix = "-DP";
    }
#ifdef DARTS_DEBUG
    cout << "CSR Stream Scalar Kernel\n";
#endif

    int nStream = 32;
    cudaStream_t *stream;
    cudaEvent_t  *cuEvent;
    stream = new cudaStream_t[nStream];
    cuEvent = new cudaEvent_t[nStream];

    int chunk = numRows/nStream;
    int *sNumRows = new int[nStream];
    int *sNumNonZeroes = new int[nStream];
    int *svcStart = new int[nStream];
    int *srStart   = new int[nStream];
    int *ssStart   = new int[nStream];
    
    // Setup thread configuration
    int *nBlocksScalar = new int[nStream]; 

    for(int i=0; i<nStream; ++i){
        //CUDA_SAFE_CALL(cudaStreamCreateWithFlags(&stream[i],cudaStreamNonBlocking));
        CUDA_SAFE_CALL(cudaStreamCreate(&stream[i]));
        //CUDA_SAFE_CALL(cudaEventCreate(&cuEvent[i]));
        CUDA_SAFE_CALL(cudaEventCreateWithFlags(&cuEvent[i],cudaEventDisableTiming));
    
        sNumRows[i]     = (i==(nStream-1))?(numRows-i*chunk):(chunk);
        sNumNonZeroes[i]= (i==(nStream-1))?(h_rowDelimiters[numRows]-h_rowDelimiters[i*chunk]):(h_rowDelimiters[(i+1)*chunk]-h_rowDelimiters[i*chunk]);
        svcStart[i]    = h_rowDelimiters[i*chunk]-h_rowDelimiters[0] ;
        srStart[i]     = i*chunk;
        ssStart[i]     = h_rowDelimiters[i*chunk];
        nBlocksScalar[i] = (int) ceil((floatType) sNumRows[i] / BLOCK_SIZE);

    }

    for(int i=0; i<nStream; ++i){
#ifdef DARTS_DEBUG
        std::cout<<"stream: "<<i<<std::endl;
        std::cout<<"svcStart["<<i<<"] = "<<svcStart[i]<<",svcStart["<<i<<"] = "<<svcStart[i]<<",srStart["<<i<<"] = "<<srStart[i]<<",sNumRows["<<i<<"] = "<<sNumRows[i]<<",sNumNonZeroes["<<i<<"] = "<<sNumNonZeroes[i]<<std::endl;
        std::cout<<"d_val addr:  "<<d_val<<std::endl; 
#endif
        CUDA_SAFE_CALL(cudaMemcpyAsync(d_val + svcStart[i] , h_val+svcStart[i],   sNumNonZeroes[i] * sizeof(floatType),cudaMemcpyHostToDevice,stream[i]));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d_cols+ svcStart[i], h_cols+svcStart[i], sNumNonZeroes[i] * sizeof(int),cudaMemcpyHostToDevice,stream[i]));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d_rowDelimiters + srStart[i], h_rowDelimiters+srStart[i],(sNumRows[i]+1) * sizeof(int), cudaMemcpyHostToDevice,stream[i]));

        if(suffix == "-DP"){
            shoc_spmv_csr_scalar_section_kernel<floatType, texReaderDP><<<nBlocksScalar[i], BLOCK_SIZE, 0, stream[i]>>>
        (d_val+svcStart[i], d_cols+svcStart[i], d_rowDelimiters+srStart[i], sNumRows[i], d_out+srStart[i],ssStart[i]);
        }else{
            shoc_spmv_csr_scalar_section_kernel<floatType, texReaderSP><<<nBlocksScalar[i], BLOCK_SIZE,0, stream[i]>>>
        (d_val+svcStart[i], d_cols+svcStart[i], d_rowDelimiters+srStart[i], sNumRows[i], d_out+srStart[i],ssStart[i]);
        }
        CUDA_SAFE_CALL(cudaMemcpyAsync(h_out+srStart[i], d_out+srStart[i], sNumRows[i] * sizeof(floatType),cudaMemcpyDeviceToHost,stream[i]));


    }

    //CUDA_SAFE_CALL(cudaThreadSynchronize());
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    //std::cout<<"h_out addr: "<<h_out<<",h_out[1]="<<h_out[1]<<std::endl; 
    
    delete [] sNumRows ;
    delete [] sNumNonZeroes ;
    delete [] svcStart;
    delete [] srStart; 
    delete [] ssStart; 
    delete [] nBlocksScalar ; 

    for(int i=0; i<nStream; ++i){
        CUDA_SAFE_CALL(cudaStreamDestroy(stream[i]));
        CUDA_SAFE_CALL(cudaEventDestroy(cuEvent[i]));
    }
    delete [] stream;
    delete [] cuEvent;


#ifdef CUDA_RECORD

#endif
}
    

template <typename floatType>
void CuSparse_csrTest(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

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
    
    
    /* cusparse APIs */
    //cusparseStatus_t status;
    int devId;
    cudaDeviceProp prop;
    CUDA_SAFE_CALL(cudaGetDevice(&devId));
    CUDA_SAFE_CALL(cudaGetDeviceProperties( &prop, devId)) ;

    /* initialize cusparse library */
    cusparseHandle_t handle=0;
    CUSPARSE_SAFE_CALL(cusparseCreate(&handle));

#if CUDA_V10
    cusparseSpMatDescr_t matA;
    cusparseDnVecDescr_t vecX,vecY;
    void* dBuffer  = NULL;
    size_t bufferSize = 0;

    cudaDataType cuValueType;//for cuda v10
    if (suffix == "-SP"){
        cuValueType = CUDA_R_32F;
    }else if(suffix == "-DP"){
        cuValueType = CUDA_R_64F;
    }
    /*create sparse matrix A in CSR format */
    CUSPARSE_SAFE_CALL(cusparseCreateCsr(&matA,numRows,numRows,numNonZeroes,d_rowDelimiters,d_cols,d_val,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,cuValueType));
    /*create dense vector X */
    CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecX,numRows,d_vec,cuValueType));
    /*create dense vector Y */
    CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecY,numRows,d_out,cuValueType));

#elif CUDA_V9
    /* create and setup matrix descriptor */
    cusparseMatDescr_t descr = 0;
    CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descr));
    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

#endif

    const floatType alpha = 1;
    const floatType beta  = 0;

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
    cout << "CSR (cuSparse) Scalar Kernel\n";
#endif
    
#ifdef CUDA_RECORD
        CUDA_SAFE_CALL(cudaEventRecord(start, 0));
#endif

#if CUDA_V10
       // /*execuse SpMV */ 
       CUDA_SAFE_CALL(cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,&alpha,maxA,vecX,&beta,vecY,cuValueType, CUSPARSE_MV_ALG_DEFAULT,dBuffer)); 

#elif CUDA_V9       
        if(suffix == "-DP"){
            CUSPARSE_SAFE_CALL(cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,numRows,numRows,numNonZeroes,(const double *)&alpha,descr,(const double *)d_val,d_rowDelimiters,d_cols,(const double*)d_vec,(const double*)&beta,(double*)d_out));
        }else{
            CUSPARSE_SAFE_CALL(cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE,numRows,numRows,numNonZeroes,(const float *)&alpha,descr,(const float *)d_val,d_rowDelimiters,d_cols,(const float*)d_vec,(const float*)&beta,(float*)d_out));
        }
#endif

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

        /* destroy handle */
        CUSPARSE_SAFE_CALL(cusparseDestroy(handle));

#if CUDA_V10 
        /* destroy matrix/vec descriptor */
        CUSPARSE_SAFE_CALL(cusparseDestroySpMat(matA));
        CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
        CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY));

#elif CUDA_V9
        /* destroy matrix descriptor */
        CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descr));

#endif

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
void CuSparse_csrStreamTest(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice ){

    int deviceStart = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+deviceStart;
    int secStart = h_rowDelimiters[0] ;
    floatType *h_val = csrHost->getVal()+secStart;
    int *h_cols = csrHost->getCols()+secStart;
    floatType *h_vec = csrHost->getVec();
    floatType *h_out = csrHost->getOut()+deviceStart;

    int numRows = csrDevice->getNumRows();
    int numCols = csrDevice->getNumCols();
    
    //numRows = 5000;
    //int numNonZeroes = h_rowDelimiters[numRows]-secStart;
#ifdef DARTS_DEBUG 
    std::cout<<"deviceStart: "<<deviceStart<<",secStart: "<<secStart<<",numRows: "<<numRows<<std::endl;
#endif
    floatType *d_val = csrDevice->getVal();
    int *d_cols = csrDevice->getCols();
    int *d_rowDelimiters = csrDevice->getRowDelimiters();
    floatType *d_vec = csrDevice->getVec();
    floatType *d_out = csrDevice->getOut();


    // Bind texture for position
    string suffix;
    if (sizeof(floatType) == sizeof(float)){
        suffix = "-SP";
    }else {
        suffix = "-DP";
    }
#ifdef DARTS_DEBUG
    cout << "CSR Stream Scalar Kernel\n";
#endif

    int nStream = 32;// less than MAXSTREAM
    cudaStream_t *stream;
    cudaEvent_t  *cuEvent;
    stream = new cudaStream_t[nStream];
    cuEvent = new cudaEvent_t[nStream];

    int chunk = numRows/nStream;
    int *sNumRows = new int[nStream];
    int *sNumNonZeroes = new int[nStream];
    int *svcStart = new int[nStream];
    int *srStart   = new int[nStream];
    int *soStart   = new int[nStream];
    int *ssStart   = new int[nStream];
   
    floatType **d_val_sec               = new floatType*[nStream];
    floatType **d_out_sec               = new floatType*[nStream]; 
    int **d_cols_sec                    = new int*[nStream];
    int **d_rowDelimiters_sec           = new int*[nStream];

    floatType **h_val_sec               = new floatType*[nStream];
    floatType **h_out_sec               = new floatType*[nStream]; 
    int **h_cols_sec                    = new int*[nStream];
    int **h_rowDelimiters_sec           = new int*[nStream];

    int *nBlocks_sec                    = new int[nStream]; 

    /* cusparse APIs */
 
#ifdef DARTS_DEBUG
    int devId;
    cudaDeviceProp prop;
    CUDA_SAFE_CALL(cudaGetDevice(&devId));
    CUDA_SAFE_CALL(cudaGetDeviceProperties( &prop, devId)) ;
    std::cout<<"cuda device async Eng count: "<<prop.asyncEngineCount<<std::endl;
#endif

#if CUDA_V10
    /* create matrix and vec descriptor */
    cusparseSpMatDescr_t *matA = new cusparseSpMatDescr_t [nStream];
    cusparseDnVecDescr_t *vecY = new cusparseDnVecDescr_t [nStream];
#endif
    
    cusparseHandle_t *handle = new cusparseHandle_t [nStream] ;
    for(int i=0; i<nStream; ++i){
        //CUDA_SAFE_CALL(cudaStreamCreateWithFlags(&stream[i],cudaStreamNonBlocking));
        //CUDA_SAFE_CALL(cudaEventCreateWithFlags(&cuEvent[i],cudaEventDisableTiming));
        CUDA_SAFE_CALL(cudaStreamCreate(&stream[i]));
        CUDA_SAFE_CALL(cudaEventCreate(&cuEvent[i]));
    
        sNumRows[i]     = (i==(nStream-1))?(numRows-i*chunk):(chunk);
        if(sNumRows[i] == 0) continue;
        sNumNonZeroes[i]= (i==(nStream-1))?(h_rowDelimiters[numRows]-h_rowDelimiters[i*chunk]):(h_rowDelimiters[(i+1)*chunk]-h_rowDelimiters[i*chunk]);
        svcStart[i]    = h_rowDelimiters[i*chunk]-h_rowDelimiters[0] ;
        srStart[i]     = i*(chunk+1);
        soStart[i]     = i*chunk;
        ssStart[i]     = h_rowDelimiters[i*chunk];


        d_val_sec[i]              = d_val + svcStart[i] ;
        d_cols_sec[i]             = d_cols + svcStart[i];
        d_rowDelimiters_sec[i]   = d_rowDelimiters + srStart[i] ;
        d_out_sec[i]              = d_out + soStart[i]; 

        h_val_sec[i]              = h_val + svcStart[i];
        h_cols_sec[i]             = h_cols + svcStart[i];
        h_rowDelimiters_sec[i]    = h_rowDelimiters + soStart[i];
        h_out_sec[i]              = h_out + soStart[i]; 
        nBlocks_sec[i]            = ceil(sNumRows[i]/(float)BLOCK_SIZE) ; 

        CUSPARSE_SAFE_CALL(cusparseCreate(&handle[i]));
        CUSPARSE_SAFE_CALL(cusparseSetStream(handle[i],stream[i]));

#if CUDA_V10
        /*create sparse matrix A in CSR format */
        CUSPARSE_SAFE_CALL(cusparseCreateCsr(&matA[i],sNumRows[i],numCols,sNumNonZeroes[i],d_rowDelimiters_sec[i],d_cols_sec[i],d_val_sec[i],CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I,CUSPARSE_INDEX_BASE_ZERO,cuValueType));
        /*create dense vector Y */
        CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecY[i],sNumRows[i],d_out_sec[i],cuValueType));

#endif
    }

#if CUDA_V10
    void* dBuffer  = NULL;
    cudaDataType cuValueType;
    if (suffix == "-SP"){
        cuValueType = CUDA_R_32F;
    }else if(suffix == "-DP"){
        cuValueType = CUDA_R_64F;
    }

    /*create dense vector X */
    cusparseDnVecDescr_t vecX;
    CUSPARSE_SAFE_CALL(cusparseCreateDnVec(&vecX,numRows,d_vec,cuValueType));

#endif

    const floatType alpha = 1;
    const floatType beta  = 0;

    ///* initialize cusparse library */
    //cusparseHandle_t handle ;
    //CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
    
    cusparseMatDescr_t  descr;
    CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&descr));
    CUSPARSE_SAFE_CALL(cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL));
    CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO));
    int sNumRowDelimiter = 0;
    for(int i=0; i<nStream; ++i){
     
        //CUSPARSE_SAFE_CALL(cusparseSetStream(handle,stream[i])); 
#ifdef DARTS_DEBUG
        std::cout<<"stream: "<<i<<std::endl;
        std::cout<<"svcStart["<<i<<"] = "<<svcStart[i]<<",srStart["<<i<<"] = "<<srStart[i]<<",soStart["<<i<<"] = "<<soStart[i]<<",sNumRows["<<i<<"] = "<<sNumRows[i]<<",sNumNonZeroes["<<i<<"] = "<<sNumNonZeroes[i]<<std::endl;
#endif

        CUDA_SAFE_CALL(cudaMemcpyAsync(d_val + svcStart[i] , h_val+svcStart[i],   sNumNonZeroes[i] * sizeof(floatType),cudaMemcpyHostToDevice,stream[i]));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d_cols+ svcStart[i], h_cols+svcStart[i], sNumNonZeroes[i] * sizeof(int),cudaMemcpyHostToDevice,stream[i]));
        CUDA_SAFE_CALL(cudaMemcpyAsync(d_rowDelimiters + srStart[i], h_rowDelimiters+soStart[i],(sNumRows[i]+1) * sizeof(int), cudaMemcpyHostToDevice,stream[i]));
        
        minusVal<<<nBlocks_sec[i],BLOCK_SIZE,0,stream[i]>>>(d_rowDelimiters_sec[i],ssStart[i],sNumRows[i]+1); 


#if CUDA_V10
       // /*execuse SpMV */ 
       CUDA_SAFE_CALL(cusparseSpMV(handle[i],CUSPARSE_OPERATION_NON_TRANSPOSE,&alpha,maxA[i],vecX,&beta,vecY[i],cuValueType, CUSPARSE_MV_ALG_DEFAULT,dBuffer)); 

#elif CUDA_V9

        if(suffix == "-DP"){
            CUSPARSE_SAFE_CALL(cusparseDcsrmv(handle[i],CUSPARSE_OPERATION_NON_TRANSPOSE,sNumRows[i],numCols,sNumNonZeroes[i],(const double *)&alpha,descr,(const double *)(d_val_sec[i]),d_rowDelimiters_sec[i],d_cols_sec[i],(const double*)d_vec,(const double*)&beta,(double*)(d_out_sec[i])));
        }else{
            CUSPARSE_SAFE_CALL(cusparseScsrmv(handle[i],CUSPARSE_OPERATION_NON_TRANSPOSE,sNumRows[i],numCols,sNumNonZeroes[i],(const float *)&alpha,descr,(const float *)(d_val_sec[i]),d_rowDelimiters_sec[i],d_cols_sec[i],(const float*)d_vec,(const float*)&beta,(float*)(d_out_sec[i])));
        }
#endif

#ifdef DARTS_DEBUG 
        ////CUDA_SAFE_CALL(cudaDeviceSynchronize());
        //int num = 20;
        //printVal<floatType><<<nBlocks,BLOCK_SIZE,0, stream[i]>>>(d_val+svcStart[i], num);
        //printVal<int><<<nBlocks,BLOCK_SIZE,0,stream[i]>>>(d_cols+svcStart[i], num);
        //printVal<int><<<nBlocks,BLOCK_SIZE,0,stream[i]>>>(d_rowDelimiters+srStart[i], num);
        //
        //for(int j=0; j<num; ++j){
        //    printVal<floatType><<<nBlocks,BLOCK_SIZE,0,stream[i]>>>(d_vec+(h_cols+svcStart[i])[j] , 1);
        //}
        //printVal<floatType><<<nBlocks,BLOCK_SIZE,0,stream[i]>>>(d_out+soStart[i], num);

#endif
        //CUDA_SAFE_CALL(cudaEventRecord(cuEvent[i],stream[i]));
        //CUDA_SAFE_CALL(cudaStreamWaitEvent(stream[i],cuEvent[i],0));
        
        CUDA_SAFE_CALL(cudaMemcpyAsync(h_out_sec[i], d_out_sec[i], sNumRows[i] * sizeof(floatType),cudaMemcpyDeviceToHost,stream[i]));

#ifdef DARTS_DEBUG 
        //CUDA_SAFE_CALL(cudaDeviceSynchronize());
        //for(int j=0; j<10; ++j){
        //    std::cout<<"h_out["<<j<<"] = "<<(h_out+soStart[i])[j]<<std::endl;
        //}
#endif

    }

    //CUDA_SAFE_CALL(cudaThreadSynchronize());
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

#ifdef DARTS_DEBUG 
    //std::cout<<"h_out[0] = "<< h_out[0]<<",h_out[1] = "<< h_out[1]<<", h_out[2] = "<<h_out[2]<<std::endl;
#endif

#if CUDA_V10 
    /* destroy matrix/vec descriptor */
    CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecX));
    for(int i =0; i<nStream;++i){ 
        CUSPARSE_SAFE_CALL(cusparseDestroySpMat(matA[i]));
        CUSPARSE_SAFE_CALL(cusparseDestroyDnVec(vecY[i]));
    }
    delete [] matA;
    delete [] vecY;
#elif CUDA_V9
    /* destroy matrix descriptor */
    CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(descr));
    
#endif

    /* destroy handle */
    //CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    for(int i =0; i<nStream;++i){ 
        CUSPARSE_SAFE_CALL(cusparseDestroy(handle[i]));
    }
    delete [] handle;
    delete [] sNumRows ;
    delete [] sNumNonZeroes ;
    delete [] svcStart;
    delete [] srStart; 
    delete [] ssStart; 
    delete [] d_val_sec;           
    delete [] d_cols_sec;          
    delete [] d_rowDelimiters_sec;
    delete [] d_out_sec;          

    delete [] h_val_sec;           
    delete [] h_cols_sec;           
    delete [] h_rowDelimiters_sec;
    delete [] h_out_sec;           
    delete [] nBlocks_sec; 
    
    for(int i=0; i<nStream; ++i){
        CUDA_SAFE_CALL(cudaStreamDestroy(stream[i]));
        CUDA_SAFE_CALL(cudaEventDestroy(cuEvent[i]));
    }
    delete [] stream;
    delete [] cuEvent;


#ifdef CUDA_RECORD

#endif
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
// Function: shoc_spmv_csr_scalar_kernel
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
shoc_spmv_csr_scalar_kernel(const fpType * __restrict__ val,
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
shoc_spmv_csr_scalar_section_kernel(const fpType * __restrict__ val,
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
            //if(threadIdx.x <10&&blockIdx.x ==0){
            //    printf("val[%ld]=%g, vecTexReader(%ld)=%g\n",j,val[j],col,vecTexReader(col));
            //}
#endif
        }
        out[myRow] = t;

#ifdef DARTS_DEBUG
            //if(threadIdx.x <10&&blockIdx.x ==0){
            //    printf("out[%ld]=%g\n",out[myRow]);
            //}
#endif
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


template void SHOC_csrTestScalar<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void SHOC_csrTestScalar<float>(ResultDatabase* resultDB, OptionParser* op, CSRMM<float> *csrHost, CSRMM<float> *csrDevice );


template void SHOC_csrStreamTestScalar<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void SHOC_csrStreamTestScalar<float>(ResultDatabase* resultDB, OptionParser* op, CSRMM<float> *csrHost, CSRMM<float> *csrDevice );




template void CuSparse_csrTest<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void CuSparse_csrTest<float>(ResultDatabase* resultDB, OptionParser* op, CSRMM<float> *csrHost, CSRMM<float> *csrDevice );

template void CuSparse_csrStreamTest<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
//template void CuSparse_csrStreamTest<float>(ResultDatabase* resultDB, OptionParser* op, CSRMM<float> *csrHost, CSRMM<float> *csrDevice );


template void csrTestVector<double>(ResultDatabase* resultDB, OptionParser* op, CSRMM<double> *csrHost, CSRMM<double> *csrDevice );
template void csrTestVector<float>(ResultDatabase* resultDB, OptionParser* op, CSRMM<float> *csrHost, CSRMM<float> *csrDevice );
