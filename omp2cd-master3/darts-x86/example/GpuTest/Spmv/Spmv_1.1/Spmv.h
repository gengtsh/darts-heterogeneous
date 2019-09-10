#ifndef SPMV_H_
#define SPMV_H_


#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

#include "cudacommon.h"
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "tbb/concurrent_vector.h"

#include <cassert>
#include <pthread.h>

#include "OptionParser.h"
#include "ResultDatabase.h"
#include "RecordDatabase.h"
#include "util.h"
#include "conf.h"
#include "CSRMM.h"
#include "GroupCSR.h"

using namespace std;

extern pthread_mutex_t mutex2;
extern pthread_mutex_t mutex3;


enum kernelType{CSR_SCALAR, CSR_VECTOR, ELLPACKR};



template <typename floatType>
void spmvCPU(CSRMM<floatType> *csr);

template <typename floatType>
void spmv_cpu_kernel(CSRMM<floatType> *csr,int startRow, int endRow);

// ****************************************************************************
// Function: spmvCpu
//
// Purpose:
//   Runs sparse matrix vector multiplication on the CPU
//
// Arguements:
//   val: array holding the non-zero values for the matrix
//   cols: array of column indices for each element of A
//   rowDelimiters: array of size dim+1 holding indices to rows of A;
//                  last element is the index one past the last
//                  element of A
//   vec: dense vector of size dim to be used for multiplication
//   dim: number of rows/columns in the matrix
//   out: input - buffer of size dim
//        output - result from the spmv calculation
//
// Programmer: Lukasz Wesolowski
// Creation: June 23, 2010
// Returns:
//   nothing directly
//   out indirectly through a pointer
// ****************************************************************************

template <typename floatType>
void spmvCPU(CSRMM<floatType> *csr)
{
    const floatType *val = csr->getVal();
    const int *cols = csr->getCols();
    const int *rowDelimiters = csr->getRowDelimiters();
	const floatType *vec = csr->getVec(); 
    int dim = csr->getNumRows();
    floatType *out = csr->getOut();

    
    for (int i=0; i<dim; i++)
    {
       // if (i == 209752){
       //     std::cout<<"cpu: row "<<i<<", total "<< rowDelimiters[i+1]-rowDelimiters[i]<<"cols"<<std::endl;
       // 
       //  
       //     floatType t = 0;
       //     for (int j = rowDelimiters[i]; j < rowDelimiters[i + 1]; j++)
       //     {
       //         int col = cols[j];
       //         t += val[j] * vec[col];
       //         std::cout<<"val["<<j<<"] = "<<val[j]<<",vec["<<cols[j]<<"] ="<<vec[col]<< ", tmp = "<<val[j]*vec[col]<<std::endl;
       //     }
       //     std::cout<<"out["<<i<<"]="<<t<<std::endl;
       // }
        
        floatType t = 0;
        for (int j = rowDelimiters[i]; j < rowDelimiters[i + 1]; j++)
        {
            int col = cols[j];
            t += val[j] * vec[col];
             
        }
        out[i] = t;

    }
}



template <typename floatType>
void spmv_cpu_kernel(CSRMM<floatType> *csr,int startRow, int endRow)
{
    const floatType *val = csr->getVal();
    const int *cols = csr->getCols();
    const int *rowDelimiters = csr->getRowDelimiters();
	const floatType *vec = csr->getVec(); 
    int dim = csr->getNumRows();
    floatType *out = csr->getOut();
    
    for (int i=startRow; i<endRow; i++)
    {
        //if (i == (startRow+10)){
        //    std::cout<<"cpu: row "<<i<<", total "<< rowDelimiters[i+1]-rowDelimiters[i]<<" cols"<<std::endl;
        //
        // 
        //    floatType t = 0;
        //    for (int j = rowDelimiters[i]; j < rowDelimiters[i + 1]; j++)
        //    {
        //        int col = cols[j];
        //        t += val[j] * vec[col];
        //        std::cout<<"val["<<j<<"] = "<<val[j]<<",vec["<<cols[j]<<"] ="<<vec[col]<< ", tmp = "<<val[j]*vec[col]<<std::endl;
        //    }
        //    std::cout<<"out["<<i<<"]="<<t<<std::endl;
        //}
        
        floatType t = 0;
        for (int j = rowDelimiters[i]; j < rowDelimiters[i + 1]; j++)
        {
            int col = cols[j];
            t += val[j] * vec[col];
             
        }
        out[i] = t;

    }
}



template <typename floatType>
void spmv_cpu_group_kernel(floatType *vec,GroupCSR<floatType> &gcsr)
{
    const floatType *val = gcsr.val;
    const int *cols = gcsr.cols;
    const int *rowDelimiters = gcsr.rowDelimiters;
    int dim = gcsr.rend-gcsr.rstart;
    tbb::concurrent_vector<floatType> *out = &gcsr.results;
#ifdef DARTS_DEBUG
    
    pthread_mutex_lock(&mutex2);
    std::cout<<"gcsr.rstart: "<<gcsr.rstart<<" ,gcsr.rend: "<<gcsr.rend<<", gcsr.vcnum: "<<gcsr.vcnum<<std::endl;
    std::cout<<"gcsr: rowDelimiters[0]:"<<rowDelimiters[0]<<" ,gcsr.rend: "<<gcsr.rend<<", gcsr.vcnum: "<<gcsr.vcnum<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif
    
    //if(gcsr.rstart == 4641 || gcsr.rstart == 4642){
    //    int rst = gcsr.rstart;
    //    int red = gcsr.rend;
    //    //string csvv = "testv"+to_string(rst)+".csv";
    //    //string csvc = "testc"+to_string(rst)+".csv";
    //    string csvall = "testall"+to_string(rst)+".csv";
    //    ofstream out;
    //    out.open(csvall, std::ofstream::out);

    //    out<<"val:"<<"\n";
    //    //for(int i=rowDelimiters[0]; i<rowDelimiters[red-rst]; ++i){
    //    for(int i=rowDelimiters[0]; i<rowDelimiters[0]+20; ++i){
    //        out<<val[i]<<",";
    //    }
    //    out<<"cols: "<<"\n";
    //    //for(int i=rowDelimiters[0]; i<rowDelimiters[red-rst]; ++i){
    //    for(int i=rowDelimiters[0]; i<rowDelimiters[0]+20; ++i){
    //        out<<cols[i]<<",";
    //    }
    //    out<<endl;
    //    out.close();
    //}


    for (int i=0; i<dim; i++)
    {
        
        floatType t = 0;
        for (int j = rowDelimiters[i]; j < rowDelimiters[i + 1]; j++)
        {
        
#ifdef DARTS_DEBUG
            pthread_mutex_lock(&mutex2);
            
            if(gcsr.rstart == 7212){
                std::cout<<"cols["<<j<<"] = "<<cols[j]<<std::endl; 
                std::cout<<"val["<<j<<"] = "<<val[j]<<std::endl; 
                std::cout<<"vec["<<cols[j]<<"] = "<<vec[cols[j]]<<std::endl; 
            }
            pthread_mutex_unlock(&mutex2);
#endif

            int col = cols[j];
            t += val[j] * vec[col];
            
        }
        out->push_back(t);

    }
}



// ****************************************************************************
// Function: verifyResults
//
// Purpose:
//   Verifies correctness of GPU results by comparing to CPU results
//
// Arguments:
//   refResults: array holding the CPU result vector

// ****************************************************************************
// Function: verifyResults
//
// Purpose:
//   Verifies correctness of GPU results by comparing to CPU results
//
// Arguments:
//   refResults: array holding the CPU result vector
//   testResults: array hodling the GPU result vector
//   size: number of elements per vector
//   pass: optional iteration number
//
// Programmer: Lukasz Wesolowski
// Creation: June 23, 2010
// Returns:
//   nothing
//   prints "Passed" if the vectors agree within a relative error of
//   MAX_RELATIVE_ERROR and "FAILED" if they are different
// ****************************************************************************
template <typename floatType>
bool verifyResults(const floatType *refResults, const floatType *testResults,
                   const int size, const int pass = -1)
{
    bool passed = true;
    for (int i = 0; i < size; i++)
    {
        if (fabs(refResults[i] - testResults[i]) / refResults[i]
            > MAX_RELATIVE_ERROR)
        {
            cout << "Mismatch at i: "<< i << " ref: " << refResults[i] <<
                " dev: " << testResults[i] << endl;
            passed = false;
            return passed;
        }
    }
    if (pass != -1)
    {
        cout << "Pass "<<pass<<": ";
    }
    if (passed)
    {
        cout << "Passed" << endl;
    }
    else
    {
        cout << "---FAILED---" << endl;
    }
    return passed;
}

template <typename floatType>
bool calcStaticNext(int &totalNumRows, int &numRowsLeft, CSRMM<floatType> *csrHost,CSRMM<floatType> *csrDevice, string Id ){
    
    int *rowDelimiters = csrHost->getRowDelimiters();
    int newRowStart = totalNumRows - numRowsLeft;

    int cpuNumRows = csrHost->getNumRows();
    int gpuNumRows = csrDevice->getNumRows();
    
    if(Id == "cpu"){
        cpuNumRows = (numRowsLeft>cpuNumRows)?cpuNumRows:numRowsLeft;
        csrHost->assignStartPoint(newRowStart);
        csrHost->assignNumRows(cpuNumRows);

        numRowsLeft -=cpuNumRows; 

    }else if(Id == "gpu"){

        gpuNumRows = (numRowsLeft>gpuNumRows)?gpuNumRows:numRowsLeft;
        csrDevice->assignStartPoint(newRowStart);
        csrDevice->assignNumRows(gpuNumRows);

        numRowsLeft -=gpuNumRows; 
    }

    return true;

}


template <typename floatType>
void adjustCSR(GroupCSR<floatType> &gcsr,CSRMM<floatType> *csr, int rIdxStart, int rIdxEnd){

    floatType *val = csr->getVal();
    int *cols = csr->getCols();
    int *rowDelimiters = csr->getRowDelimiters();

    int numRows = csr->getNumRows();
   
    int vcNum      = rowDelimiters[numRows]; 
    int vcIdxStart = rowDelimiters[rIdxStart];
    int vcIdxEnd   = rowDelimiters[rIdxEnd];
    
    int vcdiff      = vcIdxEnd-vcIdxStart;
    
    floatType *gval = gcsr.val;
    int *gcols      = gcsr.cols;
    int *growDelimiters = gcsr.rowDelimiters;


#ifdef DARTS_DEBUG
    std::cout<<"numRows = "<<numRows<<", vcNum = "<<vcNum<<std::endl;
    std::cout<<"rIdxStart: "<<rIdxStart<<", rIdXEnd "<<rIdxEnd<<", rowDelimiters["<<rIdxStart<<"] = " <<rowDelimiters[rIdxStart]<<", rowDelimiters["<<rIdxEnd<<"] = "<<rowDelimiters[rIdxEnd]<<std::endl;
#endif

#ifdef DARTS_DEBUG
    if(rIdxStart == 7212 ){
        for(int i=rIdxStart;i<=rIdxStart+20;++i){
            std::cout<<"before: rowDelimiters["<<i<<"] = "<<rowDelimiters[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxStart];i<=rowDelimiters[rIdxStart]+20;++i){
            std::cout<<"before: cols["<<i<<"] = "<<cols[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxEnd];i<=rowDelimiters[rIdxEnd]+20;++i){
            std::cout<<"before: cols["<<i<<"] = "<<cols[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxStart];i<=rowDelimiters[rIdxStart]+20;++i){
            std::cout<<"before: val["<<i<<"] = "<<val[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxEnd];i<=rowDelimiters[rIdxEnd]+20;++i){
            std::cout<<"before: val["<<i<<"] = "<<val[i]<<std::endl;
        }
    }

#endif

#ifdef DARTS_DEBUG
    std::cout<<"memcpy val/cols from "<<vcIdxStart<<" to "<<vcIdxEnd<<" , to groupcsr val/cols, vcnum: "<<gcsr.vcnum<<std::endl;
#endif
    std::copy( val+vcIdxStart, val+vcIdxEnd, gval);
    std::copy(cols+vcIdxStart,cols+vcIdxEnd, gcols);
    //memcpy(gval, val+vcIdxStart, vcdiff*sizeof(floatType));
    //memcpy(gcols,cols+vcIdxStart,vcdiff*sizeof(int));

#ifdef DARTS_DEBUG
    std::cout<<"memcpy val/cols from "<<vcIdxEnd<<", sz: "<<vcNum-vcIdxEnd<<" to "<<vcIdxStart<<std::endl;
#endif

#ifdef DARTS_DEBUG
    if(rIdxStart == 7212 ){
        
        for(int i = vcIdxEnd;i<vcIdxEnd+20; ++i){
            std::cout<<"before: memcpy col["<<i<<"]="<<cols[i]<<std::endl;
        }
        for(int i = vcIdxEnd;i<vcIdxEnd+20; ++i){
            std::cout<<"before: memcpy val["<<i<<"]="<<val[i]<<std::endl;
        }

    }


#endif

#ifdef DARTS_DEBUG
    std::cout<<"memcpy val : src addr:   "<<val +vcIdxEnd<<", dst addr: "<<val +vcIdxStart<<" ,num: "<<vcNum-vcIdxEnd<<std::endl;
    std::cout<<"memcpy cols: src addr:   "<<cols+vcIdxEnd<<", dst addr: "<<cols+vcIdxStart<<" ,num: "<<vcNum-vcIdxEnd<<std::endl;
#endif

    std::copy(val +vcIdxEnd, val+vcNum, val+vcIdxStart);
    std::copy(cols+vcIdxEnd,cols+vcNum,cols+vcIdxStart);
    //memcpy(val +vcIdxStart,val +vcIdxEnd, (vcNum-vcIdxEnd)*sizeof(floatType));
    //memcpy(cols+vcIdxStart,cols+vcIdxEnd, (vcNum-vcIdxEnd)*sizeof(int));
 
#ifdef DARTS_DEBUG
    if(rIdxStart == 7212){
        
        for(int i = vcIdxStart;i<vcIdxStart+20; ++i){
            std::cout<<"after: memcpy col["<<i<<"]="<<cols[i]<<std::endl;
        }
        for(int i = vcIdxStart;i<vcIdxStart+20; ++i){
            std::cout<<"after memcpy val["<<i<<"]="<<val[i]<<std::endl;
        }

    }

#endif

    int rval = rowDelimiters[rIdxStart];
   
#ifdef DARTS_DEBUG
    std::cout<<"copy rowDelimiter "<<rIdxStart<<" to "<<rIdxEnd<<" to groupcsr.rowDelimiter "<<rval<<std::endl;
#endif
        
    int rsz = rIdxEnd-rIdxStart;

    std::transform(rowDelimiters+rIdxStart,rowDelimiters+rIdxEnd+1,growDelimiters, [vcIdxStart](int x){return x-vcIdxStart;});

#ifdef DARTS_DEBUG
    std::cout<<"memset rowDelimiters idx from "<<rIdxStart+1<<" to "<<rIdxEnd<<" value: "<<rval<<std::endl;
#endif
    //for(int i=rIdxStart+1; i<=rIdxEnd;++i){
    //   rowDelimiters[i] = rval;
    //}
    std::fill(rowDelimiters+rIdxStart,rowDelimiters+rIdxEnd+1,rval);

#ifdef DARTS_DEBUG
    std::cout<<"rowDelimiters after "<<rIdxEnd+1<<" will minux  "<<vcdiff<<" until to "<<numRows<<std::endl;
#endif

    //for(int i=rdIdxEnd+1; i<numRows+1;++i){
    //    rowDelimiters[i]-=vcdiff;
    //}
    std::transform(rowDelimiters+rIdxEnd+1,rowDelimiters+numRows+1, rowDelimiters+rIdxEnd+1,[vcdiff](int x){return x-vcdiff;});


#ifdef DARTS_DEBUG
    //for(int i=rIdxStart;i<=rIdxEnd;++i){
    if(rIdxStart == 7212){
        for(int i=rIdxStart;i<=rIdxStart+20;++i){
            std::cout<<"after: rowDelimiters["<<i<<"] = "<<rowDelimiters[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxStart];i<=rowDelimiters[rIdxStart]+20;++i){
            std::cout<<"after: cols["<<i<<"] = "<<cols[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxEnd];i<=rowDelimiters[rIdxEnd]+20;++i){
            std::cout<<"after: cols["<<i<<"] = "<<cols[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxStart];i<=rowDelimiters[rIdxStart]+20;++i){
            std::cout<<"after: val["<<i<<"] = "<<val[i]<<std::endl;
        }

        for(int i=rowDelimiters[rIdxEnd];i<=rowDelimiters[rIdxEnd]+20;++i){
            std::cout<<"after: val["<<i<<"] = "<<val[i]<<std::endl;
        }
    }
#endif

#ifdef DARTS_DEBUG
    std::cout<<"rowDelimiters["<<numRows<<"]= "<<rowDelimiters[numRows]<<std::endl;
#endif

}

template<typename T>
void addGroupResultsBack(T *dst, vector<GroupCSR<double>> *gsCSR ){
#ifdef DARTS_DEBUG 
    std::cout<<"addGroupResultsBack"<<std::endl;
#endif

    int gsz = gsCSR->size();
    for (int i=0; i<gsz; ++i){
       
        int rstart = gsCSR->at(i).rstart;
        int rend   = gsCSR->at(i).rend;
        tbb::concurrent_vector<T>  r = gsCSR->at(i).results;
//        std::cout<<"rstart: "<<rstart<<std::endl;
        for(int j=0;j<r.size();++j){
            dst[rstart+j] = r[j];
//            std::cout<<r[j]<<",";
        }
//        std::cout<<std::endl;
    }
}


void preProcessNonZeroes(vector<pair<int,int>> &dst,int *src,int st, int ed);

void preProcessNonZeroes(int *dst,int *src,int st, int ed);

void preProcessNonZeroes(vector<pair<int,int>> &dstVec,int* dstArray,int *src,int st, int ed);

struct GroupAttr{
    int start;
    int end;
    int groupId;
};


void calcGroups(int *refArray,vector<vector<pair<int,int>>> &groups,vector<pair<int,int>> &srcVec, int num, double pt,int thd =0);


void memsetVal(int *src, int val, int st, int ed);
int buildGroups(vector<GroupAttr> &groups,int *src, int st, int ed, int topn, double pt, int thd=0);
int buildGroups_v2(vector<GroupAttr> &groups,int *src, int st, int ed, int topn, double pt, int thd=0);

void addBenchmarkSpecOptions(OptionParser &op);

template <typename floatType>
void csrTestScalar(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice );


template <typename floatType>
void csrTestVector(ResultDatabase* resultDB, OptionParser* op, CSRMM<floatType> *csrHost, CSRMM<floatType> *csrDevice );

template <typename floatType>
void memcpyHostToDevice(floatType *dst, floatType *src, int size );

template <typename floatType>
void memcpyDeviceTexture(const void* devPtr, size_t size );

#endif // SPMV_H_
