#ifndef CSRMM_H_
#define CSRMM_H_

#include <iostream>
#include <algorithm>
#include <math.h>
#include <strings.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <cstring>

#include "cudacommon.h"
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "conf.h"

using namespace std;

template<typename T>
class CSRMM
{
private:
    T *val;
    int *cols;
    int *rowDelimiters;
    T *vec;
    T *out;
    std::string kind;//host or device
    std::string config;
    int numNonZeroes ;
    int numRows ;
    int numCols ;
    int startPoint;
public:

    CSRMM()
    {
        //kind   = "host";
        //config = "cpu";
        assignKindToHost();
        assignConfigToCPU();
        numNonZeroes = 0;
        numRows = 0;
        numCols = 0;
        startPoint = 0;
    }
    
    CSRMM(const std::string& kind)
    :kind(kind)
    {
        //config = "cpu";
        assignConfigToCPU();
        numNonZeroes = 0;
        numRows = 0;
        numCols = 0;
        startPoint= 0;
    }
    
    CSRMM(const std::string &kind,const std::string &config)
    :kind(kind)
    ,config(config)
    {
        numNonZeroes = 0;
        numRows = 0;
        numCols = 0;
        startPoint = 0;
    }
    
    
    CSRMM(string kind,string config, int numNonZeroes, int numRows)
    :kind(kind)
    ,config(config)
    ,numNonZeroes(numNonZeroes)
    ,numRows(numRows)
    {
        startPoint = 0;
        numCols = numRows;
    }

    CSRMM(string kind,string config, int numNonZeroes, int numRows,int startPoint)
    :kind(kind)
    ,config(config)
    ,numNonZeroes(numNonZeroes)
    ,numRows(numRows)
    ,startPoint(startPoint)
    {
        numCols = numRows;
    }


    CSRMM(string kind,string config, int numNonZeroes, int numRows,int numCols,int startPoint)
    :kind(kind)
    ,config(config)
    ,numNonZeroes(numNonZeroes)
    ,numRows(numRows)
    ,numCols(numCols)
    ,startPoint(startPoint)
    {
    }


    virtual ~CSRMM(){
#ifdef DARTS_DEBUG
        //std::cout<<" ~CSRMM! kind: "<<kind<<std::endl;
#endif
       // if(kind == "host"){
       //     freeHost();
       // }else{
       //     freeDevice();
       // }
    }
    
    T* getVal(void);
    int* getCols(void);
    int* getRowDelimiters(void);
    T* getVec(void);
    T* getOut(void);
    
    T* getNewVal(void);
    int* getNewCols(void);
    int* getNewRowDelimiters(void);
    T* getNewVec(void);
    T* getNewOut(void);
    
    int getNumNonZeroes(void);
    int getNumRows(void);
    int getNumCols(void);
    int getStartPoint(void);
    string getKind(void);
    string getConfig(void);
    void assignNumNonZeroes(int nItems);
    void assignNumRows(int nRows);
    void assignNumCols(int nCols);
    void assignStartPoint(int st);
    void assignConfig(string name);
    void assignConfigToCPU(void);
    void assignConfigToGPU(void);
    void assignConfigToHybrid(void);

    void assignKind(string name);
    void assignKindToHost(void);
    void assignKindToDevice(void);


    void allocateMem(void);
    void allocateHost(void);
    void allocateDevice(void);
    void allocateDevice(int numNonZeroes,int numRows);
    void allocateRowDelimiters(int num);
    void reallocateRowDelimiters(int num);
    
    void freeHost(void);
    void freeDevice(void);
    void print(void);        

    void freeMem(void);

    void MemsetAllToZero(void);
    void MemsetHostToZero(void);
    void MemsetDeviceToZero(void);

    template <typename T2>
    void Memset(T2* ptr, int constVal, int size);
    

    //template <typename T2>
    //void Memset(T2* ptr, int constVal, int size){

    //    if(config == "cpu"){
    //        memset(ptr,constVal, size * sizeof(T2));
    //    }else{
    //        CUDA_SAFE_CALL(cudaMemset(ptr,constVal,size * sizeof(T2)));
    //    }
    //}


};




#endif

