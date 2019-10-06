#include "CSRMM.h"


template<typename T>
T* CSRMM<T>::getVal(void){
    return val;
}

template<typename T>
int* CSRMM<T>::getCols(void){
    return cols;
}

template<typename T>
int* CSRMM<T>::getRowDelimiters(void){
    return rowDelimiters;
}

template<typename T>
T* CSRMM<T>::getVec(void){
    return vec;
}


template<typename T>
T* CSRMM<T>::getOut(void){
    return out;
}

template<typename T>
T* CSRMM<T>::getNewVal(void){
    return val+rowDelimiters[startPoint];
}

template<typename T>
int* CSRMM<T>::getNewCols(void){
    return cols+rowDelimiters[startPoint];
}

template<typename T>
int* CSRMM<T>::getNewRowDelimiters(void){
    return rowDelimiters+startPoint;
}

template<typename T>
T* CSRMM<T>::getNewVec(void){
    return vec;
}


template<typename T>
T* CSRMM<T>::getNewOut(void){
    return out+startPoint;
}



template<typename T>
int CSRMM<T>::getNumNonZeroes(void){
    return numNonZeroes;
}

template<typename T>
int CSRMM<T>::getNumRows(void){
    return numRows;
}


template<typename T>
int CSRMM<T>::getNumCols(void){
    return numCols;
}


template<typename T>
int CSRMM<T>::getStartPoint(void){
    return startPoint;
}

template<typename T>
string CSRMM<T>::getKind(void){
    return kind;
}

template<typename T>
string CSRMM<T>::getConfig(void){
    return config;
}
template<typename T>
void CSRMM<T>::assignNumNonZeroes(int nItems){
    numNonZeroes = nItems;
}

template<typename T>
void CSRMM<T>::assignNumRows(int nRows){
    numRows = nRows;
}

template<typename T>
void CSRMM<T>::assignNumCols(int nCols){
    numCols = nCols;
}

template<typename T>
void CSRMM<T>::assignStartPoint(int st){
    startPoint = st ;
}

template<typename T>
void CSRMM<T>::assignConfig(string name){
    config = name;
}

template<typename T>
void CSRMM<T>::assignConfigToCPU(void){
    config = "cpu";
}

template<typename T>
void CSRMM<T>::assignConfigToGPU(void){
    config = "gpu";
}

template<typename T>
void CSRMM<T>::assignConfigToHybrid(void){
    config = "hybrid";
}

template<typename T>
void CSRMM<T>::assignKind(string name){
    kind = name;
}

template<typename T>
void CSRMM<T>::assignKindToHost(void){
    kind = "host";
}

template<typename T>
void CSRMM<T>::assignKindToDevice(void){
    kind = "device";
}

template<typename T>
void CSRMM<T>::allocateHost(){
    if (config == "cpu"){
        val = new T [numNonZeroes*sizeof(T)] ;
        cols= new int [numNonZeroes*sizeof(int)];
        rowDelimiters = new int [(numRows+1)*sizeof(int)];
        vec = new T [numRows*sizeof(T)];
        out = new T [numRows*sizeof(T)];
    }else{
        CUDA_SAFE_CALL(cudaMallocHost(&val, numNonZeroes * sizeof(T)));
        CUDA_SAFE_CALL(cudaMallocHost(&cols, numNonZeroes * sizeof(int)));
        CUDA_SAFE_CALL(cudaMallocHost(&rowDelimiters, (numRows + 1) * sizeof(int)));
        CUDA_SAFE_CALL(cudaMallocHost(&vec, numRows * sizeof(T)));
        CUDA_SAFE_CALL(cudaMallocHost(&out, numRows * sizeof(T)));
    
    }
}

template<typename T>
void CSRMM<T>::allocateDevice(void){

      // Allocate device (GPU) memory
      CUDA_SAFE_CALL(cudaMalloc(&val,  numNonZeroes * sizeof(T)));
      CUDA_SAFE_CALL(cudaMalloc(&cols, numNonZeroes * sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc(&rowDelimiters, (NSTREAM+numRows+1) * sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc(&vec,  numRows * sizeof(T)));
      CUDA_SAFE_CALL(cudaMalloc(&out,  numRows * sizeof(T)));
}


template<typename T>
void CSRMM<T>::allocateDevice(int numNonZeroes, int numRows){

      // Allocate device (GPU) memory
      CUDA_SAFE_CALL(cudaMalloc(&val,  numNonZeroes * sizeof(T)));
      CUDA_SAFE_CALL(cudaMalloc(&cols, numNonZeroes * sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc(&rowDelimiters, (NSTREAM+numRows+1) * sizeof(int)));
      CUDA_SAFE_CALL(cudaMalloc(&vec,  numRows * sizeof(T)));
      CUDA_SAFE_CALL(cudaMalloc(&out,  numRows * sizeof(T)));
}

template<typename T>
void CSRMM<T>::allocateRowDelimiters(int num){
      CUDA_SAFE_CALL(cudaMalloc(&rowDelimiters, (num) * sizeof(int)));
}

template<typename T>
void CSRMM<T>::reallocateRowDelimiters(int num){
    CUDA_SAFE_CALL(cudaFreeHost(rowDelimiters));
    CUDA_SAFE_CALL(cudaMalloc(&rowDelimiters, (num) * sizeof(int)));
}

template<typename T>
void CSRMM<T>::allocateMem(void){
    if(kind == "host"){
        allocateHost();
    }else{
        allocateDevice();
    }
}

template<typename T>
void CSRMM<T>::freeHost(void){
    if(config == "cpu"){
        delete [] val;
        delete [] cols;
        delete [] rowDelimiters;
        delete [] vec;
        delete [] out;

    }else{

      // Free device memory
      CUDA_SAFE_CALL(cudaFreeHost(val));
      CUDA_SAFE_CALL(cudaFreeHost(cols));
      CUDA_SAFE_CALL(cudaFreeHost(rowDelimiters));
      CUDA_SAFE_CALL(cudaFreeHost(vec));
      CUDA_SAFE_CALL(cudaFreeHost(out));

    }
}

template<typename T>
void CSRMM<T>::freeDevice(void){

      // Free device memory
      CUDA_SAFE_CALL(cudaFree(val));
      CUDA_SAFE_CALL(cudaFree(cols));
      CUDA_SAFE_CALL(cudaFree(rowDelimiters));

      CUDA_SAFE_CALL(cudaFree(vec));
      CUDA_SAFE_CALL(cudaFree(out));

}


template<typename T>
void CSRMM<T>::freeMem(void){
    
       if(kind == "host"){
           freeHost();
       }else{
           freeDevice();
       }

}

template<typename T>
void CSRMM<T>::MemsetHostToZero(void){

    if(config == "cpu"){
        memset(val,0,numNonZeroes * sizeof(T));
        memset(cols,0,numNonZeroes * sizeof(int));
        memset(rowDelimiters,0,(numRows+1) * sizeof(int));
        memset(vec,0,numRows * sizeof(T));
        memset(out,0,numRows * sizeof(T));
    }

}


template<typename T>
void CSRMM<T>::MemsetDeviceToZero(void){
    CUDA_SAFE_CALL(cudaMemset(val,0,numNonZeroes * sizeof(T)));
    CUDA_SAFE_CALL(cudaMemset(cols,0,numNonZeroes * sizeof(int)));
    CUDA_SAFE_CALL(cudaMemset(rowDelimiters,0,(numRows+1) * sizeof(int)));
    CUDA_SAFE_CALL(cudaMemset(vec,0,numRows * sizeof(T)));
    CUDA_SAFE_CALL(cudaMemset(out,0,numRows * sizeof(T)));

}

template<typename T>
void CSRMM<T>::MemsetAllToZero(void){
       if(kind == "host"){
           MemsetHostToZero();
       }else{
           MemsetDeviceToZero();
       }
}


template<typename T>
template<typename T2>
void CSRMM<T>::Memset(T2* ptr, int constVal, int size){
    if(config == "cpu"){
        memset(ptr,constVal, size * sizeof(T2));
    }else{
        CUDA_SAFE_CALL(cudaMemset(ptr,constVal,size * sizeof(T2)));
    }
}

template<typename T>
void CSRMM<T>::print(void){
    std::cout<<"kind: "<<kind<<", config: "<<config<<std::endl;
    std::cout<<"numNonZeroes: "<<numNonZeroes<<std::endl;
    std::cout<<"NumRows: "<<numRows<<std::endl;
    std::cout<<"startPoint: "<<startPoint<<std::endl;
    std::cout<<"allocal val: "<<numNonZeroes<<"* sizeof "<<sizeof(T)<<std::endl;
    std::cout<<"allocate cols: "<<numNonZeroes<<"*sizeof "<<sizeof(int)<<std::endl;
    std::cout<<"allocate rowDelimiters: "<<numRows+1<<"*sizeof "<<sizeof(int)<<std::endl;
    std::cout<<"allocate vec: "<<numRows<<"*sizeof "<<sizeof(T)  <<std::endl;
    std::cout<<"allocate out: "<<numRows<<"*sizeof "<< sizeof(T) <<std::endl;
    std::cout<<"total allocate (MB): "<<(numNonZeroes*sizeof(T) + numNonZeroes*sizeof(int)+(numRows+1)*sizeof(int) + numRows*sizeof(T)*2)/(1024*1024)<<std::endl;
}

template class CSRMM<double>;
template void CSRMM<double>::Memset<double>(double *ptr, int constVal, int size);
template void CSRMM<double>::Memset<int>(int *ptr, int constVal, int size);

template class CSRMM<float>;
template void CSRMM<float>::Memset<float>(float *ptr, int constVal, int size);
template void CSRMM<float>::Memset<int>(int *ptr, int constVal, int size);
