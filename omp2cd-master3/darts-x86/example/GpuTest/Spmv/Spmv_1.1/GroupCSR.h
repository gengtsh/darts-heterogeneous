#ifndef GROUP_CSR_H_
#define GROUP_CSR_H_

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
#include "tbb/concurrent_vector.h"

#include "conf.h"

using namespace std;

template<typename T>
class GroupCSR
{
public:
    T *val;
    int *cols;
    int *rowDelimiters;

    int rstart;
    int rend;
    int vcnum;
    tbb::concurrent_vector<T> results;
   
    GroupCSR()
    {
        
#ifdef DARTS_DEBUG
      //  std::cout<<"GroupCSR empty!"<<std::endl;
#endif
        val = NULL;
        cols = NULL;
        rowDelimiters=NULL;

        rstart=0;
        rend=0;
        vcnum=0;
    }
    GroupCSR(int rstart,int rend, int vcnum)
    :rstart(rstart)
    ,rend(rend)
    ,vcnum(vcnum)
    {
#ifdef DARTS_DEBUG
        std::cout<<"GroupCSR normal!!!!!"<<std::endl;
#endif
        int rnum = rend - rstart;
        val = new T[vcnum];
        cols= new int[vcnum];
        rowDelimiters = new int[rnum+1];
    }
    
    GroupCSR& operator=(const GroupCSR &rhs){
        if(this != &rhs){
#ifdef DARTS_DEBUG
            std::cout<<"GroupCSR assign!!!!!"<<std::endl;
#endif
            rstart =rhs.rstart ; 
            rend   =rhs.rend   ;
            vcnum  =rhs.vcnum;
            results=rhs.results;
            
            int rnum = rend-rstart;
           
            val = new T[vcnum];
            cols= new int[vcnum];
            rowDelimiters = new int[rnum+1];
            std::copy(rhs.val,rhs.val+vcnum,val);
            //memcpy(val, rhs.val, vcnum*sizeof(T));
            memcpy(cols, rhs.cols, vcnum*sizeof(int));
            memcpy(rowDelimiters, rhs.rowDelimiters, (rnum+1)*sizeof(int));

        }

        return *this;
    }
    
    GroupCSR(const GroupCSR &rhs)
    :rstart (rhs.rstart )
    ,rend   (rhs.rend   )
    ,vcnum  (rhs.vcnum  )
    ,results(rhs.results)
    {

#ifdef DARTS_DEBUG
    //    std::cout<<"GroupCSR copy!!!!!"<<std::endl;
#endif
        int rnum = rend-rstart;
       
        //std::cout<<"allocate news"<<std::endl;
        val = new T[vcnum];
        cols= new int[vcnum];
        rowDelimiters = new int[rnum+1];
        std::copy(rhs.val,rhs.val+vcnum,val);
        //memcpy(val, rhs.val, vcnum*sizeof(T));
        memcpy(cols, rhs.cols, vcnum*sizeof(int));
        memcpy(rowDelimiters, rhs.rowDelimiters, (rnum+1)*sizeof(int));


#ifdef DARTS_DEBUG
      //  std::cout<<"GroupCSR copy finish!!!!!"<<std::endl;
#endif
    }

    virtual ~GroupCSR(){
#ifdef DARTS_DEBUG
     //   std::cout<<"~GroupCSR!!!!!"<<std::endl;
#endif
        if(val != NULL){
            delete []val;
        }
        if(cols != NULL){
            delete []cols;
        }
        if(rowDelimiters != NULL){
            delete []rowDelimiters;
        }

        val  = NULL;
        cols = NULL;
        rowDelimiters = NULL;
    }
    


};




#endif

