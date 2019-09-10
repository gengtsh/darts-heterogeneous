#ifndef ARRAY_STATISTIC_H
#define ARRAY_STATISTIC_H

#include <iostream>
#include <algorithm>
#include <math.h>
#include <strings.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <string>

#include "conf.h"
#include "util.h"

using namespace std;

template<typename T>
class ArrayStats{

    public: 
        ArrayStats()
        {
            data   = NULL;
            info   = "";
            size   = 0 ;     
            sum    = 0 ;
            max    = 0 ; 
            maxIdx = 0 ; 
            min    = FLT_MAX ; 
            minIdx = 0 ; 
            mean   = 0 ; 
            md     = 0 ; 
            sd     = 0 ; 
            cv     = 0 ; 
            range  = 0 ;

        };
        
        
        ArrayStats(string nm)
        :info(nm)
        {
            data   = NULL;
            size   = 0 ;
            sum    = 0 ;
            max    = 0 ; 
            maxIdx = 0 ; 
            min    = FLT_MAX ; 
            minIdx = 0 ; 
            mean   = 0 ; 
            md     = 0 ; 
            sd     = 0 ; 
            cv     = 0 ; 
            range  = 0 ;

        };

        ArrayStats(string nm, T * da)
        :info(nm)
        ,data(da)
        {
            size   = 0 ; 
            sum    = 0 ;
            max    = 0 ; 
            maxIdx = 0 ; 
            min    = FLT_MAX ; 
            minIdx = 0 ; 
            mean   = 0 ; 
            md     = 0 ; 
            sd     = 0 ; 
            cv     = 0 ; 
            range  = 0 ;

        };
       
        T       *data;
        string  info;
        int     size;
        int64_t sum;
        T       max;
        int     maxIdx;
        T       min;
        int     minIdx;
        double  mean;
        double  md;//mean deviation
        double  sd;//standard deviation
        double  cv;//coefficient of variation
        T       range;
         
        void calcArrayStats(int start,int end);
        void print(void)const;

        ~ArrayStats()
        {

        };
};




#endif

