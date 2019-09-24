#include "ArrayStats.h"


template<typename T>
void ArrayStats<T>::calcArrayStats(int start, int end){
    if(data == NULL){
        std::cout<<"forget to initial data to arraystats"<<std::endl;
        exit(0); 
    
    }else{

        int64_t   tmp = 0;
        int64_t   tmp2= 0;
        size = end-start; 
        for (int i=start; i<end; ++i){
            T tt = data[i]; 
            sum +=tt;
            if(max < tt){
                max = tt;
                maxIdx = i;
            }
            if(min > tt){
                min = tt;
                minIdx = i;
            }
        }
        range = max - min;
        mean = sum/size;
        for (int i=start; i<end; ++i){
            T tt = data[i];
            tmp  +=pow(tt-mean,2);
            tmp2 +=abs(tt-mean);

        }
        sd = sqrt(tmp/size);
        cv = sd/mean;
        md = tmp2/size;
    }
}

template <typename T>
void ArrayStats<T>::print(void) const{
    
    std::cout   <<info  <<","
                <<size  <<","
                <<sum   <<","
                <<range <<"," 
                <<max   <<","
                <<maxIdx<<","
                <<min   <<","
                <<minIdx<<","
                <<mean  <<","
                <<md    <<","
                <<sd    <<","
                <<cv    <<","
                <<std::endl;
}



template class ArrayStats<double>;
template class ArrayStats<int>;

