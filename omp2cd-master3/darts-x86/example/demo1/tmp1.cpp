#include "tmp1.h"

#include <iostream>
#include <cstdint>

void TMP::fire(void){

    LOAD_FRAME(TMP_TP); 
    int Id = getID();
    RESET(aa[Id]);

    int tt1 = FRAME(tt1);
    int tt2 = FRAME(tt2);
    
    int tt3 = tt1+tt2;
    
    std::cout<<"aa: "<<Id<<std::endl;
    
    SYNC(ss);
    EXIT_TP();
}

void SN::fire(void){

    LOAD_FRAME(TMP_TP);
    RESET(ss);
    
    if(FRAME(k)-- >0){

        for(int i=0;i<4;++i){
            SYNC(aa[i]);
        }
    }else{
        SIGNAL(signalUp);
    }
    std::cout<<"SN!"<<std::endl;
    EXIT_TP();
}
