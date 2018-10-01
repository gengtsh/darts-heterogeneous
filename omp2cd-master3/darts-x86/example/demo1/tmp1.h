#include "DARTS.h"
DEF_CODELET_ITER(TMP,0,SHORTWAIT);
DEF_CODELET(SN,0,SHORTWAIT);

DEF_TP(TMP_TP)
{

    int  tt1;
    int  tt2;

    TMP *aa;
    TMP bb1;
    TMP bb2;
    SN ss;
    int k;

    Codelet *signalUp;
    
    TMP_TP(int tt1, int tt2,Codelet *up)
        :tt1(tt1)
        ,tt2(tt2)
        ,ss(4,4,this,LONGWAIT)
        ,signalUp(up)
    
    {
        k=3;
        aa= new TMP[4];
        for(int i=0;i<4; ++i){
            aa[i] = TMP{0,1,this,SHORTWAIT,i};
            add(aa+i);
        }

    }

    ~TMP_TP(){

        delete []aa;
        std::cout<<"descructor!"<<std::endl;
    }
};
