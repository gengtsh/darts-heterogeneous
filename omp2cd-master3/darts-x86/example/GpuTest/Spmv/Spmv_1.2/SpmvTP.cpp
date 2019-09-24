#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include <cassert>
#include <pthread.h>
pthread_mutex_t mutex1;
pthread_mutex_t mutex2;

pthread_mutex_t mutex3;
//#include <sstream>
#include <iostream>

#include "SpmvTP.h"

//	pthread_mutex_lock(&mutex);
//	pthread_mutex_unlock(&mutex);

uint64_t ttStart, ttEnd;


void 
SpmvCpuCD::fire(void) 
{

//#ifdef DARTS_RECORD 
    uint64_t start, end;
//#endif

//#ifdef DARTS_RECORD 
    start = getTime();
//#endif
    
    LOAD_FRAME(SpmvTP);
    RESET(spmvCpu);
#ifdef DARTS_DEBUG
	std::cout<<"Invoke SpmvCpu"<<std::endl;	
#endif

    
    CSRMM<double>*csrHost =&FRAME(csrHost);
    //csrHost->print();


    spmvCPU(csrHost);


    
    double *refOut = FRAME(refOut);
    double *out = csrHost->getOut();
    int numRows = csrHost->getNumRows();
    memcpy(refOut,out, (numRows)*sizeof(double));

    csrHost->Memset<double>(out,0,numRows);

//#ifdef DARTS_RECORD 
    end = getTime();
//#endif    

//#ifdef DARTS_RECORD

    //int numNonZeroes = csrHost->getNumNonZeroes();
    //Record rd("ref", numRows,numNonZeroes,start,end);
    //
    //RecordDatabase *recordDB = FRAME(recordDB);
    //recordDB->assignRef(rd);
    //recordDB->printRef();
    std::cout<<"ref: "<<end-start<<std::endl;
//#endif


#ifdef DARTS_DEBUG_VAL  
    double *val = csrHost->getVal();
    int *cols = csrHost->getCols();
    int *rowDelimiters = csrHost->getRowDelimiters();
    double *vec = csrHost->getVec();
    for (size_t i=0;i<20 ; ++i){
        std::cout<<"refVal["<<i<<"] = "<<val[i]<<std::endl;
    }
  
    for (size_t i=0;i<20 ; ++i){
        std::cout<<"refCols["<<i<<"] = "<<cols[i]<<std::endl;
    }

    for (size_t i=0;i<20 ; ++i){
        std::cout<<"rowDelimiters["<<i<<"] = "<<rowDelimiters[i]<<std::endl;
    }
   
    for (size_t i=0;i<20 ; ++i){
        std::cout<<"refVec["<<i<<"] = "<<vec[i]<<std::endl;
    }

    for (size_t i=0;i<20 ; ++i){
        std::cout<<"refOut["<<i<<"] = "<<refOut[i]<<std::endl;
    }
   
    int n = rowDelimiters[2000];
    std::cout<<"rowDelimiters[2000] = "<<n<<std::endl;
    int i,j;
    for (i= n;i<n+20 ; ++i){
        std::cout<<"refVal["<<i<<"] = "<<val[i]<<std::endl;
    }

    for (i= n;i<n+20 ; ++i){
        std::cout<<"refCols["<<i<<"] = "<<cols[i]<<std::endl;
    }

    for (size_t i=n;i<n+20 ; ++i){
        std::cout<<"rowDelimiters["<<i<<"] = "<<rowDelimiters[i]<<std::endl;
    }

    for(size_t i=2000; i< 2020; ++i){
        std::cout<<"refOut["<<i<<"] = "<<refOut[i]<<std::endl;
    }
#endif

    SYNC(spmvCpuSync);
    EXIT_TP();
}


void 
SpmvCpuSyncCD::fire(void) 
{
	LOAD_FRAME(SpmvTP);
    RESET(spmvCpuSync);
#ifdef DARTS_DEBUG
	std::cout<<"Invoke SpmvCpuSync"<<std::endl;	
#endif
    CSRMM<double>*csrHost =&FRAME(csrHost);
    CSRMM<double>*csrDevice =&FRAME(csrDevice);

    string config = FRAME(config);
    
    int *rowDelimiters = csrHost->getRowDelimiters();

    int totalNumRows = FRAME(totalNumRows);
    int nCnt = FRAME(nCnt); 
    FRAME(numRowsLeft) = totalNumRows ;
    
    //========== configue csrHost and csrDevice ===========//
    if(config == "gpu" ){
        CSRMM<double>*csrDevice = &FRAME(csrDevice);
        //int gpuNumRows = FRAME(gpuNumRows) ;
        int gpuNumRows = FRAME(totalNumRows);
        
        //int gpuNumRows = 10000;
        //FRAME(gpuNumRows) = gpuNumRows;
        csrDevice->assignNumRows(gpuNumRows);
        
        __sync_fetch_and_sub(&FRAME(numRowsLeft),gpuNumRows); 
        if(FRAME(IsStream)==false){
            SYNC(spmvGpu);
        }else{
            SYNC(spmvGpuStream);
        }
    }else if(config == "hybrid"){

        //adjust csr based on the groups
        vector<GroupAttr> *groupsVec = &FRAME(groupsVec);
        vector<GroupCSR<double>> *gsCSR = &FRAME(gsCSR);
        for(int i=0; i<groupsVec->size();++i){
            int rstart = groupsVec->at(i).start;
            int rend   = groupsVec->at(i).end;
            int vcnum = rowDelimiters[rend]-rowDelimiters[rstart];

            GroupCSR<double> gcsr(rstart,rend,vcnum);
            adjustCSR<double>(gcsr, csrHost,rstart,rend);
            gsCSR->push_back(gcsr);
        }


        for(int i=0; i<groupsVec->size();++i){
            // invoke cpu group loop 
            SYNC(spmvCpuGroupLoop[i]); 
        }

        //========== configue hybrid ===========//
        
        int cpuNumRows  = FRAME(cpuNumRows);    
        int gpuNumRows  = FRAME(gpuNumRows);    
        int hostStart   = gpuNumRows;
        int deviceStart = 0;
        
        csrHost->assignStartPoint(hostStart);
        csrHost->assignNumRows(cpuNumRows);
        csrDevice->assignStartPoint(deviceStart);
        csrDevice->assignNumRows(gpuNumRows);
  
        __sync_fetch_and_sub(&FRAME(numRowsLeft),gpuNumRows); 
       
        if(FRAME(IsStream)==false){
            SYNC(spmvGpu);
        }else{
            SYNC(spmvGpuStream);
        }
        
        if(FRAME(avCnt)>0){

            __sync_fetch_and_sub(&FRAME(numRowsLeft),cpuNumRows); 
        }

        // invoke cpuloop if possible
        int nCnt = FRAME(nCnt); 
        for(int i=0; i<nCnt; ++i){
            SYNC(spmvCpuLoop[i]);
        }

    }else if(config == "cpu"){
        
        //int cpuNumRows = FRAME(cpuNumRows);    
        int cpuNumRows = FRAME(totalNumRows);
        csrHost->assignNumRows(cpuNumRows);

        __sync_fetch_and_sub(&FRAME(numRowsLeft),cpuNumRows); 
        for (size_t i=0; i<nCnt; ++i){
            SYNC(spmvCpuLoop[i]);
        }
    
    }
   

    SYNC(midSync);
   
#ifdef DARTS_RECORD
    int cpuNumRows = FRAME(cpuNumRows);    
    int cpuNumNonZeroes = rowDelimiters[cpuNumRows];
    RecordDatabase *recordDB = FRAME(recordDB);
    recordDB->getLastCpuRecord()->assignThreeParams("cpuloop", cpuNumRows,cpuNumNonZeroes);

#endif
    
    ttStart = getTime();

    EXIT_TP();
}


extern "C"
void SpmvCpuGroupLoopCD::fire(void){

#ifdef DARTS_RECORD
    uint64_t rdStart = getTime();
#endif
	LOAD_FRAME(SpmvTP);
	uint64_t Id = getID();
	RESET(spmvCpuGroupLoop[Id]);	

#ifdef DARTS_DEBUG

    pthread_mutex_lock(&mutex2);
	std::cout<<"invoke spmvCpuGroupLoop["<<Id<<"] !"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif

    CSRMM<double>*csrHost =&FRAME(csrHost);
    int numRows = csrHost->getNumRows();
   
    
    vector<GroupAttr> *groupsVec = &FRAME(groupsVec); 
    int start = groupsVec->at(Id).start;
    int end   = groupsVec->at(Id).end;
    int gId = groupsVec->at(Id).groupId;
    
    vector<GroupCSR<double>> *gsCSR = &FRAME(gsCSR);

    spmv_cpu_group_kernel(csrHost->getVec(),gsCSR->at(Id));

    SYNC(spmvCpuGroupLoopSync[gId]);
#ifdef DARTS_DEBUG
    std::cout<<"spmvCpuGroupLoop["<<Id<<"] sync groupLoopSync["<<gId<<"]"<<std::endl;
#endif

#ifdef DARTS_RECORD
    uint64_t rdEnd = getTime();
    RecordDatabase *recordDB = FRAME(recordDB);
    Record *rd = recordDB->getLastCpuRecord();
    rd->calcAndAddToValues(rdStart,rdEnd); 

#endif

#ifdef DARTS_RECORD
    recordDB->addToAllProc(rdEnd-rdStart);
#endif

#ifdef DARTS_DEBUG_CE
    pthread_mutex_lock(&mutex2);
	std::cout<<"spmvCpuGroupLoop["<<Id<<"] !"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif


    EXIT_TP();
}


void SpmvCpuGroupLoopSyncCD::fire(void)
{
    LOAD_FRAME(SpmvTP);

	uint64_t Id = getID();

#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
    std::cout<<"invoke spmvCpuGroupLoopSync["<<Id<<"]!"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif
    
    RESET(spmvCpuGroupLoopSync[Id]);	

    vector<pair<int,int>> *glsVec = &FRAME(glsVec);
    int nCnt = FRAME(nCnt);
    int sz = glsVec->at(Id).second; 
    
    __sync_fetch_and_sub(&FRAME(glCnt),sz); 
    if(__sync_bool_compare_and_swap(&FRAME(IsBlockCpuLoop), FRAME(glCnt)<nCnt, false )){
        SYNC(spmvCpuLoopSync);
#ifdef DARTS_DEBUG
	    pthread_mutex_lock(&mutex2);
        std::cout <<"SpmvCpuLoopGroupSync "<<Id<<" invoke SpmvCpuLoopSync, glCnt: "<<FRAME(glCnt)<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif
    }   
    SYNC(midSync);

    EXIT_TP();
}



extern "C"
void SpmvCpuLoopCD::fire(void){

#ifdef DARTS_RECORD
    uint64_t rdStart = getTime();
#endif

    LOAD_FRAME(SpmvTP);
	uint64_t Id = getID();
	RESET(spmvCpuLoop[Id]);	

#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"invoke spmvCpuLoop["<<Id<<"]!"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif

    CSRMM<double>*csrHost =&FRAME(csrHost);
    int numRows = csrHost->getNumRows();
   
    
    
    int avCnt = FRAME(avCnt);
    int chunk = numRows/avCnt;
    int secStart = csrHost->getStartPoint();
    int start = chunk*Id ;
    int end   = (Id ==(avCnt-1))?numRows:(chunk*(Id+1)) ;

    start +=secStart;
    end   +=secStart;

    spmv_cpu_kernel(csrHost,start,end);

    SYNC(spmvCpuLoopSync);


#ifdef DARTS_RECORD
    uint64_t rdEnd = getTime();
    RecordDatabase *recordDB = FRAME(recordDB);
    //Record *rd = recordDB->getLastCpuRecord();
    //rd->calcAndAddToValues(rdStart,rdEnd); 

#endif

#ifdef DARTS_RECORD
    //uint64_t suStart = getTime();
    //IntArrayStats rowInfo("cpuloop start at "+to_string(start),FRAME(nonZeroesPerRow));
    //rowInfo.calcArrayStats(start,end);
    //recordDB->addArrayStats(rowInfo);
    //uint64_t suEnd = getTime();

    //recordDB->addTime(rdEnd-rdStart,suEnd-suStart);
    recordDB->addToAllProc(rdEnd-rdStart);
#endif

#ifdef DARTS_DEBUG_CE
    pthread_mutex_lock(&mutex2);
	std::cout<<"spmvCpuLoop["<<Id<<"]!"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif

    EXIT_TP();
}


void
SpmvCpuLoopSyncCD::fire(void)
{
#ifdef DARTS_DEBUG
	
    pthread_mutex_lock(&mutex2);
    std::cout<<"invoke spmvCpuLoopSync!"<<std::endl;
    pthread_mutex_unlock(&mutex2);

#endif
    LOAD_FRAME(SpmvTP);
    RESET(spmvCpuLoopSync);

    CSRMM<double>*csrHost   = &FRAME(csrHost);
    CSRMM<double>*csrDevice = &FRAME(csrDevice);
    RecordDatabase *recordDB = FRAME(recordDB);

#ifdef DARTS_RECORD
    Record *cpuRd = recordDB->getLastCpuRecord();
    cpuRd->value = cpuRd->calcMax() ;
    recordDB->addRecord("cpu");
    recordDB->cpuCntPlus1();
    //cpuRd->print();
#endif

    string config = FRAME(config);
    bool IsNext = false;
   
    int sz = csrHost->getNumRows();
    
    if(config == "cpu"){
        int numRowsLeft = FRAME( numRowsLeft      );
        if(numRowsLeft == 0){
            SYNC(midSync);
        }else{
            
           IsNext = calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft), csrHost,csrDevice, "cpu");

        }

    }else{
        //====(config == "hybrid")=========//

        
#ifdef DARTS_DEBUG
        pthread_mutex_lock(&mutex2);
        std::cout<<"cpu sync before: numRowsLeft = "<< FRAME(numRowsLeft)<<", numNonZeroesLeft: "<<FRAME(numNonZeroesLeft)<<", cpu rowStartPoint: "<<csrHost->getStartPoint()<<std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
   
        int nCnt = FRAME(nCnt);
        int avCnt = FRAME(avCnt);


#ifdef DARTS_DEBUG
	    pthread_mutex_lock(&mutex2);
        std::cout <<"SpmvCpuLoopSync glCnt: "<<FRAME(glCnt)<<", avCnt:"<<FRAME(avCnt)<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif
        if(avCnt < nCnt){
            avCnt = nCnt - FRAME(glCnt);
            FRAME(avCnt) = avCnt;
        }

#ifdef DARTS_DEBUG
	    pthread_mutex_lock(&mutex2);
        std::cout <<"SpmvCpuLoopSync glCnt: "<<FRAME(glCnt)<<", avCnt:"<<FRAME(avCnt)<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif

        pthread_mutex_lock(&mutex1);
        __sync_synchronize();
        if(FRAME(numRowsLeft) == 0){
            SYNC(midSync);
            IsNext = false;
        }else{

            IsNext = calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft), csrHost,csrDevice, "cpu");

        }
        
        __sync_synchronize();
        pthread_mutex_unlock(&mutex1);
        
    
    }

#ifdef DARTS_RECORD
    cpuRd->clearRecord();
    cpuRd->assignThreeParams("cpuloop", csrHost->getNumRows(),csrHost->getNumNonZeroes());
#endif


    if(IsNext == true){
        int avCnt = FRAME(avCnt);
        FRAME(spmvCpuLoopSync).getSyncSlot()->initSyncSlot(avCnt,avCnt);
        for (size_t i=0; i<avCnt; ++i){
            SYNC(spmvCpuLoop[i]);
        }
    }
    
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
    std::cout<<"cpu sync after: numRowsLeft = "<< FRAME(numRowsLeft)<<", numNonZeroesLeft: "<<FRAME(numNonZeroesLeft)<<", cpu rowStartPoint: "<<csrHost->getStartPoint()<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif

//#ifdef DARTS_RECORD
//    int cpuNumRows = FRAME(cpuNumRows);    
//    int cpuNumNonZeroes = rowDelimiters[cpuNumRows];
//    Record rd("cpuloop", cpuNumRows,cpuNumNonZeroes);
//    RecordDatabase *recordDB = FRAME(recordDB);
//    recordDB->assignLastCpuRecord(rd);
//
//#endif


    EXIT_TP();
}

extern "C"
void SpmvGpuCD::fire(void){
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"Invoke SpmvGpu!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif
	LOAD_FRAME(SpmvTP);
	RESET(spmvGpu);
    
    RecordDatabase *recordDB = FRAME(recordDB);
    ResultDatabase *resultDB = FRAME(resultDB);
    OptionParser *op = FRAME(op);
    CSRMM<double>*csrHost   = &FRAME(csrHost);
    CSRMM<double>*csrDevice = &FRAME(csrDevice);

    
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
    std::cout<<"startPoint: "<<csrDevice->getStartPoint()<<", numRows: "<<csrDevice->getNumRows()<<", numNonZeroes: "<<csrDevice->getNumNonZeroes()<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif
   
#ifdef DARTS_RECORD
    uint64_t rdStart = getTime();
#endif


    csrTestScalar(resultDB,op,csrHost,csrDevice); 


#ifdef DARTS_RECORD
    uint64_t rdEnd = getTime();
#endif


#ifdef DARTS_DEBUG_VAL  
    int startPoint = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+startPoint;
    double *h_val = csrHost->getVal()+h_rowDelimiters[0];
    int *h_cols= csrHost->getCols()+h_rowDelimiters[0];
    double *h_out = csrHost->getOut()+startPoint;
    std::cout<<"h_val addr = "<<h_val<<std::endl;
    std::cout<<"h_cols addr = "<<h_cols<<std::endl;
    std::cout<<"h_out addr = "<<h_out<<std::endl;
    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_val["<<i<<"] = "<<h_val[i]<<std::endl;
    }

    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_cols["<<i<<"] = "<<h_cols[i]<<std::endl;
    }

    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_rowDelimiters["<<i<<"] = "<<h_rowDelimiters[i]<<std::endl;
    }
    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_out["<<i<<"] = "<<h_out[i]<<std::endl;
    }
#endif

    SYNC(spmvGpuCheck);

#ifdef DARTS_RECORD
    recordDB->addToAllProc(rdEnd-rdStart);
    std::cout<<"gpu time: "<<rdEnd-rdStart<<std::endl;
#endif

#ifdef DARTS_DEBUG_CE
    pthread_mutex_lock(&mutex2);
	std::cout<<"SpmvGpu!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif

    EXIT_TP();

}

extern "C"
void SpmvGpuCheckCD::fire(void){
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"Invoke SpmvGpuCheck!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif
	LOAD_FRAME(SpmvTP);
	RESET(spmvGpuCheck);

    CSRMM<double>*csrHost   = &FRAME(csrHost);
    CSRMM<double>*csrDevice = &FRAME(csrDevice);
   
    int sz = csrDevice->getNumRows();

    string config = FRAME(config);
    if(config == "gpu" ){

        int numRowsLeft         =FRAME( numRowsLeft      );
        if(numRowsLeft == 0){
            SYNC(midSync);
            EXIT_TP();
        }
        //int totalNumRows        =FRAME( totalNumRows     );
        //int nRows = (numRowsLeft > FRAME(gpuNumRows))?FRAME(gpuNumRows): numRowsLeft;
        //int newStart = totalNumRows - numRowsLeft;

#ifdef DARTS_DEBUG
        
        int oldStart = csrDevice->getStartPoint();
        int oldNumRows = csrDevice->getNumRows();
        std::cout<<"totalNumRows: "<<FRAME(totalNumRows)<<std::endl;
        std::cout<<"numRowsLeft: "<<FRAME(numRowsLeft)<<std::endl;
        std::cout<<"old start: "<< oldStart<<", oldNumRows: "<<oldNumRows<<std::endl;
#endif
        //csrDevice->assignStartPoint(newStart);
        //int *rowDelimiters = csrHost->getRowDelimiters();
        //csrDevice->assignNumRows(nRows);
        //
        //__sync_fetch_and_sub(&FRAME(numRowsLeft),nRows);

        calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft), csrHost,csrDevice, "gpu");

#ifdef DARTS_DEBUG
        std::cout<<"device: new start: "<< csrDevice->getStartPoint()<<", newNumRows: "<<csrDevice->getNumRows()<<std::endl;
        std::cout<<"numRowsLeft: "<<FRAME(numRowsLeft)<<std::endl;
#endif
   
        SYNC(spmvGpu);

    }else {
        //=========(config == "hybrid")=======//
        bool IsNext = false;
        pthread_mutex_lock(&mutex1);
        __sync_synchronize();
       
#ifdef DARTS_DEBUG
        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu sync before: numRowsLeft = "<< FRAME(numRowsLeft)<<", numNonZeroesLeft: "<<FRAME(numNonZeroesLeft)<<", gpu rowStartPoint: "<<csrDevice->getStartPoint()<< std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
        if(FRAME(numRowsLeft) == 0){
            SYNC(midSync);
            IsNext = false;
        }else{
            IsNext = calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft),csrHost,csrDevice, "gpu");
            
        }
        
        __sync_synchronize();
        pthread_mutex_unlock(&mutex1);
        
#ifdef DARTS_DEBUG
        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu IsNext = "<<IsNext<<std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
        if (IsNext == true){
            SYNC(spmvGpu);
        }

#ifdef DARTS_DEBUG

        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu sync after: numRowsLeft = "<< FRAME(numRowsLeft)<<", gpu rowStartPoint: "<<csrDevice->getStartPoint()<< std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
    }
   
    EXIT_TP();

}


extern "C"
void SpmvGpuStreamCD::fire(void){
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"Invoke SpmvGpuStream!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif
	LOAD_FRAME(SpmvTP);
	RESET(spmvGpuStream);
    
    RecordDatabase *recordDB = FRAME(recordDB);
    ResultDatabase *resultDB = FRAME(resultDB);
    OptionParser *op = FRAME(op);
    CSRMM<double>*csrHost   = &FRAME(csrHost);
    CSRMM<double>*csrDevice = &FRAME(csrDevice);


#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
    std::cout<<"startPoint: "<<csrDevice->getStartPoint()<<", numRows: "<<csrDevice->getNumRows()<<", numNonZeroes: "<<csrDevice->getNumNonZeroes()<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif
   
#ifdef DARTS_RECORD
    uint64_t rdStart = getTime();
#endif

    //csrTestScalar(resultDB,op,csrHost,csrDevice); 
    csrStreamTestScalar(resultDB,op,csrHost,csrDevice); 


#ifdef DARTS_RECORD
    uint64_t rdEnd = getTime();
#endif


#ifdef DARTS_DEBUG_VAL  
    int startPoint = csrDevice->getStartPoint();
    int *h_rowDelimiters = csrHost->getRowDelimiters()+startPoint;
    double *h_val = csrHost->getVal()+h_rowDelimiters[0];
    int *h_cols= csrHost->getCols()+h_rowDelimiters[0];
    double *h_out = csrHost->getOut()+startPoint;
    std::cout<<"h_val addr = "<<h_val<<std::endl;
    std::cout<<"h_cols addr = "<<h_cols<<std::endl;
    std::cout<<"h_out addr = "<<h_out<<std::endl;
    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_val["<<i<<"] = "<<h_val[i]<<std::endl;
    }

    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_cols["<<i<<"] = "<<h_cols[i]<<std::endl;
    }

    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_rowDelimiters["<<i<<"] = "<<h_rowDelimiters[i]<<std::endl;
    }
    for(size_t i=0; i< 20; ++i){
        std::cout<<"h_out["<<i<<"] = "<<h_out[i]<<std::endl;
    }
#endif

    SYNC(spmvGpuStreamCheck);


#ifdef DARTS_RECORD
    recordDB->addToAllProc(rdEnd-rdStart);
    std::cout<<"gpu stream time: "<<rdEnd-rdStart<<std::endl;
#endif

#ifdef DARTS_DEBUG_CE
    pthread_mutex_lock(&mutex2);
	std::cout<<"SpmvGpuStream!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif

    EXIT_TP();

}

extern "C"
void SpmvGpuStreamCheckCD::fire(void){
#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"Invoke SpmvGpuStreamCheck!"<<std::endl;	
    pthread_mutex_unlock(&mutex2);
#endif
	LOAD_FRAME(SpmvTP);
	RESET(spmvGpuStreamCheck);

    CSRMM<double>*csrHost   = &FRAME(csrHost);
    CSRMM<double>*csrDevice = &FRAME(csrDevice);
   
    int sz = csrDevice->getNumRows();

    string config = FRAME(config);
    if(config == "gpu" ){

        int numRowsLeft         =FRAME( numRowsLeft      );
        if(numRowsLeft == 0){
            SYNC(midSync);
            EXIT_TP();
        }
        //int totalNumRows        =FRAME( totalNumRows     );
        //int nRows = (numRowsLeft > FRAME(gpuNumRows))?FRAME(gpuNumRows): numRowsLeft;
        //int newStart = totalNumRows - numRowsLeft;


#ifdef DARTS_DEBUG
        
        int oldStart = csrDevice->getStartPoint();
        int oldNumRows = csrDevice->getNumRows();
        std::cout<<"totalNumRows: "<<FRAME(totalNumRows)<<std::endl;
        std::cout<<"numRowsLeft: "<<FRAME(numRowsLeft)<<std::endl;
        std::cout<<"old start: "<< oldStart<<", oldNumRows: "<<oldNumRows<<std::endl;
#endif
        //csrDevice->assignStartPoint(newStart);
        //int *rowDelimiters = csrHost->getRowDelimiters();
        //csrDevice->assignNumRows(nRows);
        //
        //__sync_fetch_and_sub(&FRAME(numRowsLeft),nRows);

        calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft), csrHost,csrDevice, "gpu");

#ifdef DARTS_DEBUG
        std::cout<<"device: new start: "<< csrDevice->getStartPoint()<<", newNumRows: "<<csrDevice->getNumRows()<<std::endl;
        std::cout<<"numRowsLeft: "<<FRAME(numRowsLeft)<<std::endl;
#endif
   
        SYNC(spmvGpuStream);

    }else {
        //=========(config == "hybrid")=======//
        bool IsNext = false;
        pthread_mutex_lock(&mutex1);
        __sync_synchronize();
       
#ifdef DARTS_DEBUG
        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu sync before: numRowsLeft = "<< FRAME(numRowsLeft)<<", numNonZeroesLeft: "<<FRAME(numNonZeroesLeft)<<", gpu rowStartPoint: "<<csrDevice->getStartPoint()<< std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
        if(FRAME(numRowsLeft) == 0){
            SYNC(midSync);
            IsNext = false;
        }else{
            IsNext = calcStaticNext(FRAME(totalNumRows), FRAME(numRowsLeft),csrHost,csrDevice, "gpu");
            
        }
        
        __sync_synchronize();
        pthread_mutex_unlock(&mutex1);
        
#ifdef DARTS_DEBUG
        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu IsNext = "<<IsNext<<std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
        if (IsNext == true){
            SYNC(spmvGpuStream);
        }

#ifdef DARTS_DEBUG

        pthread_mutex_lock(&mutex2);
        std::cout<<"gpu sync after: numRowsLeft = "<< FRAME(numRowsLeft)<<", gpu rowStartPoint: "<<csrDevice->getStartPoint()<< std::endl;
        pthread_mutex_unlock(&mutex2);
#endif
    }
   
    EXIT_TP();

}


void
MidSyncCD::fire(void)
{
    
    ttEnd = getTime();
    std::cout<<"total time: "<<ttEnd-ttStart<<std::endl;

#ifdef DARTS_DEBUG
    pthread_mutex_lock(&mutex2);
	std::cout<<"invoke midSync!"<<std::endl;
    pthread_mutex_unlock(&mutex2);
#endif
	LOAD_FRAME(SpmvTP);
    RESET(midSync);
   
    string config = FRAME(config);

    CSRMM<double>*csrHost =&FRAME(csrHost);
    double *refOut = FRAME(refOut);
    csrHost->assignStartPoint(0);
    double *out = csrHost->getOut();
    int numRows = FRAME(totalNumRows);
    vector<GroupCSR<double>> *gsCSR = &FRAME(gsCSR);
    if(config == "hybrid"){
        addGroupResultsBack(out, gsCSR);
    }

    verifyResults(refOut, out,numRows);

    csrHost->freeMem();

    
    CSRMM<double>*csrDevice =&FRAME(csrDevice);
    if(csrDevice->getKind()=="device"){
        csrDevice->freeMem();
    }

    SYNC(sync);
    EXIT_TP();
}

void
SyncCD::fire(void)
{
#ifdef DARTS_DEBUG
	std::cout<<"invoke Sync!"<<std::endl;
#endif
	LOAD_FRAME(SpmvTP);
	SIGNAL(signalUp);
    EXIT_TP();
}
