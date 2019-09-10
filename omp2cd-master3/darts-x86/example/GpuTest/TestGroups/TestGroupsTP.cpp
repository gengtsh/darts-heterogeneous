#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "TestGroupsTP.h"
#include <cassert>
#include <pthread.h>
#include <sstream>
#include <iostream>
#include <unistd.h>
#include <thread>
#include <chrono>

pthread_mutex_t mutex1;
pthread_mutex_t mutex2;

void 
GroupCD::fire(void) 
{
	LOAD_FRAME(TestGroupsTP);

	int64_t Id	= getID();	

	RESET(group[Id]);
	
	
////	pthread_mutex_lock(&mutex);
////	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<std::endl;
////	pthread_mutex_unlock(&mutex);
////	std::stringstream ss;
////	ss <<"Group"<<Id<<" invoked!"<<"\n";
////	std::cout <<ss.str();

//	pthread_mutex_lock(&mutex2);
    std::cout <<"Group "<<Id<<" invoked!"<<std::endl;
//	pthread_mutex_unlock(&mutex2);


	int64_t nCnt = FRAME(nCnt);

    if((Id == 0)||(Id == 1)){
        //sleep(5);
        //std::this_thread::sleep_for(std::chrono::seconds(5));
        for(int i=0; i< 100000;++i){
    
            for(int j=0; j< 10000;++j){

            }
        }
    }else if(Id == 2|| Id==3 || Id ==4){
        //sleep(2);
        //std::this_thread::sleep_for(std::chrono::seconds(2));
    
        for(int i=0; i< 10000;++i){

            for(int j=0; j< 10000;++j){

            }
        }
    }else if(Id == 6|| Id==7 ){
        //sleep(2);
        //std::this_thread::sleep_for(std::chrono::seconds(2));
    
        for(int i=0; i< 1000;++i){

            for(int j=0; j< 1000;++j){

            }
        }
    }


    int syncId = FRAME(idxSync)[Id];
	SYNC(groupSync[syncId]);

	//pthread_mutex_lock(&mutex1);
    std::cout <<"Group "<<Id<<" sync Id:  "<<syncId<<std::endl;
	//pthread_mutex_unlock(&mutex1);

	EXIT_TP();
}

void
GroupSyncCD::fire(void)
{
	LOAD_FRAME(TestGroupsTP);
	int64_t Id = getID();	
	
	RESET(groupSync[Id]);

////	pthread_mutex_lock(&mutex);
////	std::cout <<"S"<<Id<<",ts:"<<ts[Id]<<std::endl;
////	pthread_mutex_unlock(&mutex);
////	std::stringstream ss;
////	ss <<"GroupSync "<<Id<<" invoked!"<<"\n";
////	std::cout<<ss.str();

	//pthread_mutex_lock(&mutex2);
    std::cout <<"GroupSync "<<Id<<" invoked!"<<std::endl;
    //pthread_mutex_unlock(&mutex2);

   
    
    int nCnt = FRAME(nCnt);
    int currSZ = FRAME(groups)[Id].size();
	//pthread_mutex_lock(&mutex1);
    //__sync_synchronize(); 
    //FRAME(currGroupCnt) = FRAME(currGroupCnt)-currSZ;
    //if(FRAME(IsGroup) == true && FRAME(currGroupCnt <nCnt)){
    //  SYNC(cpuLoopSync);
    //}
    //__sync_synchronize();
    //pthread_mutex_unlock(&mutex1);
    
    __sync_fetch_and_sub(&FRAME(currGroupCnt),currSZ); 
    if(__sync_bool_compare_and_swap(&FRAME(IsGroup), FRAME(currGroupCnt)<nCnt, false )){
           SYNC(cpuLoopSync);

	    //pthread_mutex_lock(&mutex1);
        std::cout <<"GroupSync "<<Id<<" invoke cpuLoopSync, currGroupCnt: "<<FRAME(currGroupCnt)<<std::endl;
	    //pthread_mutex_unlock(&mutex1);
    }   

    

	//pthread_mutex_lock(&mutex1);
    std::cout <<"GroupSync "<<Id<<",currGroupCnt: "<<FRAME(currGroupCnt)<<std::endl;
	//pthread_mutex_unlock(&mutex1);

    
    SYNC(sync);


	EXIT_TP();
	
}


void 
CpuLoopCD::fire(void) 
{

	LOAD_FRAME(TestGroupsTP);
	int64_t Id	= getID();	

	RESET(cpuLoop[Id]);
	
////	pthread_mutex_lock(&mutex);
////	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<std::endl;
////	pthread_mutex_unlock(&mutex);
////	std::stringstream ss;
////	ss <<"cpuLoop"<<Id<<" invoked!"<<"\n";
////	std::cout <<ss.str();

//	pthread_mutex_lock(&mutex2);
    std::cout <<"cpuLoop "<<Id<<" invoked!"<<std::endl;
//	pthread_mutex_unlock(&mutex2);

    for(int i=0; i< 10000; ++i){
        for(int j=0; j<5000; ++j){

        }
    }

	int64_t nCnt = FRAME(nCnt);

    SYNC(cpuLoopSync);

	EXIT_TP();
}


void
CpuLoopSyncCD::fire(void)
{
    
//	pthread_mutex_lock(&mutex2);
    std::cout<<"cpuLoopSync!"<<std::endl;
//	pthread_mutex_unlock(&mutex2);
    LOAD_FRAME(TestGroupsTP);
    
    int nCnt = FRAME(nCnt);
    int avCnt = FRAME(avCnt);

//	pthread_mutex_lock(&mutex2);
    std::cout<<"cpuLoopSync: avCnt: "<<avCnt<<", currGroupCnt: "<<FRAME(currGroupCnt)<<std::endl;
//	pthread_mutex_unlock(&mutex2);
    
    
    if(avCnt < nCnt ){
        avCnt = nCnt - FRAME(currGroupCnt);
        FRAME(avCnt) = avCnt;
    }
//	pthread_mutex_lock(&mutex2);
    std::cout<<"cpuLoopSync: avCnt: "<<avCnt<<", currGroupCnt: "<<FRAME(currGroupCnt)<<std::endl;
//	pthread_mutex_unlock(&mutex2);
    --FRAME(timeStep);
    if(FRAME(timeStep) == 0){
        SYNC(sync);
    }else{
        FRAME(cpuLoopSync).getSyncSlot()->initSyncSlot(avCnt,avCnt);
        RESET(cpuLoopSync);
        std::cout<<" invoke "<<avCnt<<" cpuLoop!"<<std::endl;
        for(int i=0; i<avCnt; ++i){
            SYNC(cpuLoop[i]);
        } 
    }


    EXIT_TP();
}

void
SyncCD::fire(void)
{
    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(TestGroupsTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





