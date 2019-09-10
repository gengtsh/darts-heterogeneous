#ifndef TESTFGTP_H
#define TESTFGTP_H

#include <vector>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
//#include "tbb/concurrent_vector.h"

using namespace std;

#include "DARTS.h"

#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES

#define SHORTWAIT 0x4

DEF_CODELET_ITER(GroupCD,0,SHORTWAIT);
DEF_CODELET_ITER(GroupSyncCD,0,SHORTWAIT);

DEF_CODELET_ITER(CpuLoopCD,0,SHORTWAIT);
DEF_CODELET(CpuLoopSyncCD,0,SHORTWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(TestGroupsTP)
{
	uint64_t timeStep;
    
	int nCnt;
    int currGroupCnt;
    int avCnt;
    //vector<int> Idx;
    vector<vector<int>> groups;
    //tbb::concurrent_vector<int> idxSync;
    int *idxSync;
    GroupCD *group;
    GroupSyncCD *groupSync;
    CpuLoopCD *cpuLoop;
    CpuLoopSyncCD cpuLoopSync;
	SyncCD	sync;
    Codelet *signalUp;
    bool IsGroup;    

	TestGroupsTP( uint64_t ts, Codelet *up)
	:timeStep(ts)
    ,nCnt(N_CORES-1) //1 GPU core, (N_CORES-1) CPU cores
	,signalUp(up)
	{
        vector<int> Idx1;
        Idx1.push_back(1);
        Idx1.push_back(2);
        groups.push_back(Idx1);
        Idx1.clear();
        Idx1.push_back(3);
        Idx1.push_back(4);
        Idx1.push_back(5);
        groups.push_back(Idx1);
        Idx1.clear();
        for(int i = 6; i<20; ++i){
            Idx1.push_back(i);
        }
        groups.push_back(Idx1);
        Idx1.clear();
        
        for(int i = 20; i<25; ++i){
            Idx1.push_back(i);
        }
        groups.push_back(Idx1);

        Idx1.clear();

        idxSync = new int[100];
        int ii=0;
        for (int i=0; i< groups.size();++i){
            for (int j=0;j<groups[i].size();++j){
                std::cout<<groups[i][j]<<",";
                idxSync[ii]=i;
                ++ii;
            }
            std::cout<<std::endl;
        }
        group = new GroupCD[ii];
        groupSync = new GroupSyncCD[groups.size()];
        
        currGroupCnt = 0;
        for (int i=0; i< groups.size();++i){
            int tt = groups[i].size();
            std::cout<<"groups["<<i<<"]: "<<groups[i].size()<<std::endl; 
            for (int j=0;j<tt;++j){
                group[currGroupCnt] = GroupCD {0,1,this,SHORTWAIT,currGroupCnt};
                add(group + currGroupCnt);
            
                currGroupCnt++;
            }
            groupSync[i] = GroupSyncCD{tt,tt,this, SHORTWAIT,i };
            add(groupSync + i);
        
        }
        avCnt = nCnt-currGroupCnt;
        cpuLoop = new CpuLoopCD[nCnt];
        std::cout<<"avCnt: "<<avCnt<<std::endl; 
        for(int i=0; i<nCnt; ++i){
            cpuLoop[i] = CpuLoopCD{1,1,this, SHORTWAIT,i};
            if(i<avCnt){
                cpuLoop[i].decDep();
            }

        }
        
        IsGroup = (avCnt >0)? false: true; 
        int tmpCnt = (avCnt>0)?avCnt:1;
        cpuLoopSync = CpuLoopSyncCD{tmpCnt,tmpCnt,this,SHORTWAIT};
        
        std::cout<<"nCnt: "<<nCnt<<",currGroupCnt: "<<currGroupCnt<<std::endl;
    
        int sCnt = groups.size()+1;
        sync = SyncCD{sCnt,sCnt, this, SHORTWAIT};

	}
	
	virtual ~TestGroupsTP(){
		delete []group;
		delete []groupSync;
        delete []idxSync;
	}
};

#endif
