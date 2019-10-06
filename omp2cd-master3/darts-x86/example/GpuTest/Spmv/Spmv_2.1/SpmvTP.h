
#ifndef SPMVTP_H
#define SPMVTP_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

#include "OptionParser.h"
#include "ResultDatabase.h"
#include "RecordDatabase.h"
#include "Record.h"
#include "util.h"
#include "CSRMM.h"
#include "conf.h"
#include "Spmv.h"
#include "ArrayStats.h"
#include "GroupCSR.h"

#include "DARTS.h"


using namespace darts;
using namespace std;

#define GPUMETA 0x4

#define N_CU_CORES TOTAL_NUM_CU

DEF_CODELET(SpmvCpuCD,2,LONGWAIT);
DEF_CODELET(SpmvCpuSyncCD,2,LONGWAIT);
DEF_CODELET_ITER(SpmvCpuLoopCD,0,SHORTWAIT);
DEF_CODELET(SpmvCpuLoopSyncCD,0,SHORTWAIT);

DEF_CODELET_ITER(SpmvCpuGroupLoopCD,0,SHORTWAIT);
DEF_CODELET_ITER(SpmvCpuGroupLoopSyncCD,0,SHORTWAIT);


DEF_CODELET(SpmvGpuSHOCCsrCD,0,SHORTWAIT);
DEF_CODELET(SpmvGpuSHOCCsrCheckCD,0,SHORTWAIT);

DEF_CODELET(SpmvGpuSHOCCsrStreamCD,0,SHORTWAIT);
DEF_CODELET(SpmvGpuSHOCCsrStreamCheckCD,0,SHORTWAIT);

DEF_CODELET(SpmvGpuCuSparseCsrCD,0,SHORTWAIT);
DEF_CODELET(SpmvGpuCuSparseCsrCheckCD,0,SHORTWAIT);

DEF_CODELET(SpmvGpuCuSparseCsrStreamCD,0,SHORTWAIT);
DEF_CODELET(SpmvGpuCuSparseCsrStreamCheckCD,0,SHORTWAIT);


DEF_CODELET(MidSyncCD,2,LONGWAIT);
DEF_CODELET(SyncCD,2,LONGWAIT);

DEF_TP(SpmvTP)
{
    RecordDatabase *recordDB;
    ResultDatabase *resultDB;
    OptionParser *op;

    string config;
    CSRMM<double> csrHost;
    CSRMM<double> csrDevice;
    double *refOut;
    double *h_vec;
    int nCnt;
    int glCnt;
    int avCnt;
    bool IsRefReq;
    bool IsCpuParallel;
    bool IsGpuSplit;
    bool IsBlockCpuLoop;
    bool IsStream;
    bool IsCuSparse;
    bool IsSHOC;
    int totalNumRows;
    int totalNumNonZeroes;
    int numRowsLeft;
    int numNonZeroesLeft;
    int cpuNumRows;
    int gpuNumRows;
    int *nonZeroesPerRow;
    //vector<pair<int,int>>nonZeroesPerRowVec;
    //vector<vector<pair<int,int>>>groupsVec;
    vector<GroupAttr>groupsVec;
    vector<pair<int,int>>glsVec;
    vector<GroupCSR<double>> gsCSR;
    SpmvCpuCD   spmvCpu;
    SpmvCpuSyncCD   spmvCpuSync;
    MidSyncCD   midSync;
    SpmvCpuLoopCD *spmvCpuLoop = NULL;
    SpmvCpuLoopSyncCD spmvCpuLoopSync;
    
    SpmvCpuGroupLoopCD *spmvCpuGroupLoop = NULL;
    SpmvCpuGroupLoopSyncCD *spmvCpuGroupLoopSync=NULL;
    
    SpmvGpuSHOCCsrCD spmvGpuSHOCCsr;
    SpmvGpuSHOCCsrCheckCD spmvGpuSHOCCsrCheck;
    SpmvGpuSHOCCsrStreamCD spmvGpuSHOCCsrStream;
    SpmvGpuSHOCCsrStreamCheckCD spmvGpuSHOCCsrStreamCheck;
    
    SpmvGpuCuSparseCsrCD spmvGpuCuSparseCsr;
    SpmvGpuCuSparseCsrCheckCD spmvGpuCuSparseCsrCheck;
   
    SpmvGpuCuSparseCsrStreamCD spmvGpuCuSparseCsrStream;
    SpmvGpuCuSparseCsrStreamCheckCD spmvGpuCuSparseCsrStreamCheck;

    SyncCD	sync;
    Codelet *signalUp;


    SpmvTP(RecordDatabase *recordDB,ResultDatabase *resultDB, OptionParser *op,Codelet *up)
        :recordDB(recordDB)
        ,resultDB(resultDB)
        ,op(op)
        ,csrHost("host")
        //,spmvCpu(0,1, this, LONGWAIT)
	    ,sync(1,1,this,LONGWAIT)
        ,signalUp(up)
        {
#ifdef DARTS_DEBUG
            std::cout<<"SpmvTP begin!"<<std::endl;
#endif 
            //op->print();    
            // This benchmark either reads in a matrix market input file or generates a random matrix
            string inFileName = op->getOptionString("mm_filename");
            
            config = op->getOptionString("config");
            csrHost.assignConfig(config);

            if (inFileName == "random"){

            }else{
                char filename[FIELD_LENGTH];
                strcpy(filename, inFileName.c_str());
                readMarketMatrixToCSR<double>(filename, &csrHost); 
            }

#ifdef DARTS_DEBUG
            csrHost.print();        
#endif
            int numRows = csrHost.getNumRows();
            int numCols = csrHost.getNumCols();
            totalNumRows = csrHost.getNumRows();
            totalNumNonZeroes = csrHost.getNumNonZeroes();
           
            refOut = new double[numRows];
            int *rowDelimiters = csrHost.getRowDelimiters();

#ifdef DARTS_DEBUG
            //std::cout<<"process nonZeroesPerRowVec"<<std::endl;
#endif
            
            //nonZeroesPerRow = new int[numRows];
            //preProcessNonZeroes(nonZeroesPerRowVec,nonZeroesPerRow,rowDelimiters,0,numRows);
            //calcGroups(rowDelimiters,HroupsVec, nonZeroesPerRowVec,100, 0.8 ); 
            //vector<pair<int,int>> refcsv;
            //refcsv.push_back(make_pair(4641,4642));
            //refcsv.push_back(make_pair(4642,4643));
            //refcsv.push_back(make_pair(4643,4644));
            //for(int i=0; i<refcsv.size();++i){
            //    int rst = refcsv[i].first;
            //    int red = refcsv[i].second;
            //    //string csvr = "refr"+to_string(rst)+".csv";
            //    //string csvv = "refv"+to_string(rst)+".csv";
            //    //string csvc = "refc"+to_string(rst)+".csv";
            //    string csvall = "refall"+to_string(rst)+".csv";
            //    
            //    ofstream out;

            //    out.open(csvall, std::ofstream::out);

            //    out<<"rowDelimiters:\n";
            //    for(int i=rowDelimiters[rst]; i<rowDelimiters[rst]+20; ++i){
            //        out<<i<<",";
            //    }
            //    out<<"\n";

            //    out<<"cols:\n";
            //    for(int i=rowDelimiters[rst]; i<rowDelimiters[rst]+20; ++i){
            //        out<<csrHost.getCols()[i]<<",";
            //    }
            //    out<<"\n";

            //    out<<"val:\n";
            //    for(int i=rowDelimiters[rst]; i<rowDelimiters[rst]+20; ++i){
            //        out<<csrHost.getVal()[i]<<",";
            //    }
            //    out<<endl;
            //    out.close();
            //    //out.open(csvr, std::ofstream::out);
            //    ////out<<"rowDelimiters:\n";
            //    //for(int i=rowDelimiters[rst]; i<rowDelimiters[red]; ++i){
            //    //    out<<i<<",";
            //    //}
            //    //out<<endl;
            //    //out.close();
            //    //
            //    //out.open(csvv, std::ofstream::out);
            //    ////out<<"val:\n";
            //    //for(int i=rowDelimiters[rst]; i<rowDelimiters[red]; ++i){
            //    //    out<<csrHost.getVal()[i]<<"\n";
            //    //}
            //    //out<<endl;
            //    //out.close();
            //    
            //    //out.open(csvc, std::ofstream::out);
            //    //////out<<"cols: \n";
            //    //for(int i=rowDelimiters[rst]; i<rowDelimiters[red]; ++i){
            //    //    out<<csrHost.getCols()[i]<<"\n";
            //    //}
            //    //out<<endl;
            //    //out.close();
            //}

            //============= initial vec =============// 
            h_vec = csrHost.getVec(); 
            fill(h_vec, numRows, op->getOptionFloat("maxval"));

            //============= Compute reference solution =============// 
            //spmvCPU(&csrHost);
            //double *out = csrHost.getOut();
            //int numRows = csrHost.getNumRows();
            //memcpy(refOut,out, (numRows)*sizeof(double));
            //csrHost->Memset<double>(out,0,numRows);
            
            //============= initial csrDevice =============// 
            if(config != "cpu"){
                int numNonZeroes = csrHost.getNumNonZeroes();
                csrDevice.assignKind("device");
                csrDevice.assignConfig(config);
                csrDevice.assignNumNonZeroes(numNonZeroes);
                csrDevice.assignNumRows(numRows);
                csrDevice.assignNumCols(numCols);
                csrDevice.allocateDevice();

                double *d_vec = csrDevice.getVec();
                memcpyHostToDevice(d_vec,csrHost.getVec(),numRows);
                memcpyDeviceTexture<double>(d_vec,numRows);
            }

            //op->print();



            IsRefReq = true;
            IsCpuParallel = true;
            IsGpuSplit = true;
           
            int dep = (IsRefReq==true)?1:0; 
            nCnt = (IsCpuParallel == true)?N_CU_CORES:0;
            
            IsStream   = true;
            IsSHOC     = false; 
            IsCuSparse = not(IsSHOC);



            //============= reference codelet =============// 
#ifdef DARTS_DEBUG
            std::cout<<"add spmvCpu!"<<std::endl;
#endif
            spmvCpu = SpmvCpuCD{0,1,this,SHORTWAIT};
            add(&spmvCpu);
#ifdef DARTS_DEBUG
            std::cout<<"add spmvCpuSync!"<<std::endl;
#endif
            spmvCpuSync = SpmvCpuSyncCD{1,1,this,SHORTWAIT};


            //========= gpu codelet (gpu and hybrid)======// 
            if(config != "cpu"){

                if(IsSHOC == true){
                    if(IsStream == false){
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuSHOCCsr!"<<std::endl;
#endif
                        spmvGpuSHOCCsr = SpmvGpuSHOCCsrCD{dep,1,this,GPUMETA};
                        add(&spmvGpuSHOCCsr);
                
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuSHOCCsrCheck!"<<std::endl;
#endif
                        spmvGpuSHOCCsrCheck = SpmvGpuSHOCCsrCheckCD{1,1,this,GPUMETA};
                    }else{

#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuSHOCCsrStream!"<<std::endl;
#endif
                        spmvGpuSHOCCsrStream = SpmvGpuSHOCCsrStreamCD{dep,1,this,GPUMETA};
                        add(&spmvGpuSHOCCsrStream);
                
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuSHOCCsrStreamCheck!"<<std::endl;
#endif
                        spmvGpuSHOCCsrStreamCheck = SpmvGpuSHOCCsrStreamCheckCD{1,1,this,GPUMETA};

                    }
                }else if(IsCuSparse == true){

                    if(IsStream == false){
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuCuSparseCsr!"<<std::endl;
#endif
                        spmvGpuCuSparseCsr = SpmvGpuCuSparseCsrCD{dep,1,this,GPUMETA};
                        add(&spmvGpuCuSparseCsr);
                
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuCuSparseCsrCheck!"<<std::endl;
#endif
                        spmvGpuCuSparseCsrCheck = SpmvGpuCuSparseCsrCheckCD{1,1,this,GPUMETA};
                    }else{

#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuCuSparseCsrStream!"<<std::endl;
#endif
                        spmvGpuCuSparseCsrStream = SpmvGpuCuSparseCsrStreamCD{dep,1,this,GPUMETA};
                        add(&spmvGpuCuSparseCsrStream);
                
#ifdef DARTS_DEBUG
                        std::cout<<"add spmvGpuCuSparseCsrStreamCheck!"<<std::endl;
#endif
                        spmvGpuCuSparseCsrStreamCheck = SpmvGpuCuSparseCsrStreamCheckCD{1,1,this,GPUMETA};

                    }

                }
            }


            //=============hybrid  case : cpu  =============//
            int glsNum = 0;
            if( config == "hybrid"){
#ifdef DARTS_DEBUG
                std::cout<<"build groups!"<<std::endl;
#endif
                //glsNum = buildGroups(groupsVec,rowDelimiters, 0, numRows, 50, 0.8);
                glsNum = buildGroups_v2(groupsVec,rowDelimiters, 0, numRows, 50, 0.7,100);

                int glNum = groupsVec.size(); 
#ifdef DARTS_DEBUG
                std::cout<<"group loop Num = "<<glNum<<",group loop sync Num: "<<glsNum<<std::endl;
#endif
                spmvCpuGroupLoop = new SpmvCpuGroupLoopCD[glNum]; 
                spmvCpuGroupLoopSync = new SpmvCpuGroupLoopSyncCD[glsNum]; 
                glCnt = glNum; 
                
                int sc = 0;
                int pgId = groupsVec[0].groupId;
                for(int i=0; i<glNum; ++i){
                    spmvCpuGroupLoop[i] = SpmvCpuGroupLoopCD{dep,1,this,SHORTWAIT,i};
                    if(pgId == groupsVec[i].groupId){
                        ++sc;
                        if(i == glNum-1){
                            spmvCpuGroupLoopSync[pgId] = SpmvCpuGroupLoopSyncCD{sc,sc,this,SHORTWAIT,pgId};
                            glsVec.push_back(make_pair(pgId,sc));
#ifdef DARTS_DEBUG
                            std::cout<<"group Id: "<<pgId<<", dep: "<<sc<<std::endl; 
#endif
                        }
                    }else{
                        spmvCpuGroupLoopSync[pgId] = SpmvCpuGroupLoopSyncCD{sc,sc,this,SHORTWAIT,pgId};
                        glsVec.push_back(make_pair(pgId,sc));
#ifdef DARTS_DEBUG
                        std::cout<<"group Id: "<<pgId<<", dep: "<<sc<<std::endl; 
#endif
                        sc = 1;
                        pgId = groupsVec[i].groupId;
                    }
                }
           
                avCnt = nCnt - glCnt;
                spmvCpuLoop = new SpmvCpuLoopCD[nCnt]; 
                for (int i = 0; i<nCnt; ++i ){
                    spmvCpuLoop[i] = SpmvCpuLoopCD{dep+1,1,this, SHORTWAIT,i};
                    if(i<avCnt){
                        spmvCpuLoop[i].decDep();
                    }
                }
                IsBlockCpuLoop = (avCnt >0)? false: true;
                int tmpCnt = (avCnt>0)?avCnt: 1;
                spmvCpuLoopSync = SpmvCpuLoopSyncCD{tmpCnt,tmpCnt,this, SHORTWAIT};
#ifdef DARTS_DEBUG
                std::cout<<"nCnt: "<<nCnt<<",glCnt: "<<glCnt<<",avCnt: "<<avCnt<<std::endl;
#endif               

            }
            

            //============= pure cpu case =============// 
#ifdef DARTS_DEBUG
            std::cout<<"nCnt = "<<nCnt<<std::endl;
#endif            
            if((config == "cpu")&& (nCnt !=0)){
                avCnt = nCnt; 
                spmvCpuLoop = new  SpmvCpuLoopCD[avCnt];
                for (size_t i =0; i< avCnt; ++i){
                    spmvCpuLoop[i] = SpmvCpuLoopCD{dep,1,this, SHORTWAIT,i};
                    add(spmvCpuLoop +i);
                }
                spmvCpuLoopSync = SpmvCpuLoopSyncCD{avCnt,avCnt,this,SHORTWAIT};
            }

            
            
            //=============mid sync codelet =============// 
            //cpu: one from ref, one from loop
            //gpu: one from ref, one from gpu
            //hybrid: one from ref, one from cpu, one from gpu
            int depMidSync =0;
            if((config == "cpu")&&(IsCpuParallel==false )){
                depMidSync = 1;
            }else if(config == "hybrid" && IsRefReq == true){
                depMidSync = 3 + glsNum;
            }
            else{
                depMidSync = 2;
            }
            midSync = MidSyncCD{depMidSync,depMidSync,this,SHORTWAIT};


            //=========== set init val ================//
            cpuNumRows = 50000;
            gpuNumRows = 5000;
            numRowsLeft = totalNumRows; 
        
        
        }

    virtual ~SpmvTP(){

#ifdef DARTS_DEBUG
        std::cout<<"SpmvTP finish!"<<std::endl;
#endif 
        if(config == "cpu"){
            delete [] spmvCpuLoop; 
        }
        delete [] refOut;
        //delete [] nonZeroesPerRow;
        delete [] spmvCpuGroupLoop;
        delete [] spmvCpuGroupLoopSync; 
    }

};


#endif
