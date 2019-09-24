#ifndef RECORDDATABASE_H_
#define RECORDDATABASE_H_


#include <cfloat>
#include <unistd.h>
#include <climits>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>

#include <string>
#include <vector>
#include "tbb/concurrent_vector.h"

#include "Record.h"
#include "ArrayStats.h"

using namespace std;

typedef ArrayStats<int> IntArrayStats;

class RecordDatabase
{
    
    public:
        RecordDatabase()
        {
            cpuCnt = 0;
            gpuCnt = 0;
        }
    
        void assignRef(Record &rd);
        void assignLastCpuRecord(Record &rd);
        void assignLastGpuRecord(Record &rd);
        void cpuCntPlus1(void);
        void gpuCntPlus1(void);
        int  getCpuCnt(void) const;
        int  getGpuCnt(void) const;
        Record* getLastCpuRecord(void) const;
        Record* getLastGpuRecord(void) const;
        
        void addRecord(string config);

        void clearAllRecords();
        //void dumpCsv(string fileName);
        
        void dumpCpuRecords(ostream&);
        void dumpGpuRecords(ostream&);


        void addArrayStats(IntArrayStats &info);
        void dumpArrayStatsToCsv(string fileName);
        
        //void addTime(int64_t t1, int64_t t2);
        //void dumpTimingsToCsv(string fileName);

        void addToAllProc(uint64_t t);
        void dumpAllProcToCsv(string fileName);

        bool IsFileEmpty(string fileName);
        void printRef(void) ; 
        

        ~RecordDatabase()
        {

        }

    protected:
        int cpuCnt;
        int gpuCnt;
        vector<Record> cpuRecords;
        vector<Record> gpuRecords;
        Record lastCpuRecord;
        Record lastGpuRecord;
        Record refRecord;
        //
        tbb::concurrent_vector<IntArrayStats> IntArrayStatsVec;
        
        //tbb::concurrent_vector<pair<uint64_t, uint64_t>> timings;

        tbb::concurrent_vector<uint64_t> allProc ;

};




#endif
