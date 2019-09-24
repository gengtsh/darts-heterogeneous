#ifndef RECORD_H_
#define RECORD_H_


#include <cfloat>
#include <unistd.h>
#include <climits>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <iterator>
#include <string>
#include <vector>
#include "tbb/concurrent_vector.h"
#include "util.h"
#include "conf.h"


using namespace std;

class Record
{
    public:
        Record()=default;
        Record(const Record&);
        Record(const string test, const int numRows, const int numNonZeroes);
        Record(const string test, const int numRows, const int numNonZeroes, uint64_t start );
        Record(const string test, const int numRows, const int numNonZeroes, uint64_t start, uint64_t end );
       
        Record& operator=(const Record&);
        
        Record(Record&&) noexcept;
        //Record& operator=(Record &&rhs) noexcept;
        
        void assignTest(string test);
        void assignStart(uint64_t st);
        void assignEnd(uint64_t ed);
        void assignValue(uint64_t val);
        void assignThreeParams(string test,int numRows, int numNonZeroes);
        
        void assignFiveParams(string tt,int nRows, int nNonZeroes, uint64_t st, uint64_t ed);

        void calcAvgNonZeroesPerRow(void);
        void calcValue(void);
        void addToValues(uint64_t val);
        void addToValues(void);
        void calcAndAddToValues(uint64_t st,uint64_t ed);
        string getTest(void) const;
        int getNumRows(void) const;
        int getNumNonZeroes(void) const;
        uint64_t getStart(void) const;
        uint64_t getEnd(void) const;
        double   getAvgNonZeroesPerRow(void) const;
        
        double calcPercentile(double q) ;
        uint64_t calcMax(void);
        uint64_t calcMin(void);
        double calcMean(void);
        double calcMedian(void);
        double calcStdDev(void); 
        
        uint64_t getMax(void)const;
        uint64_t getMin(void) const;
        double   getMean(void) const;
        
        double   getMedian(void) const;
        double   getStdDev(void) const;
        uint64_t getValue(void) const;
        //vector<uint64_t> *getValues(void) const;
        
        tbb::concurrent_vector<uint64_t> *getValues(void) const;
        void     clearRecord(void);
       
        void clearValues(void);

        void print() const;

        ~Record()
        {

        }

    public:

        string test;
        int numRows;
        int numNonZeroes;
    
    protected:

        //string test;
        //int numRows;
        //int numNonZeroes;
        uint64_t start;
        uint64_t end;
        double avgNonZeroesPerRow;
        uint64_t maxRd;
        uint64_t minRd;
        double   meanRd;
        double   medianRd;
        double   stdDevRd;
    public:
        uint64_t value;
        //vector <uint64_t> values;
        tbb::concurrent_vector <uint64_t> values;
};




#endif
