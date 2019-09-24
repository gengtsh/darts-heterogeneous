#include "Record.h"


#include <cassert>
#include <pthread.h>
pthread_mutex_t vecMutex;
#include <iostream>

//	pthread_mutex_lock(&mutex);
//	pthread_mutex_unlock(&mutex);


using namespace std;

Record::Record(const string test, const int numRows, const int numNonZeroes )
    :test(test)
    ,numRows(numRows)
    ,numNonZeroes(numNonZeroes)
    {
        start = 0;
        end   = 0;
        calcAvgNonZeroesPerRow();
        maxRd = 0;
        minRd = 0;
        meanRd = 0;
        medianRd = 0;
        stdDevRd = 0;
        value    = 0;
    }


Record::Record(const string test, const int numRows, const int numNonZeroes, uint64_t start)
    :test(test)
    ,numRows(numRows)
    ,numNonZeroes(numNonZeroes)
    ,start(start)
    {
        end   = 0;
        calcAvgNonZeroesPerRow();

        maxRd = 0;
        minRd = 0;
        meanRd = 0;
        medianRd = 0;
        stdDevRd = 0;
        value   = 0; 
    }

Record::Record(const string test, const int numRows, const int numNonZeroes, uint64_t start, uint64_t end )
    :test(test)
    ,numRows(numRows)
    ,numNonZeroes(numNonZeroes)
    ,start(start)
    ,end(end)
    {
        calcAvgNonZeroesPerRow();
        calcValue();
        calcMean();

        maxRd = 0;
        minRd = 0;
        medianRd = 0;
        stdDevRd = 0;
    }



Record::Record(const Record& rhs)
        :test              (rhs.test              )    
        ,numRows           (rhs.numRows           )  
        ,numNonZeroes      (rhs.numNonZeroes      )  
        ,start             (rhs.start             )  
        ,end               (rhs.end               )  
        ,avgNonZeroesPerRow(rhs.avgNonZeroesPerRow)  
        ,maxRd             (rhs.maxRd             )    
        ,minRd             (rhs.minRd             )    
        ,meanRd            (rhs.meanRd            )    
        ,medianRd          (rhs.medianRd          )    
        ,stdDevRd          (rhs.stdDevRd          )    
        ,value             (rhs.value             )   
        ,values            (rhs.values            )   
        {

#ifdef DARTS_DEBUG
            std::cout<<"record copy copy copy"<<std::endl;
#endif
        }
Record::Record(Record&& rhs) noexcept
        :test              (rhs.test              )    
        ,numRows           (rhs.numRows           )  
        ,numNonZeroes      (rhs.numNonZeroes      )  
        ,start             (rhs.start             )  
        ,end               (rhs.end               )  
        ,avgNonZeroesPerRow(rhs.avgNonZeroesPerRow)  
        ,maxRd             (rhs.maxRd             )    
        ,minRd             (rhs.minRd             )    
        ,meanRd            (rhs.meanRd            )    
        ,medianRd          (rhs.medianRd          )    
        ,stdDevRd          (rhs.stdDevRd          )    
        ,value             (rhs.value             )   
        ,values            (rhs.values            )   
        {
#ifdef DARTS_DEBUG
        std::cout<<"record move move move"<<std::endl;
#endif
        }

Record& Record::operator=(const Record& rhs){
    test               = rhs.test               ;
    numRows            = rhs.numRows            ;
    numNonZeroes       = rhs.numNonZeroes       ;
    start              = rhs.start              ;
    end                = rhs.end                ;
    avgNonZeroesPerRow = rhs.avgNonZeroesPerRow ;
    maxRd              = rhs.maxRd              ;    
    minRd              = rhs.minRd              ;    
    meanRd             = rhs.meanRd             ;    
    medianRd           = rhs.medianRd           ;    
    stdDevRd           = rhs.stdDevRd           ;    
    value              = rhs.value             ;
    values             = rhs.values             ;
    
#ifdef DARTS_DEBUG
    std::cout<<"record operator ===="<<std::endl;
    std::cout<<"rhs.numRows = " << rhs.numRows<<std::endl;
#endif
    return *this;
}


//Record& Record::operator=(Record &&rhs) noexcept{
//    if( this !=&rhs){
//        test               = rhs.test               ;
//        numRows            = rhs.numRows            ;
//        numNonZeroes       = rhs.numNonZeroes       ;
//        start              = rhs.start              ;
//        end                = rhs.end                ;
//        avgNonZeroesPerRow = rhs.avgNonZeroesPerRow ;
//        maxRd              = rhs.maxRd              ;    
//        minRd              = rhs.minRd              ;    
//        meanRd             = rhs.meanRd             ;    
//        medianRd           = rhs.medianRd           ;    
//        stdDevRd           = rhs.stdDevRd           ;    
//        value              = rhs.value             ;
//        values             = rhs.values             ;
//        
//    }
//        return *this;
//}

void Record::assignTest(string tt){
    test = tt;
}
void Record::assignStart(uint64_t st){
    start = st;
}

void Record::assignEnd(uint64_t ed){
    end = ed;
}

void Record::assignValue(uint64_t val){
    value = val;
}

void Record::assignThreeParams(string tt,int nRows, int nNonZeroes){
    test = tt;
    numRows = nRows;
    nNonZeroes = nNonZeroes;
    calcAvgNonZeroesPerRow();
}


void Record::assignFiveParams(string tt,int nRows, int nNonZeroes, uint64_t st, uint64_t ed){
    test = tt;
    numRows = nRows;
    nNonZeroes = nNonZeroes;
    start = st;
    end   = ed;
    calcAvgNonZeroesPerRow();
    calcAndAddToValues(start,end);
}

void Record::calcAvgNonZeroesPerRow(){
    avgNonZeroesPerRow = (1.0*numNonZeroes)/numRows;
}

void Record::calcValue(void){
    value = end - start;
}

void Record::addToValues(uint64_t val){
    
	//pthread_mutex_lock(&vecMutex);
    values.push_back(val);
	//pthread_mutex_unlock(&vecMutex);
}



void Record::addToValues(void){
    
	//pthread_mutex_lock(&vecMutex);
    values.push_back(value);
	//pthread_mutex_unlock(&vecMutex);
}

void Record::calcAndAddToValues(uint64_t st,uint64_t ed){
    
    value = ed - st;
    addToValues();
    start = (start ==0)? st:( (start < st)? start: st);   
    end   = (end > ed)? end : ed;
    //calcMean();
}

double Record::calcPercentile(double q) 
{
    int n = values.size();
    if (n == 0)
        return FLT_MAX;
    if (n == 1)
        return values[0];

    if (q <= 0)
        return values[0];
    if (q >= 100)
        return values[n-1];

    double index = ((n + 1.) * q / 100.) - 1;

    //vector<uint64_t> sorted = static_cast<vector<uint64_t>> (values);
    tbb::concurrent_vector<uint64_t> sorted = values;
    sort(sorted.begin(), sorted.end());

    if (n == 2)
        return (sorted[0] * (1 - q/100.)  +  sorted[1] * (q/100.));

    int index_lo = int(index);
    double frac = index - index_lo;
    if (frac == 0)
        return sorted[index_lo];

    double lo = sorted[index_lo];
    double hi = sorted[index_lo + 1];
    return lo + (hi-lo)*frac;
}


string Record::getTest() const{
    return test;
}

int Record::getNumRows(void)const {
    return numRows;
}

int Record::getNumNonZeroes(void) const{
    return numNonZeroes;
}

double Record::getAvgNonZeroesPerRow(void) const{
    return avgNonZeroesPerRow;
}


uint64_t Record::getStart(void) const{
    return start;
}
uint64_t Record::getEnd(void) const{
    return end;
}

uint64_t Record::getMax(void) const{
    return maxRd;
}
uint64_t Record::getMin(void) const{
    return minRd;
}
double Record::getMean(void) const {
    return meanRd;
}

double Record::getMedian(void) const {
    return medianRd;
}

double Record::getStdDev(void) const {
    return stdDevRd;
}

uint64_t Record::getValue(void) const {
    return value;
}

//vector<uint64_t> * Record::getValues(void) const {
tbb::concurrent_vector<uint64_t> * Record::getValues(void) const {
    return &values;
}

uint64_t Record::calcMax(void) 
{
    uint64_t r = -FLT_MAX;
    for (int i=0; i<values.size(); i++)
    {
        r = max(r, values[i]);
    }
    maxRd = r;
    return r;
}

uint64_t Record::calcMin(void) 
{
    uint64_t r = FLT_MAX;
    for (int i=0; i<values.size(); i++)
    {
        r = min(r, values[i]);
    }
    minRd = r;
    return r;
}

double Record::calcMean(void) 
{
    double r = 0;
    for (int i=0; i<values.size(); i++)
    {
        r += values[i];
    }
    meanRd =  r / double(values.size());
    return meanRd;
}

double Record::calcMedian(void){
    medianRd = calcPercentile(50);
    return medianRd;
}

double Record::calcStdDev(void) 
{
    double r = 0;
    double u;
    if (meanRd == 0){
        calcMean();
    }
    u = meanRd;
    if (u == FLT_MAX)
        return FLT_MAX;
    for (int i=0; i<values.size(); i++)
    {
        r += (values[i] - u) * (values[i] - u);
    }
    r = sqrt(r / values.size());
    stdDevRd = r;
    return r;
}


void Record::print(void) const{

    std::cout
        << test << ","
        << numRows << ","
        << numNonZeroes << ","
        << avgNonZeroesPerRow << ","
        << meanRd <<","
        ;
    //    <<std::endl;
        for (int i=0; i<values.size();++i){
            std::cout<<values[i]<<",";
        }
        std::cout<<std::endl;
}

void Record::clearValues(void){
    values.clear();
}

void Record::clearRecord(void){

        start = 0;
        end   = 0;
        avgNonZeroesPerRow =0;
        maxRd = 0;
        minRd = 0;
        meanRd = 0;
        medianRd = 0;
        stdDevRd = 0;
        value    = 0;
        clearValues();
}
