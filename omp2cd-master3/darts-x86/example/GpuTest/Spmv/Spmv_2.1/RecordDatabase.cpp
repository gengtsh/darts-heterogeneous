#include "RecordDatabase.h"

void RecordDatabase::assignRef(Record &rd){
    refRecord = rd;
}

void RecordDatabase::assignLastCpuRecord(Record &rd){
    lastCpuRecord = rd;
    lastCpuRecord.print();
}


void RecordDatabase::assignLastGpuRecord(Record &rd){
    lastGpuRecord = rd;
}



void RecordDatabase::cpuCntPlus1(void){
    cpuCnt ++;
}

void RecordDatabase::gpuCntPlus1(void){
    gpuCnt ++;
}

int RecordDatabase::getCpuCnt(void) const{
    return cpuCnt;
}

int RecordDatabase::getGpuCnt(void) const{
    return gpuCnt;
}

void RecordDatabase::printRef(void) {
    refRecord.print();

}

static string RemoveAllButLeadingSpaces(const string &a)
{
    string b;
    int n = a.length();
    int i = 0;
    while (i<n && a[i] == ' ')
    {
        b += a[i];
        ++i;
    }
    for (; i<n; i++)
    {
        if (a[i] != ' ' && a[i] != '\t')
            b += a[i];
    }
    return b;
}


void RecordDatabase::addRecord(string config){

    vector<Record> *results;
    Record *rd;
    if(config == "cpu"){
        results = &cpuRecords;
        rd      = &lastCpuRecord;
    }else{
        results = &gpuRecords;
        rd      = &lastGpuRecord;
    }
    string test_orig = rd->test;
    int    numRows   = rd->numRows;
    int    numNonZeroes = rd->numNonZeroes;
    string test = RemoveAllButLeadingSpaces(test_orig);

    int index;
    for (index = 0; index < results->size(); index++)
    {
        if (results->at(index).test == test &&
            results->at(index).numRows == numRows &&
            results->at(index).numNonZeroes == numNonZeroes)
        {
            results->at(index).values.push_back(rd->value); 
            break;
        }
    }

    if (index >= results->size())
    {
        results->push_back(*rd);
    }

    //rd->clearRecord();
}


Record *RecordDatabase::getLastCpuRecord(void) const{
    return &lastCpuRecord;
}

Record *RecordDatabase::getLastGpuRecord(void) const{
    return &lastGpuRecord;
}

void RecordDatabase::clearAllRecords(void){
    cpuRecords.clear();
    gpuRecords.clear();
}

void RecordDatabase::dumpCpuRecords(ostream &out){

   for (int i=0;  i<cpuRecords.size();++i){
        Record r = cpuRecords[i];
        std::cout<<r.test<<","
            <<r.numRows<<","
            <<r.numNonZeroes<<",";
        for(int j=0;j<r.values.size();++j){
            std::cout<<r.values[j]<<",";
        }
        std::cout<<std::endl;
   }
}

void RecordDatabase::dumpGpuRecords(ostream &out){

   for (int i=0;  i<gpuRecords.size();++i){
        Record r = gpuRecords[i];
        std::cout<<r.test<<","
            <<r.numRows<<","
            <<r.numNonZeroes<<",";
        for(int j=0;j<r.values.size();++j){
            std::cout<<r.values[j]<<",";
        }
        std::cout<<std::endl;
   }
}
void RecordDatabase::addArrayStats(IntArrayStats &info){
    IntArrayStatsVec.push_back(info);

}

bool RecordDatabase::IsFileEmpty(string fileName)
{
      bool fileEmpty;

      ifstream file(fileName.c_str());

      //If the file doesn't exist it is by definition empty
      if(!file.good())
      {
        return true;
      }
      else
      {
        fileEmpty = (bool)(file.peek() == ifstream::traits_type::eof());
        file.close();
        
	    return fileEmpty;
      }
  
      //Otherwise, return false  
        return false;
}


void RecordDatabase::dumpArrayStatsToCsv(string fileName){

    //bool emptyFile;
    //emptyFile = this->IsFileEmpty(fileName);
    
    //remove(fileName.c_str());

    //Open file and append by default
    ofstream out;
    out.open(fileName.c_str(), std::ofstream::out); 
    
    out<<"info  "<<","
       <<"size  "<<","
       <<"sum   "<<","
       <<"range "<<","
       <<"max   "<<","
       <<"maxIdx"<<","
       <<"min   "<<","
       <<"minIdx"<<","
       <<"mean  "<<","
       <<"md    "<<","
       <<"sd    "<<","
       <<"cv    "<<","
       <<endl;
    for (int i=0;i<IntArrayStatsVec.size();++i){
         
        out     <<IntArrayStatsVec[i].info  <<","
                <<IntArrayStatsVec[i].size  <<","
                <<IntArrayStatsVec[i].sum   <<"," 
                <<IntArrayStatsVec[i].range <<"," 
                <<IntArrayStatsVec[i].max   <<","
                <<IntArrayStatsVec[i].maxIdx<<","
                <<IntArrayStatsVec[i].min   <<","
                <<IntArrayStatsVec[i].minIdx<<","
                <<IntArrayStatsVec[i].mean  <<","
                <<IntArrayStatsVec[i].md    <<","
                <<IntArrayStatsVec[i].sd    <<","
                <<IntArrayStatsVec[i].cv    <<","
                <<endl;
    
    
    }



    out << endl;
    out.close();

}

//void RecordDatabase::addTime(int64_t t1,int64_t t2){
//    timings.push_back(std::make_pair(t1,t2)); 
//}
//
//void RecordDatabase::dumpTimingsToCsv(string fileName){
//   
//    ofstream out;
//    out.open(fileName.c_str(), std::ofstream::out); 
//    
//    out<<"exe time  "<<","
//       <<"get stats time "<<","
//       <<"diff"<<","
//       <<endl;
//    
//    for(int i=0;i<timings.size();++i){
//        int64_t diff = timings[i].first - timings[i].second;
//        out<<timings[i].first<<","
//            <<timings[i].second<<","
//            <<diff<<","
//            <<endl;
//    
//    }
//    out << endl;
//    out.close();
//}



void RecordDatabase::addToAllProc(uint64_t t1){
    allProc.push_back(t1); 
}

void RecordDatabase::dumpAllProcToCsv(string fileName){
   
    ofstream out;
    out.open(fileName.c_str(), std::ofstream::out); 
    
    out<<"exe time  "<<","
       <<endl;
   
    uint64_t sum = 0;
    for(int i=0;i<allProc.size();++i){
        out<<allProc[i]<<","
            <<endl;
        sum +=allProc[i];
    
    }
    out << endl;
    out<<"total exe time: "<< sum <<endl; 
    out.close();
}

