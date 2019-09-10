#include <unistd.h>
#include <climits>
#include <cstdio>
#include <cctype>
#include <cerrno>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <string>
#include <strings.h>
#include "OptionParser.h"
#include "ResultDatabase.h"
#include "RecordDatabase.h"
#include "util.h"
//#include "CSRMM.h"
#include "DARTS.h"
#include "SpmvTP.h"
#include "conf.h"

using namespace std;
size_t g_nCU,g_nSU;


static inline void usage(const char *name) 
{
	std::cout << "USAGE: " << name << " input1:  random or file (matrix market format) name!" <<"  input2: cpu or gpu or hybrid"<<std::endl;
	exit(0);
}

void removeSubstrs(string& s, const string& p){
    string::size_type findpos = s.find(p);
    if(findpos != string::npos){
        s.erase(s.begin()+findpos, s.begin()+findpos+p.length());
    }
}

int main(int argc, char *argv[])
{

    const char *str_nCU  = getenv("DARTS_NUM_CU");
    const char *str_nSU  = getenv("DARTS_NUM_SU");

    g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 0;
    g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;

    string inMatrix;
    string config = "cpu";
   // double GpuRatio = 0; // 0: pure CPU, 1: pure GPU, (0,1) hybrid
    
    switch ( argc ) {
        case 3: config = argv[2];
        case 2: inMatrix = argv[1];break;
        default: usage(argv[0]);
    }
    std::for_each(config.begin(),config.end(),[](char &c){c =tolower(c); });
    // Get args
    OptionParser op;

    op.addOption("passes", OPT_INT, "1", "specify number of passes", 'n');
    op.addOption("iterations", OPT_INT, "1", "Number of SpMV iterations per pass");
    op.addOption("mm_filename", OPT_STRING, inMatrix, inMatrix +" which stores the matrix in Matrix Market format");
    op.addOption("maxval", OPT_FLOAT, "10", "Maximum value for random matrices");
    op.addOption("config", OPT_STRING, config, "code will be run in "+config+" way");
    //op.print(); 

    //addBenchmarkSpecOptions(op);

    ////CSRMM<double> csrHost("host",config);


    ////// This benchmark either reads in a matrix market input file or generates a random matrix
    ////string inFileName = op.getOptionString("mm_filename");
    ////if (inFileName == "random"){

    ////}else{
    ////    char filename[FIELD_LENGTH];
    ////    strcpy(filename, inFileName.c_str());
    ////    readMarketMatrixToCSR(filename, &csrHost); 
    ////}

    ResultDatabase resultDB;
    RecordDatabase recordDB;
    
    ThreadAffinity affin(g_nCU, g_nSU, COMPACT, 6, 4);
    //ThreadAffinity affin(g_nCU, g_nSU, SPREAD, 6, 4);

    Runtime* rt;
    if (affin.generateMask()){
        rt = new Runtime(&affin);
        
        rt->run (
                launch<SpmvTP>(&recordDB,&resultDB,&op,&Runtime::finalSignal)
                );
    }else{
		std::cerr << "Could not generate required abstract machine -- something went wrong. :(\n";
		return EXIT_FAILURE;
    }
    delete rt;
    
    // Run the benchmark
    //RunBenchmark(resultDB, op);

    //resultDB.DumpDetailedMix(cout);
    //resultDB.DumpCsvMix("1.csv");
    //recordDB.dumpCpuRecords(cout);

    //recordDB.dumpArrayStatsToCsv("rowstats.csv");
    //recordDB.dumpTimingsToCsv("timings.csv");
    
    removeSubstrs(inMatrix,".mtx");
    string outfile=inMatrix+"_"+config+"_allProc.csv";
    recordDB.dumpAllProcToCsv(outfile);
}
