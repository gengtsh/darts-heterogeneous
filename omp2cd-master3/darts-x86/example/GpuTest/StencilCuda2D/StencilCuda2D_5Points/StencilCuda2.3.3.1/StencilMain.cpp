/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * DARTS - A fine-grain dataflow-inspired runtime system.                          *
 * Copyright (C) 2011-2014  University of Delaware                                 *
 *                                                                                 *
 * This library is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU Lesser General Public                      *
 * License as published by the Free Software Foundation; either                    *
 * version 2.1 of the License, or (at your option) any later version.              *
 *                                                                                 *
 * This library is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU               *
 * Lesser General Public License for more details.                                 *
 *                                                                                 *
 * You should have received a copy of the GNU Lesser General Public                *
 * License along with this library; if not, write to the Free Software             *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA  *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <unistd.h>
#include <climits>
#include <cstdio>
#include <cerrno>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "StencilMain.h"
#include <math.h>
#include "stencil.h"
#include <vector>
using std::vector;

size_t g_nCU,g_nSU;
static inline void usage(const char *name) 
{
	std::cout << "USAGE: " << name << "<n_rows> <n_cols> <n_timesteps> <n_reps>\n";
	exit(0);
}

/**
* 
* inputs:
* n_rows -- the number of rows in <orig> and <res> 
* n_cols -- the number of columns in <orig> and <res>
* orig   -- a "true" 2D array which holds the reference results
* res    -- a "true" 2D array which holds the values to be verified
*
* output:
* True if <orig> and <res> have identical results within an epsilon. False
* otherwise.
*/
bool
result_is_correct(const size_t  n_rows, const size_t  n_cols,const double* orig,   const double* res)
{
    typedef double (*Array2D)[n_cols];
    const Array2D ORIG       = (const Array2D) orig,
                RES        = (const Array2D) res;
    bool          is_correct =   true;
    const double  epsilon    = 0.0001;

    for(size_t i = 0; is_correct && i < n_rows; ++i) 
		for(size_t j = 0; is_correct && j < n_cols; ++j) 
			if ( !( is_correct = fabs( RES[i][j] - ORIG[i][j] ) <= epsilon ) )
				printf("Values mismatch! [%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",i, j, ORIG[i][j], RES[i][j]);

	return is_correct;
}

/**
 * brief <init2Darray> takes a 2D array (assumed to be allocated in contiguous
 * memory) and initializes it with arbitrary (well-formed) double precision
 * values.
 *
 * inputs:
 * a      -- the 2D array to initialize
 * n_rows -- the number of rows in <a>
 * n_cols -- the number of columns in <a>
 *
 * output: N/A
 */
void init2Darray ( double* a, const size_t n_rows, const size_t n_cols )
{
    typedef double (*Array2D)[n_cols];
    Array2D A = (Array2D) a;
    for (size_t i = 0; i < n_rows; ++i) 
	{
		for (size_t j = 0; j < n_cols; ++j) 
		{
			A[i][j] = (1+i)*j + 1.3;
		//	std::cout<<A[i][j]<<",";
		}
		//std::cout<<"\n";
	}
}


void
StencilSeq ( double* __restrict__ dst,    double* __restrict__ src, 
               const size_t     n_rows, const size_t     n_cols,
               const size_t     n_tsteps )
{
    typedef double (*Array2D)[n_cols];
    volatile Array2D DST = (Array2D) dst,
            SRC = (Array2D) src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
		for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1])/5.5;
                //DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1])/4;
            }
        }
		SWAP_PTR(&DST,&SRC);
    }
}


int main(int argc, char *argv[])
{
    size_t      nRows    = 0, 
                nCols    = 0,
                nTmSteps = 0,
                nReps    = 0;
    const char *str_nCU  = getenv("DARTS_NUM_CU"),
               *str_nSU  = getenv("DARTS_NUM_SU");

    g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 0;
    g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;

    nTmSteps = nReps = 1;

//	size_t nCmpRep = 1;
    switch ( argc ) {
    
    //case 6: nCmpRep    = strtoul(argv[5],NULL,0);
    case 5: nReps    = strtoul(argv[4],NULL,0);
    case 4: nTmSteps = strtoul(argv[3],NULL,0); 
    case 3: nRows    = strtoul(argv[1],NULL,0); 
            nCols    = strtoul(argv[2],NULL,0);
            break;
    case 2: nRows = nCols = strtoul(argv[1],NULL,0);
            break;
    default: usage(argv[0]);
    }
    
    double GpuRatio = 0.6; // 0: pure CPU, 1: pure GPU, (0,1) hybrid

	double* OriginalMatrix = new double[nRows*nCols];

	double* InitialMatrix ; 	
    double* NewMatrix  ;
    if(GpuRatio == 0.0){
	    InitialMatrix  = new double [ nRows*nCols ];
	    NewMatrix  = new double [ nRows*nCols ];
    }else{
	    cudaMallocHost((void**)&InitialMatrix,sizeof(double)*nRows*nCols) ;
	    cudaMallocHost((void**)&NewMatrix,sizeof(double)*nRows*nCols) ;
    }
  
    init2Darray(OriginalMatrix,nRows,nCols);// initial 2D array

    uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;

	vector<uint64_t>timings;
    Runtime* rt;
	bool hard = false;
    outerStart = getTime();
	//ThreadAffinity affin(g_nCU, g_nSU, COMPACT, TPROUNDROBIN, MCSTANDARD);
	
	//ThreadAffinity affin(g_nCU, g_nSU, COMPACT, 6, 4);
    ThreadAffinity affin(g_nCU, g_nSU, SPREAD, 6, 4);
	if (affin.generateMask()){

        rt = new Runtime(&affin);
		
		for(size_t i=0; i<nReps;++i){

            std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols, InitialMatrix);
            std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols, NewMatrix);
			innerStart = getTime();//start time for kernal procedure
			rt->run(
                launch<StencilTP>(InitialMatrix,nRows,nCols,NewMatrix,nTmSteps,hard, GpuRatio,&Runtime::finalSignal)
            );
			innerStop = getTime() - innerStart;
			innerAvg += innerStop;
			timings.push_back(innerStop);
		}
	}else{
		std::cerr << "Could not generate required abstract machine -- something went wrong. :(\n";
		return EXIT_FAILURE;
	}
    delete rt;
    outerStop = getTime() - outerStart;
    outerAvg += outerStop;

    innerAvg /= nReps;
    outerAvg /= nReps;
	
//	std::sort(timings.begin(),timings.end());
//	innerAvg = 0.0;
//	for (auto j=timings.begin()+1;j<timings.end()-1;++j){
//		innerAvg +=*j;
//	
//	}
//	innerAvg /=(nReps-2);


    std::cout << nRows << " * " << nCols << ", " 
              << g_nCU                   << ", "
              << g_nSU                   << ", "
              << g_nSU*(1+g_nCU)         << ", "
              << nTmSteps                << ", "
              << nReps                   << ", "
              << std::setprecision(18)   
              << innerAvg                << ", "
//              << outerAvg                << "\n";
              << outerAvg                << ",";

	for(auto i:timings){
		std::cout<< i <<",";
	}
	std::cout<<"\n"; 

	//for verification:

#ifdef VERIFICATION
	double* SeqMatrix = new double[nRows*nCols];
    double* SeqOutMatrix = new double[nRows*nCols];
	std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols, SeqMatrix);
	std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols, SeqOutMatrix);
	
	StencilSeq ( SeqOutMatrix,SeqMatrix,nRows,nCols,nTmSteps); 

	double *seqOut, *dartsOut;
	
	if(nTmSteps%2){

		seqOut=SeqOutMatrix;
		dartsOut=NewMatrix;
	}else{

		seqOut = SeqMatrix;
		dartsOut = InitialMatrix;
	}

	if(result_is_correct(nRows,nCols,seqOut,dartsOut)){
		std::cout<<"success"<<std::endl;
	}else{
		std::cout<<"fail!"<<std::endl;
	}
#endif
	
#ifdef VERIFICATION_PRINT 
	std::cout<<std::setprecision(6)<<std::endl;
	int kk = 0;
	int ttk = 10;
	int jj = 0;
	int ttj = 10;
	std::cout<<"Original Matrix:"<<std::endl;
	for(size_t k=kk;k<kk+ttk;++k){
		for(size_t j=jj;j<jj+ttj;++j){
			std::cout<< OriginalMatrix[k*nCols+j]<<",";
		}
		std::cout<<"\n";
	}
	std::cout<<"Seq Matrix:"<<std::endl;
	for(size_t k=kk;k<kk+ttk;++k){
		for(size_t j=jj;j<jj+ttj;++j){
			std::cout<< seqOut[k*nCols+j]<<",";
		}
		std::cout<<"\n";
	}

	std::cout<<"Out Matrix:"<<std::endl;
	for(size_t k=kk;k<kk+ttk;++k){
		for(size_t j=jj;j<jj+ttj;++j){
			std::cout<< dartsOut[k*nCols+j]<<",";
		}
		std::cout<<"\n";
	}
#endif

#ifdef VERIFICATION
	delete[] SeqMatrix;
	delete[] SeqOutMatrix;
#endif
    delete[] OriginalMatrix;

    
    if(GpuRatio == 0){
    	delete[] InitialMatrix;
        delete[] NewMatrix;
    }else{
        cudaFreeHost(InitialMatrix);
        cudaFreeHost(NewMatrix);
    }


	return 0;
}
