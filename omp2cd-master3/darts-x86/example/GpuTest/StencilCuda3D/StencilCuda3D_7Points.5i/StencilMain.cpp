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
	std::cout << "USAGE: " << name << "<n_rows> <n_cols> <n_slices> <n_timesteps> <n_reps>\n";
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
    const double  epsilon    = 0.05;

    for(size_t i = 0; is_correct && i < n_rows; ++i) 
		for(size_t j = 0; is_correct && j < n_cols; ++j) 
			if ( !( is_correct = fabs( RES[i][j] - ORIG[i][j] ) <= epsilon ) )
				printf("Values mismatch! [%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",i, j, ORIG[i][j], RES[i][j]);

	return is_correct;
}



bool
result_is_correct(const size_t  n_rows, const size_t  n_cols, const size_t n_slices, const double* orig,const double* res)
{
    typedef double (*Array3D)[n_rows][n_cols];
    const Array3D ORIG       = (const Array3D) orig,
                RES        = (const Array3D) res;
    bool          is_correct =   true;
    const double  epsilon    = 0.0001;
    
    for(size_t k = 0; is_correct && k < n_slices; ++k) 
        for(size_t i = 0; is_correct && i < n_rows; ++i) 
	    	for(size_t j = 0; is_correct && j < n_cols; ++j){
	    		if ( !( is_correct = fabs( RES[k][i][j] - ORIG[k][i][j] ) <= epsilon ) )
	    			printf("Values mismatch! [%lu,%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",k,i, j,ORIG[k][i][j], RES[k][i][j]);
            }
	return is_correct;
}

void print_results(const double *results, const size_t  n_rows_st,const size_t n_rows_ed, const size_t  n_cols_st, const size_t n_cols_ed,const size_t n_slices_st,const size_t n_slices_ed,const size_t n_rows, const size_t n_cols, const size_t n_slices ){
	for(size_t s=n_slices_st;s<n_slices_ed;++s){
        std::cout<<"slice: "<<s<<std::endl;
        for(size_t k=n_rows_st;k<n_rows_ed;++k){
		    for(size_t j=n_cols_st;j<n_cols_ed;++j){
			    //std::cout<<std::setw(18)<< results[s*n_rows*n_cols+k*n_cols+j]<<",";
			    std::cout<< results[s*n_rows*n_cols+k*n_cols+j]<<",";
		    }
		    std::cout<<"\n";
	    }
		std::cout<<"\n";
    }

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


void init3Darray ( double* a, const size_t n_rows, const size_t n_cols, const size_t n_slices )
{
    typedef double (*Array3D)[n_rows][n_cols];
    Array3D A = (Array3D) a;
    for (size_t k = 0; k< n_slices; ++k){
        for (size_t i = 0; i < n_rows; ++i) 
    	{
    		for (size_t j = 0; j < n_cols; ++j) 
    		{
    			A[k][i][j] = (1+i)*j + 1.3+k+1;
        //        std::cout<<A[k][i][j]<<",";
    		}
    	//	std::cout<<"\n";
    	}
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


void
StencilSeq ( double* __restrict__ dst,    double* __restrict__ src, 
               const size_t     n_rows, const size_t     n_cols, const size_t n_slices,
               const size_t     n_tsteps )
{
    typedef double (*Array3D)[n_rows][n_cols];
    volatile Array3D DST = (Array3D) dst,
                     SRC = (Array3D) src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
        for (size_t k = 1; k< n_slices-1;++k){
	    	for (size_t i = 1; i < n_rows-1; ++i) {
                for (size_t j = 1; j < n_cols-1; ++j) {
                    DST[k][i][j] = (SRC[k][i-1][j] + SRC[k][i+1][j] + SRC[k][i][j-1] + SRC[k][i][j+1] + SRC[k][i][j]+SRC[k-1][i][j]+SRC[k+1][i][j])/7.5;
                   // if(j<5 && i<5 && k<3){
                   //     std::cout<<DST[k][i][j]<<",";
                   // }
                }

                //if(i<5 && k<3){
                //    std::cout<<"\n";
                //}
            }
            //if(k<3){
                //std::cout<<"\n"<<std::endl;
            //}
        }
	    SWAP_PTR(&DST,&SRC);
        
       //std::cout<<DST[1][1][1]<<std::endl;
       //std::cout<<SRC[1][1][1]<<std::endl;
    }
}


int main(int argc, char *argv[])
{
    size_t      nRows    = 0, 
                nCols    = 0,
                nSlices = 0,
                nTmSteps = 0,
                nReps    = 0;
    const char *str_nCU  = getenv("DARTS_NUM_CU"),
               *str_nSU  = getenv("DARTS_NUM_SU");

    g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 0;
    g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;

    nTmSteps = nReps = 1;
    nSlices = 1;
    size_t IsStream = 0;
    size_t nTiles = 0;

    switch ( argc ) {
    
    case 8: nTiles  = strtoul(argv[7],NULL,0);
    case 7: IsStream  = strtoul(argv[6],NULL,0);
    case 6: nReps    = strtoul(argv[5],NULL,0);
    case 5: nTmSteps = strtoul(argv[4],NULL,0);
    case 4: nSlices  = strtoul(argv[3],NULL,0);
    case 3: nRows    = strtoul(argv[1],NULL,0); 
            nCols    = strtoul(argv[2],NULL,0);
            break;
    case 2: nRows = nCols = strtoul(argv[1],NULL,0);
            break;
    default: usage(argv[0]);
    }
   
#ifdef CUDA_DARTS_DEBUG
    std::cout<<"nRows:"<<nRows<<",nCols:"<<nCols<<",nSlices:"<<nSlices<<",timestep:"<<nTmSteps<<",nReps:"<<nReps<<std::endl;
#endif

    bool streamming = (IsStream==0)?false:true;    
    double GpuRatio = 0.6; // 0: pure CPU, 1: pure GPU, (0,1) hybrid
    //bool    streamming = true;	
    //bool    streamming = false;	
    double* OriginalMatrix = new double[nRows*nCols*nSlices];

	double* InitialMatrix ; 	
    double* NewMatrix  ;
    if(GpuRatio == 0.0){
	    InitialMatrix  = new double [ nRows*nCols*nSlices ];
	    NewMatrix  = new double [ nRows*nCols*nSlices ];
    }else{
	    cudaMallocHost((void**)&InitialMatrix,sizeof(double)*nRows*nCols*nSlices) ;
	    cudaMallocHost((void**)&NewMatrix,sizeof(double)*nRows*nCols*nSlices) ;
    }
  
    init3Darray(OriginalMatrix,nRows,nCols,nSlices);// initial 2D array

    uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;

	vector<uint64_t>timings;
    Runtime* rt;
	bool hard = false;
    bool IsStatic = true;
    outerStart = getTime();
	//ThreadAffinity affin(g_nCU, g_nSU, COMPACT, TPROUNDROBIN, MCSTANDARD);
	ThreadAffinity affin(g_nCU, g_nSU, COMPACT, 6, 4);
    //ThreadAffinity affin(g_nCU, g_nSU, SPREAD, 6, 4);

    if (affin.generateMask()){

        rt = new Runtime(&affin);
		
		for(size_t i=0; i<nReps;++i){

            std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols*nSlices, InitialMatrix);
            std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols*nSlices, NewMatrix);
			innerStart = getTime();//start time for kernal procedure
			rt->run(
                launch<StencilTP>(InitialMatrix,nRows,nCols,nSlices,NewMatrix,nTmSteps,hard, IsStatic, GpuRatio,streamming,nTiles,&Runtime::finalSignal)
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


    
    std::cout << nRows << " * " << nCols << " * "
              << nSlices                 << ", "
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
	
#ifdef CUDA_DARTS_DEBUG
    std::cout<<"begin verify:"<<std::endl;
#endif
	double* SeqMatrix = new double[nRows*nCols*nSlices];
    double* SeqOutMatrix = new double[nRows*nCols*nSlices];
	std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols*nSlices, SeqMatrix);
	std::copy(OriginalMatrix, OriginalMatrix+nRows*nCols*nSlices, SeqOutMatrix);
	
	//StencilSeq ( SeqOutMatrix,SeqMatrix,nRows,nCols,nTmSteps); 
	StencilSeq ( SeqOutMatrix,SeqMatrix,nRows,nCols,nSlices,nTmSteps); 

#ifdef CUDA_DARTS_DEBUG
//    std::cout<<"SeqOutMatrix[1][1][1]:"<<SeqOutMatrix[1*nRows*nCols+nCols+1]<<std::endl;
//    std::cout<<"SeqMatrix[1][1][1]:"<<SeqMatrix[1*nRows*nCols+nCols+1]<<std::endl;
#endif
    
    double *seqOut, *dartsOut;
	
	if(nTmSteps%2){
		seqOut=SeqOutMatrix;
        dartsOut=NewMatrix;
    }else{
		seqOut = SeqMatrix;
        dartsOut = InitialMatrix;
	}


	//if(result_is_correct(nRows,nCols,seqOut,dartsOut)){
	if(result_is_correct(nRows,nCols,nSlices,seqOut,dartsOut)){
		std::cout<<"success"<<std::endl;
	}else{
		std::cout<<"fail!"<<std::endl;
	}
#endif
	
#ifdef VERIFICATION_PRINT
	std::cout<<std::setprecision(10)<<std::endl;
	int kk = 0;
	int ttk =5;
	int jj = 0; 
	int ttj =5;
    int ss =0;
    int tts=3;
	std::cout<<"Original Matrix:"<<std::endl;
	print_results(OriginalMatrix,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
    std::cout<<"Seq Matrix:"<<std::endl;
	print_results(seqOut,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
	std::cout<<"Out Matrix:"<<std::endl;
	print_results(dartsOut,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
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
