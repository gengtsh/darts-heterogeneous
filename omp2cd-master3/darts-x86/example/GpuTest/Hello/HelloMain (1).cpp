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
#include <math.h>

#include <vector>
using std::vector;

#include "HelloTP.h"

size_t g_nCU,g_nSU;
static inline void usage(const char *name) 
{
	std::cout << "USAGE: " << name << "<n_timesteps> <n_reps>\n";
	exit(0);
}

void init2Darray ( double* a, const size_t n_rows, const size_t n_cols )
{
    typedef double (*Array2D)[n_cols];
    Array2D A = (Array2D) a;
    for (size_t i = 0; i < n_rows; ++i) 
        for (size_t j = 0; j < n_cols; ++j) 
            A[i][j] = (1+i)*j + 1.3;
}

int main(int argc, char *argv[])
{

//    size_t      nRows    = 0, 
//                nCols    = 0;
//	size_t      nTmSteps = 0,
//                nReps    = 0;
    const char *str_nCU  = getenv("DARTS_NUM_CU"),
               *str_nSU  = getenv("DARTS_NUM_SU");

	g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 0;
	g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;
//
//    nTmSteps=1;
//	nReps = 1;
//
//    switch ( argc ) {
//    case 3: nReps    = strtoul(argv[2],NULL,0);
//    case 2: nTmSteps = strtoul(argv[1],NULL,0);
//			break;
//    default: usage(argv[0]);
//    }
	
//	nRows=1000;
//	nCols=1000;

//    switch ( argc ) {
//    case 5: nReps    = strtoul(argv[4],NULL,0);
//    case 4: nTmSteps = strtoul(argv[3],NULL,0); 
//    case 3: nRows    = strtoul(argv[1],NULL,0); 
//            nCols    = strtoul(argv[2],NULL,0);
//            break;
//    case 2: nRows = nCols = strtoul(argv[1],NULL,0);
//            break;
//    default: usage(argv[0]);
//    }

//    switch ( argc ) {
//    case 3: nRows    = strtoul(argv[1],NULL,0); 
//            nCols    = strtoul(argv[2],NULL,0);
//            break;
//    case 2: nRows = nCols = strtoul(argv[1],NULL,0);
//            break;
//    default: usage(argv[0]);
//    }

	//nReps = 1;
	//nTmSteps = 5;


    uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;

//	vector<uint64_t>timings;
    Runtime* rt;

    outerStart = getTime();
	ThreadAffinity affin(g_nCU, g_nSU, COMPACT, 6,4 );

	if (affin.generateMask()){

        rt = new Runtime(&affin);
		
			innerStart = getTime();//start time for kernal procedure
			//std::cout<<"rt->run!\n"<<std::endl;
			rt->run(
              //launch<TestS1TP>(InitialMatrix,nRows,nCols,NewMatrix,nTmSteps, &Runtime::finalSignal)
              launch<HelloTP>( &Runtime::finalSignal)
					);
			innerStop = getTime() - innerStart;
//			innerAvg += innerStop;
//			timings.push_back(innerStop);
		
	}else{
		std::cerr << "Could not generate required abstract machine -- something went wrong. :(\n";
		return EXIT_FAILURE;
	}
    delete rt;
    outerStop = getTime() - outerStart;
//    outerAvg += outerStop;

//    innerAvg /= nReps;
//    outerAvg /= nReps;

    std::cout << g_nCU                   << ", "
              << g_nSU                   << ", "
              << g_nSU*(1+g_nCU)         << ", "
              << std::setprecision(18)   
              << innerStart                << ", "
              << innerStop                << ",";

//	for(auto i:timings){
//		std::cout<< i <<",";
//	}
	std::cout<<"\n"; 
	
	return 0;

}
