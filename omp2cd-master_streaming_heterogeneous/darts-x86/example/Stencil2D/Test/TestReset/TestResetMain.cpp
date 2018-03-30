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
#include <algorithm>
#include <iomanip>
#include <math.h>
#include <vector>
using std::vector;

#include "TestReset.h"
size_t g_nCU,g_nSU;

static inline void usage(const char *name) 
{
	std::cout << "USAGE: " << name << "<n_rows> <n_cols> <n_timesteps> <n_reps>\n";
	exit(0);
}

int main(int argc, char *argv[])
{
    size_t      nIter    = 1;
	size_t		nReps	 = 1;
    const char *str_nCU  = getenv("DARTS_NUM_CU"),
               *str_nSU  = getenv("DARTS_NUM_SU");

    g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 0;
    g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;


    switch ( argc ) {
	case 3: nReps = strtoul(argv[2],NULL,0);
	case 2: nIter = strtoul(argv[1],NULL,0);
            break;
    default: usage(argv[0]);
    }


	vector<uint64_t>timings;
    
    uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;
	
	Runtime* rt;

    outerStart = getTime();
	ThreadAffinity affin(g_nCU, g_nSU, COMPACT, TPROUNDROBIN, MCSTANDARD);

	if (affin.generateMask()){

        rt = new Runtime(&affin);
		
		for(size_t i=0; i<nReps;++i){
			innerStart = getTime();//start time for kernal procedure
			rt->run(
                launch<TestResetTP>(nIter,g_nCU, &Runtime::finalSignal)
               // launch<TestResetTP>(nIter,g_nCU+g_nSU, &Runtime::finalSignal)
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

    innerAvg /= nIter;
    outerAvg /= nIter;

    std::cout << g_nCU                  << ", "
              << g_nSU                  << ", "
              << g_nSU*(1+g_nCU)        << ", "
              << nIter					<< ", "
              << std::setprecision(18)   
              << innerAvg                << ", "
//              << outerAvg                << "\n";
              << outerAvg                << ",";

	for(auto i:timings){
		std::cout<< i <<",";
	}
	std::cout<<"\n"; 

	return 0;
}
