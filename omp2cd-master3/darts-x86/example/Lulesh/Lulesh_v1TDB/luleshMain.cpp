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
#include <cstdio>
#include <cerrno>
#include <algorithm>
//#include <mkl_cblas.h>

#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>

#include "lulesh.h"
#include "luleshMain.h"
#include "LagrangeLeapFrogTP.h"

//#if _OPENMP
//# include <omp.h>
//#endif
size_t g_treeBarrier;
size_t g_nCU,g_nSU;

/******************************************/

int main(int argc, char *argv[])
{
  Domain *locDom ;
   Int_t numRanks ;
   Int_t myRank ;
   struct cmdLineOpts opts;

   numRanks = 1;
   myRank = 0;

   /* Set defaults that can be overridden by command line opts */
   //opts.its = 9999999;
   opts.its = 30;
   opts.nx  = 30;
   opts.numReg = 11;
   opts.numFiles = (int)(numRanks+10)/9;
   opts.showProg = 0;
   opts.quiet = 0;
   opts.viz = 0;
   opts.balance = 1;
   opts.cost = 1;
   opts.treeBarrier= 2;
   
   ParseCommandLineOptions(argc, argv, myRank, &opts);

	const char *str_nCU  = getenv("DARTS_NUM_CU"),
               *str_nSU  = getenv("DARTS_NUM_SU");

    g_nCU = str_nCU != 0 ? strtoul(str_nCU,NULL,0) : 1;
    g_nSU = str_nSU != 0 ? strtoul(str_nSU,NULL,0) : 1;

	g_treeBarrier = opts.treeBarrier; 

////	if ((myRank == 0) && (opts.quiet == 0)) {
////		printf("Running problem size %d^3 per domain until completion\n", opts.nx);
////		printf("Num processors: %d\n", numRanks);
//////#if _OPENMP
//////      printf("Num threads: %d\n", omp_get_max_threads());
//////#endif
////		printf("Total number of elements: %lld\n\n", (long long int)(numRanks*opts.nx*opts.nx*opts.nx));
////		printf("To run other sizes, use -s <integer>.\n");
////		printf("To run a fixed number of iterations, use -i <integer>.\n");
////		printf("To run a more or less balanced region set, use -b <integer>.\n");
////		printf("To change the relative costs of regions, use -c <integer>.\n");
////		printf("To print out progress, use -p\n");
////		printf("To write an output file for VisIt, use -v\n");
////		printf("To run tree barrier set, use -t <integer>.\n");
////		printf("See help (-h) for more options\n\n");
////   }

   // Set up the mesh and decompose. Assumes regular cubes for now
   Int_t col, row, plane, side;
   InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

   // Build the main data structure and initialize it
   locDom = new Domain(numRanks, col, row, plane, opts.nx,
                       side, opts.numReg, opts.balance, opts.cost) ;


//debug to see region sizes
//   for(Int_t i = 0; i < locDom->numReg(); i++)
//      std::cout << "region" << i + 1<< "size" << locDom->regElemSize(i) <<std::endl;

   
	uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;

	uint64_t *timings= new uint64_t [opts.its];
    Runtime* rt;

	std::cout <<"ThreadAffinity!"<<std::endl;			
//	ThreadAffinity affin(g_nCU, g_nSU, COMPACT, TPPUSHFULL, MCSTANDARD);

	ThreadAffinity affin(g_nCU, g_nSU, COMPACT, 3,1);
	timeval start;
	gettimeofday(&start, NULL) ;

	if (affin.generateMask()){

		rt= new Runtime (&affin);		
		std::cout <<"LagrangeLeapFropTP begin!"<<std::endl;			
		int its=0;

			while((locDom->time() < locDom->stoptime()) && (locDom->cycle() < opts.its)) {
	
				TimeIncrement (*locDom);
			
				innerStart = getTime();//start time for kernel procedure
				rt->run(launch<LagrangeLeapFrogTP>(locDom,&Runtime::finalSignal)); 

				innerStop = getTime() - innerStart;
				innerAvg += innerStop;
				
				timings[its]=innerStop;
				++its;

//				if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank == 0)) {
//					printf("cycle = %d, time = %e, dt=%e\n",
//							locDom->cycle(), double(locDom->time()), double(locDom->deltatime()) ) ;
//				}
			}

	}else{
		std::cerr << "Could not generate required abstract machine -- something went wrong. :(\n";
		return EXIT_FAILURE;
	}

    outerStop = getTime() - outerStart;
    outerAvg += outerStop;

    innerAvg /= opts.its;
    outerAvg /= opts.its;

    std::cout 
              << g_nCU                   << ", "
              << g_nSU                   << ", "
              << g_nSU*(1+g_nCU)         << ", "
              << opts.its                   << ", "
              << std::setprecision(18)   
              << innerAvg                << ", "
//              << outerAvg                << "\n";
              << outerAvg                << ",";

	for(Int_t i=0;i < opts.its;++i){
		printf("%lu,",timings[i]);
	}
	printf("\n");
	delete rt;
	delete[] timings;
	delete locDom;

	
	
	// Use reduced max elapsed time
   double elapsed_time;
   timeval end;
   gettimeofday(&end, NULL) ;
   elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
   double elapsed_timeG;
   elapsed_timeG = elapsed_time;
//
   // Write out final viz file */
   if (opts.viz) {
     // DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   
     DumpToVisit(myRank) ;
   }

#ifdef TRACE
	darts::printTPRecord_v4();
#endif


//  if ((myRank == 0) && (opts.quiet == 0)) {
//     VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, opts.nx, numRanks);
//  }
   return 0 ;
}



