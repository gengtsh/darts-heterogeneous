#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <unistd.h>

#include "StencilTP.h"
#include "stencil.h"
#include <cassert>
//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
#include <iostream>

void 
Stencil2D4ptGpuCD::fire(void) 
{
	LOAD_FRAME(StencilTP);

	double *src  = FRAME(Initial); //matrix pointer initial Matrix[M][n]
	double *dst = FRAME(New);
	size_t n_rows   = FRAME(nRows);   // matrix M row
	size_t n_cols   = FRAME(nCols);   // Matrix N column
	size_t ts = FRAME(timeStep);
	std::cout<<"Invoke Cuda"<<std::endl;	
	stencil2D4pt_gpu( dst, src,n_rows,n_cols,ts);
	SYNC(sync);
	EXIT_TP();
}


void
SyncCD::fire(void)
{
    std::cout<<"Sync!"<<std::endl;
	LOAD_FRAME(StencilTP);
	SIGNAL(signalUp);
    EXIT_TP();
}





