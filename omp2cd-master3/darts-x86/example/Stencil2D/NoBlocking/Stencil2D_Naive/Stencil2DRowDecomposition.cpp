#define __STDC_FORMAT_MACROS
#include <cstdint>
#include "Stencil2DRowDecomposition.h"
#include "Stencil2DKernel.h"
//#include <cmath>
#include <inttypes.h>

void 
Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);

	double        *Initial  = FRAME(Initial), //matrix pointer initial Matrix[M][n]
	              *New      = FRAME(New);
	const uint64_t n_rows   = FRAME(nRows),   // matrix M row
	               n_cols   = FRAME(nCols);   // Matrix N column
    uint64_t       nRowsCut = FRAME(nRowsCut),
                   bp       = n_cols + 1,     //block position begin in the [row1,cols1]
                   BlockM   = n_rows-2,
                   BlockN   = n_cols-2,
                   Id	    = getID();	
	const uint64_t rows_ini = BlockM / nRowsCut; // initially,the total number of rows in every nRowsCut

	RESET(compute[Id]);
//  blockM from 1 to (nRowsCut-1), is the same (nr_initial), 
//  but the last one is different with the former, it needs add the remains.
	uint64_t BlockM_final = ((Id==(nRowsCut-1))? (BlockM%nRowsCut):0) + rows_ini;
	uint64_t pos = Id*rows_ini * n_cols;   //beginning position of every nRowsCut 

	computeInner_stencil2d(bp+pos,Initial,New, BlockM_final,BlockN,n_cols);

	SYNC(syncSwap);

	EXIT_TP();
}

void
Stencil2DRowSyncSwap::fire(void)
{
	LOAD_FRAME(Stencil2DRowDecomposition);

    if ( FRAME(timeStep)-- > 0 ) {
		RESET(syncSwap);
        SWAP_PTR(&FRAME(New),&FRAME(Initial));
        uint64_t nRowsCut = FRAME(nRowsCut);	
		for(size_t i = 0; i < nRowsCut; ++i)
			SYNC(compute[i]);
    } else {
		SIGNAL(signalUp);
    }
	
	EXIT_TP();
}
