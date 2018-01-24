#define __STDC_FORMAT_MACROS
#include "Stencil2DRowDecomposition.h"
//#include <math.h>
//#include <inttypes.h>
#include <stdint.h>

//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
//#include <iostream>


void Stencil2DRowLoopCopyUp::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	
	uint64_t Id = getID();

	RESET(copyUp[Id]);	

//	uint64_t *ts = FRAME(tSteps);
//	pthread_mutex_lock(&mutex);
//	std::cout<<"U"<<Id<<",ts:"<<ts[Id]<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"P"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
//	fprintf(stderr, "U%ju,ts:%ju\n",Id,ts[Id]);

	uint64_t nRowsCut = FRAME(nRowsCut);	//rows decomposition number	
	double *initial = FRAME(initial); //matrix pointer initial Matrix[M][n]
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	uint64_t bm=n_rows - 2;//blockM
	//uint64_t bn=n_cols - 2;//blockN
	double *share = FRAME(shareRows+Id*2*n_cols);
	const uint64_t rows_ini = bm / nRowsCut; // initially, the total number of rows in every nRowsCut
	uint64_t pos = Id*rows_ini * n_cols; 
	
	copyLine_stencil2d(pos,initial,share,n_cols);//copy shared upper line between blocks (up+lower lines)
	
	if(Id==0){
		SYNC(compute[Id]);
		SYNC(compute[Id]);
	}else{
		SYNC(compute[Id]);
		SYNC(compute[Id-1]);
	}

//	bool *finalize = FRAME(finalize);
//	if( finalize[Id] == true){
//		SYNC(sync);
//		EXIT_TP();
//	}

	EXIT_TP();

}

void Stencil2DRowLoopCopyDown::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	
	uint64_t Id = getID();

	RESET(copyDown[Id]);	
	
//	uint64_t *ts = FRAME(tSteps);
//	pthread_mutex_lock(&mutex);
//	std::cout<<"D"<<Id<<",ts:"<<ts[Id]<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"D"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
//	fprintf(stderr, "D%ju,ts:%ju\n",Id,ts[Id]);

	uint64_t nRowsCut = FRAME(nRowsCut);	
	double *initial = FRAME(initial); //matrix pointer initial Matrix[M][n]
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	//uint64_t bp= n_cols + 1;//block position
	uint64_t bm=n_rows - 2;//blockM
	//uint64_t bn=n_cols - 2;//blockN
	double *share = FRAME(shareRows+Id*2*n_cols);
	const uint64_t rows_ini = bm / nRowsCut; // initially, the total number of rows in every nRowsCut
	uint64_t bm_final = ((Id==(nRowsCut-1))? (bm%nRowsCut):0) + rows_ini;
	uint64_t pos = Id*rows_ini* n_cols+(bm_final+1)*n_cols; 

	copyLine_stencil2d(pos,initial,share+n_cols,n_cols);//copy shared lines between blocks (up+lower lines)

	if((Id==nRowsCut-1)){
		SYNC(compute[Id]);
		SYNC(compute[Id]);
	}else{
		SYNC(compute[Id]);
		SYNC(compute[Id+1]);
	}

	
//	bool *finalize = FRAME(finalize);
//	if( finalize[Id] == true){
//		SYNC(sync);
//		EXIT_TP();
//	}

	EXIT_TP();


}

void Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);

	uint64_t Id = getID();	
	bool *finalize = FRAME(finalize);
	if( finalize[Id] == true){
		SYNC(sync);
		EXIT_TP();
	}
	RESET(compute[Id]);	
	uint64_t *ts = FRAME(tSteps);

//	pthread_mutex_lock(&mutex);
//	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"C"<<Id<<",ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
//	fprintf(stderr, "C%ju,ts:%ju\n",Id,ts[Id]);
	
	--(ts[Id]);
	
	if(ts[Id]==0){
		finalize[Id]=true;
	}

	uint64_t nRowsCut = FRAME(nRowsCut);	
	double *initial = FRAME(initial); //matrix pointer initial Matrix[M][n]
	const uint64_t n_rows = FRAME(nRows); // matrix M row
	const uint64_t n_cols = FRAME(nCols); // Matrix N column
	uint64_t bp= n_cols + 1;//block position
	uint64_t bm=n_rows - 2;//blockM
	uint64_t bn=n_cols - 2;//blockN
	double *share = FRAME(shareRows+Id *2*n_cols);
	const uint64_t rows_ini = bm / nRowsCut; // initially, the total number of rows in every nRowsCut
	uint64_t bm_final = ((Id ==(nRowsCut-1))? (bm%nRowsCut):0) + rows_ini;
	uint64_t pos = Id*rows_ini * n_cols; 
	
	computeInner_stencil2d(bp+pos, initial,share,bm_final,bn,n_rows,n_cols);

	if(nRowsCut==1){
		SYNC(copyUp[Id]);
		SYNC(copyDown[Id]);
	}else{
		if(Id ==0){
			SYNC(copyUp[Id]);
			SYNC(copyUp[Id+1]);
		}else if(Id==(nRowsCut-1)){
			SYNC(copyDown[Id]);
			SYNC(copyDown[Id-1]);
		}else {
			SYNC(copyUp[Id+1]);
			SYNC(copyDown[Id-1]);
		}
	}

	SYNC(copyUp[Id]);
	SYNC(copyDown[Id]);
	

	EXIT_TP();

}

void
Stencil2DRowSync::fire(void)
{
    LOAD_FRAME(Stencil2DRowDecomposition);

	ThreadedProcedure* PTP = getTP();
	PTP->setRef(1);
	SIGNAL(signalUp);
    EXIT_TP();
}

