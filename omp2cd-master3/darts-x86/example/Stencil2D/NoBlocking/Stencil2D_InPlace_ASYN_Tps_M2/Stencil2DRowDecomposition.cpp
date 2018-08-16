#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include "Stencil2DRowDecomposition.h"
#include <math.h>
#include <inttypes.h>

//#include <pthread.h>
//pthread_mutex_t mutex;
//#include <sstream>
//#include <iostream>

void Stencil2DRowCheck::fire(void){
	LOAD_FRAME(Stencil2DRowDecomposition);
	uint64_t nRowsCut = FRAME(nRowsCut);	
	uint64_t Id = getID();
	RESET(check[Id]);

	
	if(Id==0){
		SYNC(copyUp[0]);
	}else if(Id==1){
		SYNC(copyDown[nRowsCut-1]);
	}

//	uint32_t nTp = FRAME(nTp);
//	pthread_mutex_lock(&mutex);
//	std::cout<<"RC"<<Id<<",nTP:"<<nTp<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"P"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
	
	EXIT_TP();
}

void Stencil2DRowLoopCopyUp::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	
	uint64_t Id = getID();

	bool *finalize = FRAME(finalize);
	if( finalize[Id] == true){
		SYNC(sync);
		EXIT_TP();
	}

	RESET(copyUp[Id]);	

//	ThreadedProcedure* PTP = getTP();
//	PTP->incRef();
	
//	uint32_t nTp = FRAME(nTp);
//	uint64_t *ts = FRAME(tSteps);
//	pthread_mutex_lock(&mutex);
//	std::cout<<"U"<<Id<<",ts:"<<ts[Id]<<",nTP:"<<nTp<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"P"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();


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

	EXIT_TP();

}

void Stencil2DRowLoopCopyDown::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);
	
	uint64_t Id = getID();
	bool *finalize = FRAME(finalize);
	if( finalize[Id] == true){
		SYNC(sync);
		EXIT_TP();
	}
	RESET(copyDown[Id]);	

//	ThreadedProcedure* PTP = getTP();
//	PTP->incRef();

//	uint32_t nTp = FRAME(nTp);
//	uint64_t *ts = FRAME(tSteps);
//	pthread_mutex_lock(&mutex);
//	std::cout<<"D"<<Id<<",ts:"<<ts[Id]<<",nTP:"<<nTp<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"D"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
	
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

	EXIT_TP();


}

void Stencil2DRowLoop::fire(void) 
{
	LOAD_FRAME(Stencil2DRowDecomposition);

	uint64_t Id = getID();	
	bool *finalize = FRAME(finalize);
//	if( finalize[Id] == true){
//		SYNC(sync);
//		EXIT_TP();
//	}
	
	RESET(compute[Id]);	
//	ThreadedProcedure* PTP = getTP();
//	PTP->incRef();
	
	uint64_t *ts = FRAME(tSteps);
	
	uint32_t nTp = FRAME(nTp);

//	pthread_mutex_lock(&mutex);
//	std::cout<<"C"<<Id<<",ts:"<<ts[Id]<<",nTP:"<<nTp<<std::endl;
//	pthread_mutex_unlock(&mutex);
//	std::stringstream ss;
//	ss <<"C"<<Id<<"ts:"<<ts[Id]<<"\n";
//	std::cout <<ss.str();
		
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

//	printf("Line CPU %d: %" PRIu64 " - %" PRIu64 "\n",(int)getID(),(BlockPosition+pos_final),BlockM_final);
    
	computeInner_stencil2d(bp+pos, initial,share,bm_final,bn,n_rows,n_cols);

	--(ts[Id]);

	//uint32_t nTp = FRAME(nTp);

	if(ts[Id]==0){
		//SYNC(sync);
		finalize[Id]=true;
	}


	if(nRowsCut==1){
	//	SYNC(copyUp[Id]);
	//	SYNC(copyDown[Id]);
			
		if(nTp==0){
			SYNC(checkTP[nTp]);
		}
		if(nTp==(g_nSU-1)){
			SYNC(checkTP[nTp+1]);
		}
		SYNC(checkTP[nTp]);
		SYNC(checkTP[nTp+1]);
		

	}else{
		if(Id ==0){
			//SYNC(copyUp[Id]);
			SYNC(copyUp[Id+1]);
		//	if(ts[Id]!=0){
				if(nTp==0){
					SYNC(checkTP[nTp]);
				}
				SYNC(checkTP[nTp]);
		//	}
		}else if(Id==(nRowsCut-1)){
			//SYNC(copyDown[Id]);
			SYNC(copyDown[Id-1]);
		//	if(ts[Id]!=0){
				if(nTp==(g_nSU-1)){
					SYNC(checkTP[nTp+1]);
				}
				SYNC(checkTP[nTp+1]);
		//	}
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
//	pthread_mutex_lock(&mutex);
//	std::cout<<"RowSync,nTP:"<<FRAME(nTp)<<std::endl; 
//	pthread_mutex_unlock(&mutex);
	
	SIGNAL(signalUp);
	EXIT_TP();
}


