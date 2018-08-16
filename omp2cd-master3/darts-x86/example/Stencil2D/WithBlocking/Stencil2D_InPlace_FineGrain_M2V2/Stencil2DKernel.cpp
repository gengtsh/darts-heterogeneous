

#include <stdint.h>
#include <stdlib.h>
#include "Stencil2DKernel.h"
#include "Stencil2D_main.h"

pthread_mutex_t mutex2;

/**
* Naive 4pt stencil for 2D array
*/
void stencil2d_seq(double *dst,double *src,const uint64_t n_rows,const uint64_t n_cols,uint64_t n_tsteps){
    
	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
		for (size_t i = 1; i < n_rows-1; ++i) {
            for (size_t j = 1; j < n_cols-1; ++j) {
                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
			}
		}
		SWAP_PTR(&DST,&SRC);
	}

}

uint32_t computeRowDecomposition(const uint64_t n_rows, const uint64_t n_cols){

	uint64_t nTotalSize = (n_rows-2)*(n_cols-2);//the total number of pixel in inner matrix (the whole matrix - left_edge - right_edge - upper_edge - down_edge)
	uint64_t nRowsCut = nTotalSize/TOTAL_TILE_SZ ;	
	uint64_t nThreads = (g_nCU+1)*g_nSU*N_THREADS;//(nCU+1)*nSU is available cores, N_THREADS are number of threads per cores

	uint32_t nRowsCut_final = (nRowsCut<1)?1:nThreads;
//	uint32_t nRowsCut_final = (nRowsCut<1)?1:((nRowsCut<nThreads)?nRowsCut:nThreads);
	return nRowsCut_final;
}



void computeInner_stencil2d(uint64_t BlockPosition,double* src,double *dst, uint64_t BlockM, uint64_t BlockN,const uint64_t N){

    typedef double (*Array2D)[N];
    Array2D SRC = (Array2D) src,
            DST = (Array2D) dst;

	uint64_t i_begin = BlockPosition / N;
	uint64_t j_begin = BlockPosition % N;
	uint64_t i_end = i_begin + BlockM-1;
	uint64_t j_end = j_begin + BlockN-1;
	
	for(uint64_t i = i_begin;i <= i_end; i++)
		for(uint64_t j = j_begin;j <= j_end; j++)
			DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1] )/5.5;
	
	return;
}


void computeInner_stencil2dbk(double *Initial,double *share,uint64_t nRowsChunk,uint64_t nColsChunk,const uint64_t nRows, const uint64_t nCols){

	if(nRowsChunk > nRows){
		std::cout << "Inner error \n"<<std::endl;
		return;
	}

    int n_cols_bk = GRID_TILECPU_X;
    int n_rows_bk = GRID_TILECPU_Y;
    int n_cols_bk1;
    int n_rows_bk1;

    int n_cols_bk2;
    int n_rows_bk2;

    n_cols_bk1 = n_cols_bk+1;
    n_rows_bk1 = n_rows_bk+1;

    n_cols_bk2 = n_cols_bk+2;
    n_rows_bk2 = n_rows_bk+2;

    double *upper   = new double[n_cols_bk2];
	double *lower   = new double[n_cols_bk2];
	double *current =new double[n_cols_bk2];
    double *mid     =new double[n_cols_bk2];

    double *left =new double[n_rows_bk2];
    double *right=new double[n_rows_bk2];
    
    int nRowsCK;
    int nColsCK;
    int nColsCKM1;
    int64_t d_pos;
    int64_t s_pos;
    bool IsLastRBk;
    bool IsLastCBk;
    bool IsLeftCBk;
    int idxrbk;
    int idx;

#ifdef STENCIL2DDEBUG
	pthread_mutex_lock(&mutex2);
    std::cout<<"nRowsChunk: "<<nRowsChunk<<std::endl;
    std::cout<<"nColsChunk: "<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif  

    
    for(int i=0;i<nRowsChunk;i=i+n_rows_bk){

        idx = i+n_rows_bk;
        IsLastRBk = idx > nRowsChunk;
        nRowsCK = IsLastRBk? (nRowsChunk-i-1):n_rows_bk1;

#ifdef STENCIL2DDEBUG

	    pthread_mutex_lock(&mutex2);
        std::cout<<"i:"<<i<<", nRowsCK: "<<nRowsCK<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif    
        
#ifdef STENCIL2DDEBUG
//        std::cout<<"share :"<<std::endl;
//        for(int k=0;k<n_cols_bk1;k++){
//            std::cout<<share[k] <<",";
//        }
//        std::cout<<"\n"<<std::endl;
//        std::cout<<"share + nCols:"<<std::endl;
//        for(int k=0;k<n_cols_bk1;k++){
//            std::cout<<share[nCols+k] <<",";
//        }
//        std::cout<<"\n"<<std::endl;
#endif
        
        for (int j=0;j<nColsChunk;j=j+n_cols_bk){
            d_pos = i*nCols+ j;
            s_pos = j;
            IsLastCBk = j+n_cols_bk >= nColsChunk;
            IsLeftCBk = (j==0);
            nColsCK = IsLastCBk? (nColsChunk-j-1):n_cols_bk1;
            idxrbk =  IsLastRBk? 0: idx;
            nColsCKM1 = nColsCK-1;

#ifdef STENCIL2DDEBUG
            pthread_mutex_lock(&mutex2);
            std::cout<<"j:"<<j<<", nColsCKM1: "<<nColsCKM1<<std::endl;
	        pthread_mutex_unlock(&mutex2);
#endif    
	        memcpy(mid,Initial+idxrbk*nCols+j+1,sizeof(double)*(nColsCKM1));
            computeInner_stencil2dv2(IsLeftCBk,IsLastRBk,Initial+d_pos,share+s_pos,nRowsCK,nColsCK,nRows, nCols,upper,lower,current,left,right);
	        memcpy(share+j+1,mid,sizeof(double)*(nColsCKM1));
            SWAP_PTR(&left,&right);
#ifdef STENCIL2DDEBUG
//            std::cout<<"update mid:"<<std::endl;
//            for(int k=0;k<n_cols_bk1;k++){
//                std::cout<<mid[k] <<",";
//            }
//            std::cout<<"\n"<<std::endl;
//
//            std::cout<<"update share:"<<std::endl;
//            for(int k=0;k<n_cols_bk1;k++){
//                std::cout<<share[i+k] <<",";
//            }
//            std::cout<<"\n"<<std::endl;

#endif
        }
    }

	delete [] upper;
	delete [] lower;
	delete [] current;
    delete [] left;
    delete [] right;
}

void computeInner_stencil2dv2(bool IsLeftCBk,bool IsLastRBk,double *Initial,double *share,uint64_t nRowsChunk,uint64_t nColsChunk,const uint64_t nRows, const uint64_t nCols,double *upper, double *lower, double *current,double *left,double *right){

    int nColsChunkP1 = nColsChunk+1;
	memcpy(current,share,sizeof(double)*(nColsChunkP1));
	memcpy(lower,Initial+nCols,sizeof(double)*(nColsChunkP1));
    int nColsChunkM1 = nColsChunk-1;
    int st;



    for(uint64_t i=1;i<nRowsChunk;++i){
		memcpy(upper,current,sizeof(double)*(nColsChunkP1));
		memcpy(current,lower,sizeof(double)*(nColsChunkP1));
	    double *LOWER = (i==(nRowsChunk-1)&&IsLastRBk)?(share+nCols):(Initial+(i+1)*nCols);	
        memcpy(lower,LOWER,sizeof(double)*(nColsChunkP1));


#ifdef STENCIL2DDEBUG
//        if (i==nRowsChunk-1){
//            int jj =1;
//            std::cout<<"jj: "<<jj<<",left: "<<current[jj-1]<<",right: "<<current[jj+1]<<",upper: "<<upper[jj]<<",lower:"<<lower[jj]<<std::endl;
//            std::cout<<"IsLastRBk: "<<IsLastRBk<<std::endl;
//            std::cout<<"lower: "<<std::endl;
//            for(int k=1;k<nColsChunkP1;++k){
//               std::cout<<lower[k]<<",";
//            }
//            std::cout<<"\n"<<std::endl;
//
//        }
#endif

        right[i] = current[nColsChunkM1];
	    if (IsLeftCBk){
            st=1;
        }else{
            st=2;
		    Initial[i*nCols+1]=(upper[1]+left[i]+current[2]+lower[1])/5.5;
        }

        for(uint64_t j=st;j<nColsChunk;++j){
			Initial[i*nCols+j]=(upper[j]+current[j-1]+current[j+1]+lower[j])/5.5;

        }
	}
        
    
#ifdef STENCIL2DDEBUG
//        std::cout<<"Left:"<<std::endl;
//        for(int k=1;k<nRowsChunk;k++){
//            std::cout<<left[k]<<",";
//        }
//        std::cout<<"\n"<<std::endl;
//        std::cout<<"Right:"<<std::endl;
//        for(int k=1;k<nRowsChunk;k++){
//            std::cout<<right[k]<<",";
//        }
//        std::cout<<" \n"<<std::endl;
#endif
	return;
}



void computeInner_stencil2d(uint64_t BlockPosition,double *Initial,double *share,uint64_t BlockM,uint64_t BlockN,const uint64_t InitialM, const uint64_t InitialN){

	uint64_t BlockNPlus2=BlockN+2;
	double *upper = new double[BlockNPlus2];//upper is used to store current line which can be used for next line computing
	double *lower = new double[BlockNPlus2];
	double *current=new double[BlockNPlus2];
	
	if(BlockM > InitialM){
		std::cout << "Inner error \n"<<std::endl;
		return;
	}

	memcpy(current,share,sizeof(double)*(BlockNPlus2));
	memcpy(lower,Initial+BlockPosition-1,sizeof(double)*(BlockNPlus2));
	for(uint64_t i=0;i<BlockM;++i){
		memcpy(upper,current,sizeof(double)*(BlockNPlus2));
		memcpy(current,lower,sizeof(double)*(BlockNPlus2));
		double *LOWER = (i==(BlockM-1))? (share+BlockNPlus2):(Initial+BlockPosition+(i+1)*InitialN-1);
		memcpy(lower,LOWER,sizeof(double)*(BlockNPlus2));
		for(uint64_t j=0;j<BlockN;++j){
			Initial[BlockPosition+i*InitialN+j]=(upper[j+1]+current[j]+current[j+2]+lower[j+1])/5.5;
		}
	}
	delete [] upper;
	delete [] lower;
	delete [] current;
	return;
}

void computeBlock_stencil25(double *dst,double *src,size_t n_rows,size_t n_cols,size_t n_rows_ck,size_t n_cols_ck){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst;
	Array2D	SRC = (Array2D) src;

	for (size_t i = 1; i < n_rows_ck; ++i) {
        for (size_t j = 1; j < n_cols_ck; ++j) {
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
        }
    }
}


void computeInner_stencil2d_v2(double *dst,double *src,size_t rpos2,size_t n_cols){

	typedef double (*Array2D)[n_cols];
	Array2D DST = (Array2D) dst,
			SRC = (Array2D) src;

	for (size_t i = 0; i < rpos2; ++i) {
	    for (size_t j = 0; j < n_cols-2; ++j) {
	        DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1]) / 5.5;
	    }
	}
}



void computeInner_stencil2dbkv3(double *dst,double *src,size_t nRowsChunk,size_t nColsChunk,const uint64_t nRows,const uint64_t nCols){

    int n_cols_bk = GRID_TILECPU_X;
    int n_rows_bk = GRID_TILECPU_Y;
    int n_cols_bk1;
    int n_rows_bk1;


    n_cols_bk1 = n_cols_bk+1;
    n_rows_bk1 = n_rows_bk+1;

    
    int nRowsCK;
    int nColsCK;
    int64_t pos;
    bool IsLastRBk;
    bool IsLastCBk;
    bool IsLeftCBk;
    int idx;

#ifdef STENCIL2DDEBUG
	pthread_mutex_lock(&mutex2);
    std::cout<<"nRowsChunk: "<<nRowsChunk<<std::endl;
    std::cout<<"nColsChunk: "<<nColsChunk<<std::endl;
	pthread_mutex_unlock(&mutex2);
#endif  
    
    for(int i=0;i<nRowsChunk;i=i+n_rows_bk){

        idx = i+n_rows_bk;
        IsLastRBk = idx > nRowsChunk;
        nRowsCK = IsLastRBk? (nRowsChunk-i-1):n_rows_bk1;

#ifdef STENCIL2DDEBUG
	    pthread_mutex_lock(&mutex2);
        std::cout<<"i:"<<i<<", nRowsCK: "<<nRowsCK<<std::endl;
	    pthread_mutex_unlock(&mutex2);
#endif    
        
        for (int j=0;j<nColsChunk;j=j+n_cols_bk){
            pos = i*nCols+ j;
            IsLastCBk = j+n_cols_bk >= nColsChunk;
            nColsCK = IsLastCBk? (nColsChunk-j-1):n_cols_bk1;

            computeBlock_stencil25(dst+pos,src+pos,nRows,nCols,nRowsCK,nColsCK);
        }
    }

}

void StencilSeq ( double* __restrict__ dst,    double* __restrict__ src,const size_t     n_rows, const size_t     n_cols,const size_t     n_tsteps ){
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

void copyLine_stencil2d(uint64_t BlockPosition,double *Initial,double *Copy,const uint64_t InitialN){

	uint64_t i = BlockPosition/InitialN;

	memcpy(Copy,Initial+i*InitialN,sizeof(double)*InitialN);
	
	return;
}



void print_results(const double *results, const size_t  n_rows_st,const size_t n_rows_ed, const size_t  n_cols_st, const size_t n_cols_ed,const size_t n_rows, const size_t n_cols){
        for(size_t k=n_rows_st;k<n_rows_ed;++k){
		    for(size_t j=n_cols_st;j<n_cols_ed;++j){
			    //std::cout<<std::setw(18)<< results[k*n_cols+j]<<",";
			    std::cout<< results[k*n_cols+j]<<",";
		    }
		    std::cout<<"\n";
	    }
		std::cout<<"\n";
}
