#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include "stencil.h"
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "timer.h"
#include "conf.h"
//#include <thread.h>


#define SIZE2D(type,dim1,dim2) sizeof(type)*(dim1)*(dim2)


static inline void usage(const char *name) 
{
	printf("USAGE: %s,  <n_rows> <n_cols> <n_slices> <n_timesteps> <n_reps>\n",name);
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
//bool
//result_is_correct(const size_t  n_rows, const size_t  n_cols, 
//                  const double* __restrict__ orig,   const double* __restrict__ res)
//{
//    typedef double (*Array2D)[n_cols];
//    const Array2D ORIG       = (const Array2D) orig,
//                  RES        = (const Array2D) res;
//    bool          is_correct =   true;
//    const double  epsilon    = 0.0001;
//
//    for (size_t i = 0; is_correct && i < n_rows; ++i) 
//        for (size_t j = 0; is_correct && j < n_cols; ++j) 
//            if ( !( is_correct = fabs( RES[i][j] - ORIG[i][j] ) <= epsilon ) )
//                printf("Values mismatch! [%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",
//                       i, j, ORIG[i][j], RES[i][j]);
//
//    return is_correct;
//}


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
	    	for(size_t j = 0; is_correct && j < n_cols; ++j) 
	    		if ( !( is_correct = fabs( RES[k][i][j] - ORIG[k][i][j] ) <= epsilon ) )
	    			printf("Values mismatch! [%lu,%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",k,i, j,ORIG[k][i][j], RES[k][i][j]);

	return is_correct;
}

void print_results(const double *results, const size_t  n_rows_st,const size_t n_rows_ed, const size_t  n_cols_st, const size_t n_cols_ed,const size_t n_slices_st,const size_t n_slices_ed,const size_t n_rows, const size_t n_cols, const size_t n_slices ){
	for(size_t s=n_slices_st;s<n_slices_ed;++s){
        printf("slice %d:\n ",s);
        for(size_t k=n_rows_st;k<n_rows_ed;++k){
		    for(size_t j=n_cols_st;j<n_cols_ed;++j){
			    printf("%f,",results[s*n_rows*n_cols+k*n_cols+j]);
		    }
		    printf("\n");
	    }
		printf("\n");
    }

}

/**
 * brief <array2d_print> takes a 2D array (assumed to be allocated in contiguous
 * memory) and print it in the std output
 *
 * inputs:
 * a      -- the 2D array to print
 * n_rows -- the number of rows in <a>
 * n_cols -- the number of columns in <a>
 *
 * output: N/A
 */
void
array2d_print ( double* a, const size_t n_rows, const size_t n_cols )
{
    typedef double (*Array2D)[n_cols];
    Array2D A = (Array2D) a;
    for (size_t i = 0; i < n_rows; ++i) 
    {
        for (size_t j = 0; j < n_cols; ++j) 
		printf("%G\t", A[i][j]);
    	printf("\n");
    }
}
/**
 * brief <array2d_init> takes a 2D array (assumed to be allocated in contiguous
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
void
array2d_init ( double* a, const size_t n_rows, const size_t n_cols )
{
    typedef double (*Array2D)[n_cols];
    Array2D A = (Array2D) a;
    for (size_t i = 0; i < n_rows; ++i) 
        for (size_t j = 0; j < n_cols; ++j) 
            A[i][j] = (1+i)*j + 1.3;
}

/**
 * brief <stencil_init> takes a stencil structure and initializes its source and
 * destination arrays (including memory allocation) to arbitrary (well-formed)
 * double precision values
 *
 * inputs:
 * stencil -- the stencil_t data structure to initialize
 * n_rows  -- the number of rows contained in the 2D arrays
 * n_cols  -- the number of columns contained in the 2D arrays
 */
//void 
//stencil_init ( stencil_t* stencil, stencil_funptr_t code, const size_t n_rows, const size_t n_cols, size_t n_tsteps )
//{
////  typedef double (*Array2D)[n_cols];
//    double *dst = smalloc ( SIZE2D(double, n_rows, n_cols) ),
//           *src = smalloc ( SIZE2D(double, n_rows, n_cols) );
//    array2d_init ( src, n_rows, n_cols );
//    memcpy ( dst, src, sizeof(double) * n_rows * n_cols );
//    stencil->arg = smalloc ( SIZE2D( double, n_rows, n_cols ) );
//    *(stencil->arg) = (stencil_arg_t){ 
//        .dst = dst, .src = src, 
//        .n_rows = n_rows, .n_cols = n_cols, 
//        .n_tsteps = n_tsteps 
//    };
//    stencil->stencil = code;
//}

/**
 * brief <stencil_destroy> deallocates the dynamically allocated data structures
 * allocated within a stencil_t structure.
 * inputs:
 * stencil -- the stencil_t data structure to deallocate
 * output: N/A
 */
//void
//stencil_destroy( stencil_t* stencil )
//{
//    sfree ( stencil->arg->dst );
//    sfree ( stencil->arg->src );
//    sfree ( stencil->arg );
//}

/**
 * brief <stencil_arg_copy> copies the content of <sourc> into <dest>, after
 * having allocated the source and destination arrays of <dest> to match those
 * used in <sourc>.
 *
 * inputs:
 * dest  -- the target stencil_arg_t structure
 * sourc -- the source where to copy from
 *
 * output: N/A
 */
//void
//stencil_arg_copy ( stencil_arg_t* dest, stencil_arg_t* sourc )
//{
//    *dest = (stencil_arg_t) { 
//        .dst     = NULL,          .src    = NULL, 
//        .n_rows  = sourc->n_rows, .n_cols = sourc->n_cols, 
//        .n_tsteps = sourc->n_tsteps 
//    };
//    dest->src = smalloc ( SIZE2D(double,dest->n_rows,dest->n_cols) );
//    dest->dst = smalloc ( SIZE2D(double,dest->n_rows,dest->n_cols) );
//    memcpy ( dest->src, sourc->src, SIZE2D(double,dest->n_rows,dest->n_cols) );
//    memcpy ( dest->dst, sourc->dst, SIZE2D(double,dest->n_rows,dest->n_cols) );
//}

/**
 * brief <stencil_arg_dup> allocates a new argument and copies the content of
 * <sourc> into it
 * inputs:
 * sourc -- the source where to duplicate from
 *
 * output: a <stencil_arg_t*> pointer to the newly allocated memory region
 */
//stencil_arg_t*
//stencil_arg_dup ( stencil_arg_t* src )
//{
//    stencil_arg_t* dst = smalloc ( sizeof(stencil_arg_t) );
//    stencil_arg_copy ( dst, src );
//    return dst;
//}


/**
 * brief <strtoi> uses the function <strtol> in order to convert a char
 * into an integer. This in a secure way
 * inputs:
 * str -- String with the number
 * endptr -- ending pointer
 * base -- base of the input number wrote in the str
 *
 * output a <int> with the value of the input str.
 */
int strtoi (const char* str, char** endptr, int base)
{
	long long_val;

	errno = 0;
	long_val = strtol (str, endptr, base);
	if (errno)
	{
		printf("[ERROR] Something went wrong in str2int translation \n");
		exit(0);
	}
	if ((long)  long_val < INT_MIN || long_val > (long) INT_MAX)
	{
		printf("[ERROR] Limits of int exceeded\n ");
		exit(0);
	}
	return (int) long_val;
}



//void
//StencilSeq ( double* __restrict__ dst,    double* __restrict__ src, 
//               const size_t     n_rows, const size_t     n_cols,
//               const size_t     n_tsteps )
//{
//    typedef double (*Array2D)[n_cols];
//    volatile Array2D DST = (Array2D) dst,
//            SRC = (Array2D) src;
//    for (size_t ts = 0; ts < n_tsteps; ++ts) {
//		for (size_t i = 1; i < n_rows-1; ++i) {
//            for (size_t j = 1; j < n_cols-1; ++j) {
//                DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1])/5.5;
//                //DST[i][j] = (SRC[i-1][j] + SRC[i+1][j] + SRC[i][j-1] + SRC[i][j+1])/4;
//            }
//        }
//		SWAP_PTR(&DST,&SRC);
//    }
//}

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
                nSlices  = 0,
                nTmSteps = 0,
                nReps    = 0;

    
    long  nthreads     = 1;
	nTmSteps = nReps = 1;

    switch ( argc ) {
    
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

	double* OriginalMatrix = smalloc( sizeof(double) * nRows * nCols*nSlices)  ;
	double* InitialMatrix ; 	
    double* NewMatrix  ;
    cudaMallocHost((void**)&InitialMatrix,sizeof(double)*nRows*nCols*nSlices) ;
    cudaMallocHost((void**)&NewMatrix,sizeof(double)*nRows*nCols*nSlices) ;
  
    init3Darray(OriginalMatrix,nRows,nCols,nSlices);// initial 2D array

    uint64_t innerStart = 0, 
             outerStart = 0,
             innerStop  = 0,
             outerStop  = 0;
    double   innerAvg   = 0.0,
             outerAvg   = 0.0;

	uint64_t *timings= smalloc(sizeof(uint64_t)*nReps);
    outerStart = getTime();
	
		
	for(size_t i=0; i<nReps;++i){

        memcpy(InitialMatrix, OriginalMatrix, sizeof(double*)*nRows*nCols*nSlices);
        memcpy(NewMatrix, OriginalMatrix, sizeof(double*)*nRows*nCols*nSlices);

		innerStart = getTime();//start time for kernal procedure
		
		stencil3D7pt_gpu( NewMatrix, InitialMatrix, nRows,nCols,nSlices,nTmSteps);
		innerStop = getTime() - innerStart;
		innerAvg += innerStop;
		timings[i]=innerStop;
	}
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


    printf("%lu * %lu * %lu, %lu, %lu, %ld, %lu, %lu, %-18.2f, %-18.2f,",
           nRows, nCols,nSlices, 0UL, nthreads,nthreads, nTmSteps,nReps,innerAvg, outerAvg);

	for(size_t i=0;i<nReps;++i){
		printf("%lu,",timings[i]);
	}
	printf("\n");

	//for verification:




#ifdef VERIFICATION
	double* SeqMatrix = smalloc( sizeof(double) * nRows * nCols*nSlices);
    double* SeqOutMatrix = smalloc( sizeof(double) * nRows * nCols*nSlices); 
	
    memcpy(SeqMatrix, OriginalMatrix, sizeof(double*)*nRows*nCols*nSlices);
    memcpy(SeqOutMatrix, OriginalMatrix, sizeof(double*)*nRows*nCols*nSlices);
	StencilSeq ( SeqOutMatrix,SeqMatrix,nRows,nCols,nSlices,nTmSteps); 

	double *seqOut, *dartsOut;
	
	if(nTmSteps%2){
		seqOut=SeqOutMatrix;
	}else{
		seqOut = SeqMatrix;
	}


	dartsOut=NewMatrix;
	if(result_is_correct(nRows,nCols,nSlices,seqOut,dartsOut)){
		printf("success!\n");
	}else{
		printf("fail! \n");
	}

#ifdef CUDA_CUDA_DEBUG
	int kk = 0;
	int ttk = 10;
	int jj =0; 
	int ttj =10;
    int ss =0;
    int tts=4;
	printf("Original Matrix\n");
	print_results(OriginalMatrix,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
	printf("Seq Matrix\n");
	print_results(seqOut,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
	printf("Out Matrix\n");
	print_results(dartsOut,kk,ttk,jj,ttj,ss,tts,nRows,nCols,nSlices);
#endif
#endif
	

#ifdef VERIFICATION
	sfree(SeqMatrix);
	sfree(SeqOutMatrix);
#endif
    sfree(OriginalMatrix);
    
    cudaFreeHost(InitialMatrix);
    cudaFreeHost(NewMatrix);

	return 0;
}

