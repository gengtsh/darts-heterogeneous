#include "stencil.h"
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "./timer.h"

#ifdef _OPENMP
#   include <omp.h>
#endif

#define N_ROWS   1000UL
#define N_COLS   10000UL
#define N_TSTEPS  100UL
#define N_REPS     10UL

#define SIZE2D(type,dim1,dim2) sizeof(type)*(dim1)*(dim2)

static inline void usage(const char *name) 
{
    printf("USAGE: %s <n_rows> <n_cols> <n_timesteps> <n_reps> <n_CU_per_SU> <n_SU>\n", name);
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
result_is_correct(const size_t  n_rows, const size_t  n_cols, 
                  const double* restrict orig,   const double* restrict res)
{
    typedef double (*Array2D)[n_cols];
    const Array2D ORIG       = (const Array2D) orig,
                  RES        = (const Array2D) res;
    bool          is_correct =   true;
    const double  epsilon    = 0.0001;

    for (size_t i = 0; is_correct && i < n_rows; ++i) 
        for (size_t j = 0; is_correct && j < n_cols; ++j) 
            if ( !( is_correct = fabs( RES[i][j] - ORIG[i][j] ) <= epsilon ) )
                printf("Values mismatch! [%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n",
                       i, j, ORIG[i][j], RES[i][j]);

    return is_correct;
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
void 
stencil_init ( stencil_t* stencil, stencil_funptr_t code, const size_t n_rows, const size_t n_cols, size_t n_tsteps )
{
//  typedef double (*Array2D)[n_cols];
    double *dst = smalloc ( SIZE2D(double, n_rows, n_cols) ),
           *src = smalloc ( SIZE2D(double, n_rows, n_cols) );
    array2d_init ( src, n_rows, n_cols );
    memcpy ( dst, src, sizeof(double) * n_rows * n_cols );
    stencil->arg = smalloc ( SIZE2D( double, n_rows, n_cols ) );
    *(stencil->arg) = (stencil_arg_t){ 
        .dst = dst, .src = src, 
        .n_rows = n_rows, .n_cols = n_cols, 
        .n_tsteps = n_tsteps 
    };
    stencil->stencil = code;
}

/**
 * brief <stencil_destroy> deallocates the dynamically allocated data structures
 * allocated within a stencil_t structure.
 * inputs:
 * stencil -- the stencil_t data structure to deallocate
 * output: N/A
 */
void
stencil_destroy( stencil_t* stencil )
{
    sfree ( stencil->arg->dst );
    sfree ( stencil->arg->src );
    sfree ( stencil->arg );
}

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
void
stencil_arg_copy ( stencil_arg_t* dest, stencil_arg_t* sourc )
{
    *dest = (stencil_arg_t) { 
        .dst     = NULL,          .src    = NULL, 
        .n_rows  = sourc->n_rows, .n_cols = sourc->n_cols, 
        .n_tsteps = sourc->n_tsteps 
    };
    dest->src = smalloc ( SIZE2D(double,dest->n_rows,dest->n_cols) );
    dest->dst = smalloc ( SIZE2D(double,dest->n_rows,dest->n_cols) );
    memcpy ( dest->src, sourc->src, SIZE2D(double,dest->n_rows,dest->n_cols) );
    memcpy ( dest->dst, sourc->dst, SIZE2D(double,dest->n_rows,dest->n_cols) );
}

/**
 * brief <stencil_arg_dup> allocates a new argument and copies the content of
 * <sourc> into it
 * inputs:
 * sourc -- the source where to duplicate from
 *
 * output: a <stencil_arg_t*> pointer to the newly allocated memory region
 */
stencil_arg_t*
stencil_arg_dup ( stencil_arg_t* src )
{
    stencil_arg_t* dst = smalloc ( sizeof(stencil_arg_t) );
    stencil_arg_copy ( dst, src );
    return dst;
}

void
run_kernel(code_t code, stencil_t *stencil, const size_t n_reps, const char* label)
{
    loc_timer_t loc_timer = { 
        .start = {0}, .stop = {0}, .cpu = 0L, 
        .code = code, .data = (void*)stencil,
        .n_reps = n_reps 
    };

    pthread_t timer_thd;
    spthread_create ( &timer_thd, NULL /* no attributes */, timer, &loc_timer );
    if ( 0 != pthread_join ( timer_thd, NULL ) )
        fatal(EXIT_FAILURE,"pthread_join");
    timespec_diff(&loc_timer.stop, &loc_timer.start, &loc_timer.stop);
    double t = loc_timer.stop.tv_sec + (double)loc_timer.stop.tv_nsec / 1000000000.0;
    printf("%s kernel time: %.5f secs\n", label, t/n_reps);
}

int 
main(int argc, char* argv[])
{
    uint64_t    n_rows     = 0UL, 
                n_cols     = 0UL,
                n_tm_steps = 1UL,
                n_reps     = 1UL;

    switch ( argc ) {
    case 5: n_reps     = strtoul(argv[4],NULL,0);
    case 4: n_tm_steps = strtoul(argv[3],NULL,0); 
    case 3: n_rows     = strtoul(argv[1],NULL,0); 
            n_cols     = strtoul(argv[2],NULL,0);
            break;
    case 2: n_rows = n_cols = strtoul(argv[1],NULL,0);
            break;
    default: usage(argv[0]);
    }

    double *current_values = smalloc( sizeof(double) * n_rows * n_cols),
           *next_values    = smalloc( sizeof(double) * n_rows * n_cols),
           *initial_matrix = smalloc( sizeof(double) * n_rows * n_cols);

	array2d_init(initial_matrix,n_rows,n_cols); // initial 2D array

    uint64_t outer_start = 0, 
             outer_stop  = 0,
             inner_start = 0,
             inner_stop  = 0;
    double   outer_avg   = 0.0,
             inner_avg   = 0.0;

	uint64_t *timings= smalloc(sizeof(uint64_t)*n_reps);

    char *str_n_threads = getenv("OMP_NUM_THREADS");
    long  n_threads     = str_n_threads ? strtoul(str_n_threads,NULL,0) : 1L;

    outer_start = get_time();
    for (size_t i = 0; i < n_reps; ++i) {
        //memcpy(current_values, initial_matrix, sizeof(double*)*n_rows*n_cols);
        //memcpy(next_values,    initial_matrix, sizeof(double*)*n_rows*n_cols);
        init_loop(current_values,initial_matrix,n_rows,n_cols);
        init_loop(next_values,initial_matrix,n_rows,n_cols);
        inner_start = get_time();
        //stencil2D4pt_omp_v2(next_values, current_values, n_rows, n_cols, n_tm_steps);
        stencil2D4pt_omp_simd(next_values, current_values, n_rows, n_cols, n_tm_steps);
        inner_stop  = get_time() - inner_start;
        inner_avg  += inner_stop;
    
		timings[i]=inner_stop;
	}
    outer_stop  = get_time() - outer_start;
    outer_avg  += outer_stop;

    inner_avg /= n_reps;
    outer_avg /= n_reps;

    printf("%lu * %lu, %lu, %lu, %ld, %lu, %lu, %-18.2f, %-18.2f,",
           n_rows, n_cols, 0UL, n_threads, n_threads, n_tm_steps, n_reps, 
           inner_avg, outer_avg);

	for(size_t i=0;i<n_reps;++i){
		printf("%lu,",timings[i]);
	}
	printf("\n");

    sfree(current_values);
    sfree(next_values);
    sfree(initial_matrix);
	sfree(timings);
    return 0;
}
