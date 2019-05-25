#include "main.output.darts.h"
using namespace darts;
using namespace std;
static inline void usage(const char* name);
/*Function: usage, ID: 1*/
static inline void usage(const char* name)
{
    /*usage:1*/
    /*CompoundStmt:820*/
    printf("USAGE: %s <n_rows> <n_cols> <n_timesteps> <n_reps> <n_CU_per_SU> <n_SU>\n", name);
    exit(0);
}
/*Function: result_is_correct, ID: 2*/
bool result_is_correct(const size_t n_rows, const size_t n_cols, const double* restrict orig,
    const double* restrict res, const double* restrict init)
{
    /*result_is_correct:2*/
    /*CompoundStmt:823*/
    typedef double(*Array2D)[n_cols];
    const Array2D ORIG = (const Array2D)orig, RES = (const Array2D)res;
    bool is_correct = 1;
    const double epsilon = 0.10000000000000001;
    for (size_t i = 0; is_correct && i < n_rows; ++i)
        for (size_t j = 0; is_correct && j < n_cols; ++j) {
            if (!(is_correct = fabs(RES[i][j] - ORIG[i][j]) <= epsilon))
                printf("Values mismatch! [%lu,%lu]\tORIG = %5.5f != RES = %5.5f\n", i, j,
                    ORIG[i][j], RES[i][j]);
        }
    return is_correct;
}
/*Function: array2d_init, ID: 3*/
void array2d_init(double* a, const size_t n_rows, const size_t n_cols)
{
    /*array2d_init:3*/
    /*CompoundStmt:851*/
    typedef double(*Array2D)[n_cols];
    Array2D A = (Array2D)a;
    for (size_t i = 0; i < n_rows; ++i)
        for (size_t j = 0; j < n_cols; ++j)
            A[i][j] = (10 * i * i) + j * 2 + 1.3;
}
/*Function: stencil_init, ID: 4*/
void stencil_init(stencil_t* stencil, stencil_funptr_t code, const size_t n_rows,
    const size_t n_cols, size_t n_tsteps)
{
    /*stencil_init:4*/
    /*CompoundStmt:872*/
    typedef double(*Array2D)[n_cols];
    double *dst = (double *)smalloc(sizeof(double) * (n_rows) * (n_cols)),
           *src = (double *)smalloc(sizeof(double) * (n_rows) * (n_cols));
    array2d_init(src, n_rows, n_cols);
    memcpy(dst, src, sizeof(double) * n_rows * n_cols);
    stencil->arg = (stencil_arg_t*)smalloc(sizeof(double) * (n_rows) * (n_cols));
    *(stencil->arg) = (stencil_arg_t){
        .dst = dst, .src = src, .n_rows = n_rows, .n_cols = n_cols, .n_tsteps = n_tsteps
    };
    stencil->stencil = code;
}
/*Function: stencil_destroy, ID: 5*/
void stencil_destroy(stencil_t* stencil)
{
    /*stencil_destroy:5*/
    /*CompoundStmt:909*/
    sfree(stencil->arg->dst);
    sfree(stencil->arg->src);
    sfree(stencil->arg);
}
/*Function: stencil_arg_copy, ID: 6*/
void stencil_arg_copy(stencil_arg_t* dest, stencil_arg_t* sourc)
{
    /*stencil_arg_copy:6*/
    /*CompoundStmt:918*/
    *dest = (stencil_arg_t){.dst = (double*)((void*)0),
        .src = (double*)((void*)0),
        .n_rows = sourc->n_rows,
        .n_cols = sourc->n_cols,
        .n_tsteps = sourc->n_tsteps };
    dest->src = (double*)smalloc(sizeof(double) * (dest->n_rows) * (dest->n_cols));
    dest->dst = (double*)smalloc(sizeof(double) * (dest->n_rows) * (dest->n_cols));
    memcpy(dest->src, sourc->src, sizeof(double) * (dest->n_rows) * (dest->n_cols));
    memcpy(dest->dst, sourc->dst, sizeof(double) * (dest->n_rows) * (dest->n_cols));
}
/*Function: stencil_arg_dup, ID: 7*/
stencil_arg_t* stencil_arg_dup(stencil_arg_t* src)
{
    /*stencil_arg_dup:7*/
    /*CompoundStmt:976*/
    stencil_arg_t* dst = (stencil_arg_t*)smalloc(sizeof(stencil_arg_t));
    stencil_arg_copy(dst, src);
    return dst;
}
/*Function: run_kernel, ID: 8*/
void run_kernel(code_t code, stencil_t* stencil, const size_t n_reps, const char* label)
{
    /*run_kernel:8*/
    /*CompoundStmt:983*/
    loc_timer_t loc_timer = {.start = { 0 },
        .stop = { 0 },
        .cpu = 0L,
        .code = code,
        .data = (void*)stencil,
        .n_reps = n_reps };
    pthread_t timer_thd;
    spthread_create(&timer_thd, (pthread_attr_t*)((void*)0), timer, &loc_timer);
    if (0 != pthread_join(timer_thd, (void**)((void*)0)))
        fatal(1, "pthread_join");
    timespec_diff(&loc_timer.stop, &loc_timer.start, &loc_timer.stop);
    double t = loc_timer.stop.tv_sec + (double)loc_timer.stop.tv_nsec / 1.0E+9;
    printf("%s kernel time: %.5f secs\n", label, t / n_reps);
}
/*Function: main, ID: 9*/
int main(int argc, char* argv[])
{
    getOMPNumThreads();
    getOMPSchedulePolicy();
    getTPLoopThresholds();
    getNumTPs();
    affin = new ThreadAffinity(
        ompNumThreads / NUMTPS - 1, NUMTPS, COMPACT, getDARTSTPPolicy(), getDARTSMCPolicy());
    affinMaskRes = affin->generateMask();
    myDARTSRuntime = new Runtime(affin);
    RuntimeFinalCodelet = &(myDARTSRuntime->finalSignal);
    /*main:9*/
    /*CompoundStmt:1030*/
    uint64_t n_rows = 0UL, n_cols = 0UL, n_tm_steps = 1UL, n_reps = 1UL;
    switch (argc) {
    case 5:
        n_reps = strtoul(argv[4], (char**)((void*)0), 0);
    case 4:
        n_tm_steps = strtoul(argv[3], (char**)((void*)0), 0);
    case 3:
        n_rows = strtoul(argv[1], (char**)((void*)0), 0);
        n_cols = strtoul(argv[2], (char**)((void*)0), 0);
        break;
    case 2:
        n_rows = n_cols = strtoul(argv[1], (char**)((void*)0), 0);
        break;
    default:
        usage(argv[0]);
    }
    double *current_values = (double *)smalloc(sizeof(double) * n_rows * n_cols),
           *next_values = (double *)smalloc(sizeof(double) * n_rows * n_cols),
           *initial_matrix = (double *)smalloc(sizeof(double) * n_rows * n_cols);
    array2d_init(initial_matrix, n_rows, n_cols);
    uint64_t outer_start = 0, outer_stop = 0, inner_start = 0, inner_stop = 0;
    double outer_avg = 0., inner_avg = 0.;
    uint64_t* timings = (uint64_t*)smalloc(sizeof(uint64_t) * n_reps);
    char* str_n_threads = getenv("OMP_NUM_THREADS");
    long n_threads = str_n_threads ? strtoul(str_n_threads, (char**)((void*)0), 0) : 1L;
    outer_start = get_time();
    for (size_t i = 0; i < n_reps; ++i) {
        memcpy(current_values, initial_matrix, sizeof(double*) * n_rows * n_cols);
        memcpy(next_values, initial_matrix, sizeof(double*) * n_rows * n_cols);
        inner_start = get_time();
        fprintf(stderr, "rep %d\n", (int)i);
        stencil2D4pt_omp_v3(next_values, current_values, n_rows, n_cols, n_tm_steps);
        inner_stop = get_time() - inner_start;
        inner_avg += inner_stop;
        timings[i] = inner_stop;
    }
    outer_stop = get_time() - outer_start;
    outer_avg += outer_stop;
    inner_avg /= n_reps;
    outer_avg /= n_reps;
    printf("%lu * %lu, %lu, %lu, %ld, %lu, %lu, %-18.2f, %-18.2f\n",
           n_rows, n_cols, 0UL, n_threads, n_threads, n_tm_steps, n_reps, 
           inner_avg, outer_avg);

	double avg_time = 0;
	for(size_t i=0;i<n_reps;++i){
		// printf("%lu,",timings[i]);
		avg_time += timings[i];
	}
	avg_time /= n_reps;
	printf("avg_time = %E\n", avg_time);
    double* seq_matrix = (double*)smalloc(sizeof(double) * n_rows * n_cols);
    memcpy(seq_matrix, initial_matrix, sizeof(double*) * n_rows * n_cols);
    stencil2D4pt(seq_matrix, initial_matrix, n_rows, n_cols, n_tm_steps);
    if (result_is_correct(n_rows, n_cols, seq_matrix, next_values, initial_matrix))
        printf("success\n");
    sfree(seq_matrix);
    sfree(current_values);
    sfree(next_values);
    sfree(initial_matrix);
    sfree(timings);
    return 0;
}
