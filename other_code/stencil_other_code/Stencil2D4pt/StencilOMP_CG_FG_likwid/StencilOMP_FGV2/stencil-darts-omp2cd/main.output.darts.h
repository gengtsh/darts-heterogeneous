#ifndef _main_output_darts_h_
#define _main_output_darts_h_
#ifndef __DARTS_
#define __DARTS_
#endif
#include "./timer.h"
#include "TaskData.h"
#include "darts.h"
#include "ompTP.h"
#include "stencil.h"
#include "tbb/concurrent_vector.h"
#include "utils.h"
#include <limits.h>
#include <math.h>
#include <mutex>
#include <numa.h>
#include <stdbool.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
int main(int argc, char* argv[]);
void run_kernel(code_t code, stencil_t* stencil, const size_t n_reps, const char* label);
stencil_arg_t* stencil_arg_dup(stencil_arg_t* src);
void stencil_arg_copy(stencil_arg_t* dest, stencil_arg_t* sourc);
void stencil_destroy(stencil_t* stencil);
void stencil_init(stencil_t* stencil, stencil_funptr_t code, const size_t n_rows,
    const size_t n_cols, size_t n_tsteps);
void array2d_init(double* a, const size_t n_rows, const size_t n_cols);
bool result_is_correct(const size_t n_rows, const size_t n_cols, const double* restrict orig,
    const double* restrict res, const double* restrict init);
extern int DARTS_CODELETS_MULT;
extern int NUMTPS;
extern size_t numOfCUs;
extern darts::Codelet* RuntimeFinalCodelet;
extern darts::ThreadAffinity* affin;
extern bool affinMaskRes;
extern darts::Runtime* myDARTSRuntime;
extern std::vector<std::vector<void*> > threadFunctionStack;
extern size_t ompNumThreads;
extern int ompSchedulePolicy;
extern int ompScheduleChunk;
extern void omp_set_num_threads(unsigned long numThreadsToSet);
extern int omp_get_num_threads();
extern int omp_get_max_threads();
extern int omp_get_num_procs();
extern double omp_get_wtime();
extern void omp_init_lock(omp_lock_t* lock);
extern void omp_destroy_lock(omp_lock_t* lock);
extern void omp_set_lock(omp_lock_t* lock);
extern void omp_unset_lock(omp_lock_t* lock);
#endif
