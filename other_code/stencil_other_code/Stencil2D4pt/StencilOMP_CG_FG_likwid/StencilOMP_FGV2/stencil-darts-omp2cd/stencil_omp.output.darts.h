#ifndef _stencil_omp_output_darts_h_
#define _stencil_omp_output_darts_h_
#ifndef __DARTS_
#define __DARTS_
#endif
#include "TaskData.h"
#include "darts.h"
#include "ompTP.h"
#include "stencil.h"
#include "tbb/concurrent_vector.h"
#include "utils.h"
#include <limits.h>
#include <mutex>
#include <numa.h>
#include <pmmintrin.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
void stencil2D4pt_omp_v3(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps);
void stencil2D4pt_omp_v2(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps);
void* stencil_run(void* arg);
void stencil2D4pt_omp(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps);
void stencil2D4pt(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps);
class TP145;
class TP208;
class TP214;
class TP218;
/*Number of TPs to be used for the OMPFor in region TP218*/
#define NUMTPS218 NUMTPS
class TP313;
class TP315;
class TP317;
class TP319;
/*Class containing the inputs passed to task325*/
class _task325Inputs : public TaskData {
public:
    size_t chunk_darts325 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts325; /*OMP_SHARED*/
    int dCalc_outer325_size;
    int** dSwap_darts325; /*OMP_SHARED*/
    int dSwap_outer325_size;
    double* __restrict* dst_darts325; /*OMP_SHARED*/
    size_t n_cols_darts325 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts325; /*OMP_SHARED*/
    int* ts_darts325; /*OMP_SHARED*/
    double* DST_darts325; /*VARIABLE*/
    double* SRC_darts325; /*VARIABLE*/
    size_t i_darts325; /*VARIABLE*/
    size_t j_darts325; /*VARIABLE*/
    size_t pos1_darts325; /*VARIABLE*/
    size_t pos2_darts325; /*VARIABLE*/
    _task325Inputs()
        : TaskData()
    {
    }
    _task325Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts325(in_dCalc)
        , dCalc_outer325_size(in_dCalc_outer_size)
        , dSwap_darts325(in_dSwap)
        , dSwap_outer325_size(in_dSwap_outer_size)
        , dst_darts325(in_dst)
        , src_darts325(in_src)
        , ts_darts325(in_ts)
    {
        this->chunk_darts325 = *in_chunk;
        this->n_cols_darts325 = *in_n_cols;
    }
    ~_task325Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task325ListLock;
std::vector<std::list<_task325Inputs*> > task325List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task325Inputs*> task325ToDeleteVector;
class TP325;
/*Class containing the inputs passed to task371*/
class _task371Inputs : public TaskData {
public:
    size_t chunk_darts371 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts371; /*OMP_SHARED*/
    int dCalc_outer371_size;
    int** dSwap_darts371; /*OMP_SHARED*/
    int dSwap_outer371_size;
    double* __restrict* dst_darts371; /*OMP_SHARED*/
    size_t nThrdsMinus1_darts371 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts371 /*OMP_FIRSTPRIVATE*/;
    size_t n_rows_darts371 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts371; /*OMP_SHARED*/
    int* ts_darts371; /*OMP_SHARED*/
    double* DST_darts371; /*VARIABLE*/
    double* SRC_darts371; /*VARIABLE*/
    size_t i_darts371; /*VARIABLE*/
    size_t j_darts371; /*VARIABLE*/
    size_t pos1_darts371; /*VARIABLE*/
    size_t pos2_darts371; /*VARIABLE*/
    _task371Inputs()
        : TaskData()
    {
    }
    _task371Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_nThrdsMinus1,
        size_t* in_n_cols, size_t* in_n_rows, double* __restrict* in_src, int* in_ts,
        darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts371(in_dCalc)
        , dCalc_outer371_size(in_dCalc_outer_size)
        , dSwap_darts371(in_dSwap)
        , dSwap_outer371_size(in_dSwap_outer_size)
        , dst_darts371(in_dst)
        , src_darts371(in_src)
        , ts_darts371(in_ts)
    {
        this->chunk_darts371 = *in_chunk;
        this->nThrdsMinus1_darts371 = *in_nThrdsMinus1;
        this->n_cols_darts371 = *in_n_cols;
        this->n_rows_darts371 = *in_n_rows;
    }
    ~_task371Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task371ListLock;
std::vector<std::list<_task371Inputs*> > task371List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task371Inputs*> task371ToDeleteVector;
class TP371;
class TP420;
/*Class containing the inputs passed to task425*/
class _task425Inputs : public TaskData {
public:
    size_t chunk_darts425 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts425; /*OMP_SHARED*/
    int dCalc_outer425_size;
    int** dSwap_darts425; /*OMP_SHARED*/
    int dSwap_outer425_size;
    double* __restrict* dst_darts425; /*OMP_SHARED*/
    size_t k_darts425 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts425 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts425; /*OMP_SHARED*/
    int* ts_darts425; /*OMP_SHARED*/
    double* DST_darts425; /*VARIABLE*/
    double* SRC_darts425; /*VARIABLE*/
    size_t i_darts425; /*VARIABLE*/
    size_t j_darts425; /*VARIABLE*/
    size_t pos1_darts425; /*VARIABLE*/
    size_t pos2_darts425; /*VARIABLE*/
    _task425Inputs()
        : TaskData()
    {
    }
    _task425Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_k, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts425(in_dCalc)
        , dCalc_outer425_size(in_dCalc_outer_size)
        , dSwap_darts425(in_dSwap)
        , dSwap_outer425_size(in_dSwap_outer_size)
        , dst_darts425(in_dst)
        , src_darts425(in_src)
        , ts_darts425(in_ts)
    {
        this->chunk_darts425 = *in_chunk;
        this->k_darts425 = *in_k;
        this->n_cols_darts425 = *in_n_cols;
    }
    ~_task425Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task425ListLock;
std::vector<std::list<_task425Inputs*> > task425List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task425Inputs*> task425ToDeleteVector;
class TP425;
/*Class containing the inputs passed to task480*/
class _task480Inputs : public TaskData {
public:
    size_t chunk_darts480 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts480; /*OMP_SHARED*/
    int dCalc_outer480_size;
    int** dSwap_darts480; /*OMP_SHARED*/
    int dSwap_outer480_size;
    double* __restrict* dst_darts480; /*OMP_SHARED*/
    size_t n_cols_darts480 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts480; /*OMP_SHARED*/
    int* ts_darts480; /*OMP_SHARED*/
    double* DST_darts480; /*VARIABLE*/
    double* SRC_darts480; /*VARIABLE*/
    size_t i_darts480; /*VARIABLE*/
    size_t j_darts480; /*VARIABLE*/
    size_t pos1_darts480; /*VARIABLE*/
    size_t pos2_darts480; /*VARIABLE*/
    _task480Inputs()
        : TaskData()
    {
    }
    _task480Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts480(in_dCalc)
        , dCalc_outer480_size(in_dCalc_outer_size)
        , dSwap_darts480(in_dSwap)
        , dSwap_outer480_size(in_dSwap_outer_size)
        , dst_darts480(in_dst)
        , src_darts480(in_src)
        , ts_darts480(in_ts)
    {
        this->chunk_darts480 = *in_chunk;
        this->n_cols_darts480 = *in_n_cols;
    }
    ~_task480Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task480ListLock;
std::vector<std::list<_task480Inputs*> > task480List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task480Inputs*> task480ToDeleteVector;
class TP480;
/*Class containing the inputs passed to task524*/
class _task524Inputs : public TaskData {
public:
    size_t chunk_darts524 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts524; /*OMP_SHARED*/
    int dCalc_outer524_size;
    int** dSwap_darts524; /*OMP_SHARED*/
    int dSwap_outer524_size;
    double* __restrict* dst_darts524; /*OMP_SHARED*/
    size_t nThrdsMinus1_darts524 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts524 /*OMP_FIRSTPRIVATE*/;
    size_t n_rows_darts524 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts524; /*OMP_SHARED*/
    int* ts_darts524; /*OMP_SHARED*/
    double* DST_darts524; /*VARIABLE*/
    double* SRC_darts524; /*VARIABLE*/
    size_t i_darts524; /*VARIABLE*/
    size_t j_darts524; /*VARIABLE*/
    size_t pos1_darts524; /*VARIABLE*/
    size_t pos2_darts524; /*VARIABLE*/
    _task524Inputs()
        : TaskData()
    {
    }
    _task524Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_nThrdsMinus1,
        size_t* in_n_cols, size_t* in_n_rows, double* __restrict* in_src, int* in_ts,
        darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts524(in_dCalc)
        , dCalc_outer524_size(in_dCalc_outer_size)
        , dSwap_darts524(in_dSwap)
        , dSwap_outer524_size(in_dSwap_outer_size)
        , dst_darts524(in_dst)
        , src_darts524(in_src)
        , ts_darts524(in_ts)
    {
        this->chunk_darts524 = *in_chunk;
        this->nThrdsMinus1_darts524 = *in_nThrdsMinus1;
        this->n_cols_darts524 = *in_n_cols;
        this->n_rows_darts524 = *in_n_rows;
    }
    ~_task524Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task524ListLock;
std::vector<std::list<_task524Inputs*> > task524List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task524Inputs*> task524ToDeleteVector;
class TP524;
class TP571;
/*Class containing the inputs passed to task576*/
class _task576Inputs : public TaskData {
public:
    size_t chunk_darts576 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts576; /*OMP_SHARED*/
    int dCalc_outer576_size;
    int** dSwap_darts576; /*OMP_SHARED*/
    int dSwap_outer576_size;
    double* __restrict* dst_darts576; /*OMP_SHARED*/
    size_t k_darts576 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts576 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts576; /*OMP_SHARED*/
    int* ts_darts576; /*OMP_SHARED*/
    double* DST_darts576; /*VARIABLE*/
    double* SRC_darts576; /*VARIABLE*/
    size_t i_darts576; /*VARIABLE*/
    size_t j_darts576; /*VARIABLE*/
    size_t pos1_darts576; /*VARIABLE*/
    size_t pos2_darts576; /*VARIABLE*/
    _task576Inputs()
        : TaskData()
    {
    }
    _task576Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_k, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts576(in_dCalc)
        , dCalc_outer576_size(in_dCalc_outer_size)
        , dSwap_darts576(in_dSwap)
        , dSwap_outer576_size(in_dSwap_outer_size)
        , dst_darts576(in_dst)
        , src_darts576(in_src)
        , ts_darts576(in_ts)
    {
        this->chunk_darts576 = *in_chunk;
        this->k_darts576 = *in_k;
        this->n_cols_darts576 = *in_n_cols;
    }
    ~_task576Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task576ListLock;
std::vector<std::list<_task576Inputs*> > task576List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task576Inputs*> task576ToDeleteVector;
class TP576;
/*Class containing the inputs passed to task631*/
class _task631Inputs : public TaskData {
public:
    size_t chunk_darts631 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts631; /*OMP_SHARED*/
    int dCalc_outer631_size;
    int** dSwap_darts631; /*OMP_SHARED*/
    int dSwap_outer631_size;
    double* __restrict* dst_darts631; /*OMP_SHARED*/
    size_t n_cols_darts631 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts631; /*OMP_SHARED*/
    int* ts_darts631; /*OMP_SHARED*/
    double* DST_darts631; /*VARIABLE*/
    double* SRC_darts631; /*VARIABLE*/
    size_t i_darts631; /*VARIABLE*/
    size_t j_darts631; /*VARIABLE*/
    size_t pos1_darts631; /*VARIABLE*/
    size_t pos2_darts631; /*VARIABLE*/
    _task631Inputs()
        : TaskData()
    {
    }
    _task631Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts631(in_dCalc)
        , dCalc_outer631_size(in_dCalc_outer_size)
        , dSwap_darts631(in_dSwap)
        , dSwap_outer631_size(in_dSwap_outer_size)
        , dst_darts631(in_dst)
        , src_darts631(in_src)
        , ts_darts631(in_ts)
    {
        this->chunk_darts631 = *in_chunk;
        this->n_cols_darts631 = *in_n_cols;
    }
    ~_task631Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task631ListLock;
std::vector<std::list<_task631Inputs*> > task631List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task631Inputs*> task631ToDeleteVector;
class TP631;
/*Class containing the inputs passed to task677*/
class _task677Inputs : public TaskData {
public:
    size_t chunk_darts677 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts677; /*OMP_SHARED*/
    int dCalc_outer677_size;
    int** dSwap_darts677; /*OMP_SHARED*/
    int dSwap_outer677_size;
    double* __restrict* dst_darts677; /*OMP_SHARED*/
    size_t nThrdsMinus1_darts677 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts677 /*OMP_FIRSTPRIVATE*/;
    size_t n_rows_darts677 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts677; /*OMP_SHARED*/
    int* ts_darts677; /*OMP_SHARED*/
    double* DST_darts677; /*VARIABLE*/
    double* SRC_darts677; /*VARIABLE*/
    size_t i_darts677; /*VARIABLE*/
    size_t j_darts677; /*VARIABLE*/
    size_t pos1_darts677; /*VARIABLE*/
    size_t pos2_darts677; /*VARIABLE*/
    _task677Inputs()
        : TaskData()
    {
    }
    _task677Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_nThrdsMinus1,
        size_t* in_n_cols, size_t* in_n_rows, double* __restrict* in_src, int* in_ts,
        darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts677(in_dCalc)
        , dCalc_outer677_size(in_dCalc_outer_size)
        , dSwap_darts677(in_dSwap)
        , dSwap_outer677_size(in_dSwap_outer_size)
        , dst_darts677(in_dst)
        , src_darts677(in_src)
        , ts_darts677(in_ts)
    {
        this->chunk_darts677 = *in_chunk;
        this->nThrdsMinus1_darts677 = *in_nThrdsMinus1;
        this->n_cols_darts677 = *in_n_cols;
        this->n_rows_darts677 = *in_n_rows;
    }
    ~_task677Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task677ListLock;
std::vector<std::list<_task677Inputs*> > task677List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task677Inputs*> task677ToDeleteVector;
class TP677;
class TP726;
/*Class containing the inputs passed to task731*/
class _task731Inputs : public TaskData {
public:
    size_t chunk_darts731 /*OMP_FIRSTPRIVATE*/;
    int** dCalc_darts731; /*OMP_SHARED*/
    int dCalc_outer731_size;
    int** dSwap_darts731; /*OMP_SHARED*/
    int dSwap_outer731_size;
    double* __restrict* dst_darts731; /*OMP_SHARED*/
    size_t k_darts731 /*OMP_FIRSTPRIVATE*/;
    size_t n_cols_darts731 /*OMP_FIRSTPRIVATE*/;
    double* __restrict* src_darts731; /*OMP_SHARED*/
    int* ts_darts731; /*OMP_SHARED*/
    double* DST_darts731; /*VARIABLE*/
    double* SRC_darts731; /*VARIABLE*/
    size_t i_darts731; /*VARIABLE*/
    size_t j_darts731; /*VARIABLE*/
    size_t pos1_darts731; /*VARIABLE*/
    size_t pos2_darts731; /*VARIABLE*/
    _task731Inputs()
        : TaskData()
    {
    }
    _task731Inputs(size_t* in_chunk, int** in_dCalc, int in_dCalc_outer_size, int** in_dSwap,
        int in_dSwap_outer_size, double* __restrict* in_dst, size_t* in_k, size_t* in_n_cols,
        double* __restrict* in_src, int* in_ts, darts::Codelet* in_nextSyncCodelet)
        : TaskData(in_nextSyncCodelet)
        , dCalc_darts731(in_dCalc)
        , dCalc_outer731_size(in_dCalc_outer_size)
        , dSwap_darts731(in_dSwap)
        , dSwap_outer731_size(in_dSwap_outer_size)
        , dst_darts731(in_dst)
        , src_darts731(in_src)
        , ts_darts731(in_ts)
    {
        this->chunk_darts731 = *in_chunk;
        this->k_darts731 = *in_k;
        this->n_cols_darts731 = *in_n_cols;
    }
    ~_task731Inputs() {}
};
/*List for tasks created by each thread*/
std::mutex task731ListLock;
std::vector<std::list<_task731Inputs*> > task731List(MAXNUMTHREADS);
/*Vector containing all tasks created and that are pending to be deleted*/
tbb::concurrent_vector<_task731Inputs*> task731ToDeleteVector;
class TP731;
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
/*TP145: OMPParallelForDirective*/
class TP145 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets145 : public darts::Codelet {
    public:
        TP145* inputsTPParent;
        _barrierCodelets145()
            : darts::Codelet()
        {
        }
        _barrierCodelets145(uint32_t dep, uint32_t res, TP145* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    bool requestNewRangeIterations145(size_t* endRange, uint32_t codeletID);
    class _checkInCodelets146 : public darts::Codelet {
    public:
        TP145* myTP;
        TP145* inputsTPParent;
        size_t endRange;
        _checkInCodelets146()
            : darts::Codelet()
        {
        }
        _checkInCodelets146(uint32_t dep, uint32_t res, TP145* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP145* TPParent;
    TP145* controlTPParent;
    TP145* inputsTPParent;
    double** DST_darts145; /*OMP_SHARED - INPUT*/
    double** SRC_darts145; /*OMP_SHARED - INPUT*/
    size_t* n_cols_darts145 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_cols_outer145_ptr;
    size_t* n_rows_darts145 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_rows_outer145_ptr;
    size_t* n_tsteps_darts145 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_tsteps_outer145_ptr;
    size_t* i_darts145 /*VARIABLE*/;
    size_t initIteration145;
    size_t lastIteration145;
    size_t range145;
    size_t rangePerCodelet145;
    size_t minIteration145;
    size_t remainderRange145;
    _barrierCodelets145* barrierCodelets145;
    _checkInCodelets146* checkInCodelets146;
    TP145(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet,
        size_t in_initIteration, size_t in_lastIteration, double** in_DST, double** in_SRC,
        size_t* in_n_cols_outer145_ptr, size_t* in_n_rows_outer145_ptr,
        size_t* in_n_tsteps_outer145_ptr);
    ~TP145();
};
/*TP208: OMPParallelDirective*/
class TP208 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets208 : public darts::Codelet {
    public:
        TP208* inputsTPParent;
        _barrierCodelets208()
            : darts::Codelet()
        {
        }
        _barrierCodelets208(uint32_t dep, uint32_t res, TP208* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets210 : public darts::Codelet {
    public:
        TP208* myTP;
        TP208* inputsTPParent;
        _checkInCodelets210()
            : darts::Codelet()
        {
        }
        _checkInCodelets210(uint32_t dep, uint32_t res, TP208* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets214 : public darts::Codelet {
    public:
        TP208* myTP;
        TP208* inputsTPParent;
        _checkInCodelets214()
            : darts::Codelet()
        {
        }
        _checkInCodelets214(uint32_t dep, uint32_t res, TP208* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets817 : public darts::Codelet {
    public:
        TP208* myTP;
        TP208* inputsTPParent;
        _checkInCodelets817()
            : darts::Codelet()
        {
        }
        _checkInCodelets817(uint32_t dep, uint32_t res, TP208* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP208* TPParent;
    TP208* controlTPParent;
    TP208* inputsTPParent;
    double* __restrict* dst_darts208; /*OMP_SHARED - INPUT*/
    size_t* n_cols_darts208 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_cols_outer208_ptr;
    size_t* n_rows_darts208 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_rows_outer208_ptr;
    size_t* n_tsteps_darts208 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_tsteps_outer208_ptr;
    double* __restrict* src_darts208; /*OMP_SHARED - INPUT*/
    double** DST_darts208 /*VARIABLE*/;
    double** SRC_darts208 /*VARIABLE*/;
    size_t* n_ts_darts208 /*VARIABLE*/;
    unsigned int TP214_LoopCounter;
    unsigned int* TP214_LoopCounterPerThread;
    tbb::concurrent_vector<TP214*> TP214PtrVec;
    _barrierCodelets208* barrierCodelets208;
    _checkInCodelets210* checkInCodelets210;
    _checkInCodelets214* checkInCodelets214;
    _checkInCodelets817* checkInCodelets817;
    TP208(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet,
        double* __restrict* in_dst, size_t* in_n_cols_outer208_ptr, size_t* in_n_rows_outer208_ptr,
        size_t* in_n_tsteps_outer208_ptr, double* __restrict* in_src);
    ~TP208();
};
/*TP214: WhileStmt*/
class TP214 : public ompTP {
public:
    class _checkInCodelets218 : public darts::Codelet {
    public:
        TP214* myTP;
        TP208* inputsTPParent;
        _checkInCodelets218()
            : darts::Codelet()
        {
        }
        _checkInCodelets218(uint32_t dep, uint32_t res, TP214* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets274 : public darts::Codelet {
    public:
        TP214* myTP;
        TP208* inputsTPParent;
        _checkInCodelets274()
            : darts::Codelet()
        {
        }
        _checkInCodelets274(uint32_t dep, uint32_t res, TP214* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets279 : public darts::Codelet {
    public:
        TP208* inputsTPParent;
        _barrierCodelets279()
            : darts::Codelet()
        {
        }
        _barrierCodelets279(uint32_t dep, uint32_t res, TP214* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP208* TPParent;
    TP214* controlTPParent;
    TP208* inputsTPParent;
    TP214** ptrToThisTP;
    TP218** TP218Ptr;
    size_t* TP218_alreadyLaunched;
    int numTPsSet218;
    int numTPsReady218;
    size_t TPsToUse218;
    size_t codeletsPerTP218;
    size_t totalCodelets218;
    _checkInCodelets218* checkInCodelets218;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets218* firstCodelet;
#endif
    _checkInCodelets274* checkInCodelets274;
    _barrierCodelets279* barrierCodelets279;
    TP214(int in_numThreads, int in_mainCodeletID, TP208* in_TPParent, TP208* in_inputsTPParent,
        TP214** in_ptrToThisTP);
    ~TP214();
};
/*TP218: OMPForDirective*/
class TP218 : public ompTP {
public:
    bool requestNewRangeIterations218(size_t* endRange, uint32_t codeletID);
    class _checkInCodelets219 : public darts::Codelet {
    public:
        TP218* myTP;
        TP218* inputsTPParent;
        size_t endRange;
        _checkInCodelets219()
            : darts::Codelet()
        {
        }
        _checkInCodelets219(uint32_t dep, uint32_t res, TP218* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP214* TPParent;
    TP218* controlTPParent;
    TP218* inputsTPParent;
    double*** DST_darts218 /*OMP_SHARED_PRIVATE - INPUT*/;
    double*** SRC_darts218 /*OMP_SHARED_PRIVATE - INPUT*/;
    size_t* n_cols_darts218 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_cols_outer218_ptr;
    size_t* n_rows_darts218 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_rows_outer218_ptr;
    size_t* i_darts218 /*VARIABLE*/;
    size_t initIteration218;
    size_t lastIteration218;
    size_t range218;
    size_t rangePerCodelet218;
    size_t minIteration218;
    size_t remainderRange218;
    size_t readyCodelets;
    int baseNumThreads;
    int* signalNextReady;
    _checkInCodelets219* checkInCodelets219;
#if USE_SPIN_CODELETS == 0
    _checkInCodelets219* firstCodelet;
#endif
    TP218(int in_numThreads, int in_mainCodeletID, TP214* in_TPParent, size_t in_initIteration,
        size_t in_lastIteration, double** in_DST, double** in_SRC, size_t* in_n_cols_outer218_ptr,
        size_t* in_n_rows_outer218_ptr, TP218** in_ptrToThisTP);
    void inline dispatchCodelet(size_t codeletID);
    ~TP218();
};
/*TP313: OMPParallelDirective*/
class TP313 : public darts::ThreadedProcedure {
public:
    class _barrierCodelets313 : public darts::Codelet {
    public:
        TP313* inputsTPParent;
        _barrierCodelets313()
            : darts::Codelet()
        {
        }
        _barrierCodelets313(uint32_t dep, uint32_t res, TP313* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets315 : public darts::Codelet {
    public:
        TP313* myTP;
        TP313* inputsTPParent;
        _checkInCodelets315()
            : darts::Codelet()
        {
        }
        _checkInCodelets315(uint32_t dep, uint32_t res, TP313* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets315 : public darts::Codelet {
    public:
        TP313* inputsTPParent;
        _barrierCodelets315()
            : darts::Codelet()
        {
        }
        _barrierCodelets315(uint32_t dep, uint32_t res, TP313* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    darts::Codelet* nextCodelet;
    TP313* TPParent;
    TP313* controlTPParent;
    TP313* inputsTPParent;
    size_t* chunk_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* chunk_outer313_ptr;
    int** dCalc_darts313; /*OMP_SHARED - INPUT*/
    uint64_t dCalc_outer313_size;
    int** dSwap_darts313; /*OMP_SHARED - INPUT*/
    uint64_t dSwap_outer313_size;
    double* __restrict* dst_darts313; /*OMP_SHARED - INPUT*/
    size_t* nThrds_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* nThrds_outer313_ptr;
    size_t* nThrdsMinus1_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* nThrdsMinus1_outer313_ptr;
    size_t* n_cols_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_cols_outer313_ptr;
    size_t* n_rows_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_rows_outer313_ptr;
    size_t* n_tsteps_darts313 /*OMP_FIRSTPRIVATE - INPUT*/;
    size_t* n_tsteps_outer313_ptr;
    double* __restrict* src_darts313; /*OMP_SHARED - INPUT*/
    int* ts_darts313; /*OMP_SHARED - INPUT*/
    TP315* TP315Ptr;
    size_t TP315_alreadyLaunched;
    _barrierCodelets313* barrierCodelets313;
    _checkInCodelets315* checkInCodelets315;
    _barrierCodelets315* barrierCodelets315;
    TP313(int in_numThreads, int in_mainCodeletID, darts::Codelet* in_nextCodelet,
        size_t* in_chunk_outer313_ptr, int** in_dCalc, int in_dCalc_outer313_size, int** in_dSwap,
        int in_dSwap_outer313_size, double* __restrict* in_dst, size_t* in_nThrds_outer313_ptr,
        size_t* in_nThrdsMinus1_outer313_ptr, size_t* in_n_cols_outer313_ptr,
        size_t* in_n_rows_outer313_ptr, size_t* in_n_tsteps_outer313_ptr,
        double* __restrict* in_src, int* in_ts);
    ~TP313();
};
/*TP315: OMPSingleDirective*/
class TP315 : public ompOMPSingleDirectiveTP {
public:
    class _checkInCodelets317 : public darts::Codelet {
    public:
        TP315* myTP;
        TP315* inputsTPParent;
        _checkInCodelets317()
            : darts::Codelet()
        {
        }
        _checkInCodelets317(uint32_t dep, uint32_t res, TP315* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _barrierCodelets317 : public darts::Codelet {
    public:
        TP315* inputsTPParent;
        _barrierCodelets317()
            : darts::Codelet()
        {
        }
        _barrierCodelets317(uint32_t dep, uint32_t res, TP315* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP313* TPParent;
    TP315* controlTPParent;
    TP315* inputsTPParent;
    size_t chunk_darts315; /*SINGLE_THREADED - INPUT*/
    int** dCalc_darts315; /*OMP_SHARED - INPUT*/
    uint64_t dCalc_outer315_size;
    int** dSwap_darts315; /*OMP_SHARED - INPUT*/
    uint64_t dSwap_outer315_size;
    double* __restrict* dst_darts315; /*OMP_SHARED - INPUT*/
    size_t nThrdsMinus1_darts315; /*SINGLE_THREADED - INPUT*/
    size_t n_cols_darts315; /*SINGLE_THREADED - INPUT*/
    size_t n_rows_darts315; /*SINGLE_THREADED - INPUT*/
    size_t n_tsteps_darts315; /*SINGLE_THREADED - INPUT*/
    double* __restrict* src_darts315; /*OMP_SHARED - INPUT*/
    int* ts_darts315; /*OMP_SHARED - INPUT*/
    _checkInCodelets317 checkInCodelets317;
    _barrierCodelets317* barrierCodelets317;
    TP315(int in_numThreads, int in_mainCodeletID, TP313* in_TPParent,
        size_t* in_chunk_outer315_ptr, int** in_dCalc, int in_dCalc_outer315_size, int** in_dSwap,
        int in_dSwap_outer315_size, double* __restrict* in_dst,
        size_t* in_nThrdsMinus1_outer315_ptr, size_t* in_n_cols_outer315_ptr,
        size_t* in_n_rows_outer315_ptr, size_t* in_n_tsteps_outer315_ptr,
        double* __restrict* in_src, int* in_ts);
    ~TP315();
};
/*TP317: OMPTaskgroupDirective*/
class TP317 : public ompOMPTaskgroupDirectiveTP {
public:
    class _checkInCodelets320 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets320()
            : darts::Codelet()
        {
        }
        _checkInCodelets320(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets319 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets319()
            : darts::Codelet()
        {
        }
        _checkInCodelets319(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets322 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets322()
            : darts::Codelet()
        {
        }
        _checkInCodelets322(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets628 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets628()
            : darts::Codelet()
        {
        }
        _checkInCodelets628(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets631 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets631()
            : darts::Codelet()
        {
        }
        _checkInCodelets631(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets677 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets677()
            : darts::Codelet()
        {
        }
        _checkInCodelets677(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets727 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets727()
            : darts::Codelet()
        {
        }
        _checkInCodelets727(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets726 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets726()
            : darts::Codelet()
        {
        }
        _checkInCodelets726(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets729 : public darts::Codelet {
    public:
        TP317* myTP;
        TP317* inputsTPParent;
        _checkInCodelets729()
            : darts::Codelet()
        {
        }
        _checkInCodelets729(uint32_t dep, uint32_t res, TP317* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP315* TPParent;
    TP317* controlTPParent;
    TP317* inputsTPParent;
    size_t chunk_darts317; /*SINGLE_THREADED - INPUT*/
    int** dCalc_darts317; /*OMP_SHARED - INPUT*/
    uint64_t dCalc_outer317_size;
    int** dSwap_darts317; /*OMP_SHARED - INPUT*/
    uint64_t dSwap_outer317_size;
    double* __restrict* dst_darts317; /*OMP_SHARED - INPUT*/
    size_t nThrdsMinus1_darts317; /*SINGLE_THREADED - INPUT*/
    size_t n_cols_darts317; /*SINGLE_THREADED - INPUT*/
    size_t n_rows_darts317; /*SINGLE_THREADED - INPUT*/
    size_t n_tsteps_darts317; /*SINGLE_THREADED - INPUT*/
    double* __restrict* src_darts317; /*OMP_SHARED - INPUT*/
    int ts_darts317; /*SINGLE_THREADED - INPUT*/
    size_t k_darts317 /*VARIABLE*/;
    unsigned int TP319_LoopCounter;
    unsigned int* TP319_LoopCounterPerThread;
    tbb::concurrent_vector<TP319*> TP319PtrVec;
    unsigned int TP726_LoopCounter;
    unsigned int* TP726_LoopCounterPerThread;
    tbb::concurrent_vector<TP726*> TP726PtrVec;
    _checkInCodelets320 checkInCodelets320;
    _checkInCodelets319 checkInCodelets319;
    _checkInCodelets322 checkInCodelets322;
    _checkInCodelets628 checkInCodelets628;
    _checkInCodelets631 checkInCodelets631;
    _checkInCodelets677 checkInCodelets677;
    _checkInCodelets727 checkInCodelets727;
    _checkInCodelets726 checkInCodelets726;
    _checkInCodelets729 checkInCodelets729;
    TP317(int in_numThreads, int in_mainCodeletID, TP315* in_TPParent,
        size_t* in_chunk_outer317_ptr, int** in_dCalc, int in_dCalc_outer317_size, int** in_dSwap,
        int in_dSwap_outer317_size, double* __restrict* in_dst,
        size_t* in_nThrdsMinus1_outer317_ptr, size_t* in_n_cols_outer317_ptr,
        size_t* in_n_rows_outer317_ptr, size_t* in_n_tsteps_outer317_ptr,
        double* __restrict* in_src);
    ~TP317();
};
/*TP319: ForStmt*/
class TP319 : public ompTP {
public:
    class _checkInCodelets325 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets325()
            : darts::Codelet()
        {
        }
        _checkInCodelets325(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets371 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets371()
            : darts::Codelet()
        {
        }
        _checkInCodelets371(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets421 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets421()
            : darts::Codelet()
        {
        }
        _checkInCodelets421(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets420 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets420()
            : darts::Codelet()
        {
        }
        _checkInCodelets420(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets423 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets423()
            : darts::Codelet()
        {
        }
        _checkInCodelets423(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets480 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets480()
            : darts::Codelet()
        {
        }
        _checkInCodelets480(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets524 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets524()
            : darts::Codelet()
        {
        }
        _checkInCodelets524(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets572 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets572()
            : darts::Codelet()
        {
        }
        _checkInCodelets572(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets571 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets571()
            : darts::Codelet()
        {
        }
        _checkInCodelets571(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    class _checkInCodelets574 : public darts::Codelet {
    public:
        TP319* myTP;
        TP317* inputsTPParent;
        _checkInCodelets574()
            : darts::Codelet()
        {
        }
        _checkInCodelets574(uint32_t dep, uint32_t res, TP319* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP317* TPParent;
    TP319* controlTPParent;
    TP317* inputsTPParent;
    TP319** ptrToThisTP;
    unsigned int TP420_LoopCounter;
    unsigned int* TP420_LoopCounterPerThread;
    tbb::concurrent_vector<TP420*> TP420PtrVec;
    unsigned int TP571_LoopCounter;
    unsigned int* TP571_LoopCounterPerThread;
    tbb::concurrent_vector<TP571*> TP571PtrVec;
    _checkInCodelets325 checkInCodelets325;
    _checkInCodelets371 checkInCodelets371;
    _checkInCodelets421 checkInCodelets421;
    _checkInCodelets420 checkInCodelets420;
    _checkInCodelets423 checkInCodelets423;
    _checkInCodelets480 checkInCodelets480;
    _checkInCodelets524 checkInCodelets524;
    _checkInCodelets572 checkInCodelets572;
    _checkInCodelets571 checkInCodelets571;
    _checkInCodelets574 checkInCodelets574;
    TP319(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent, TP317* in_inputsTPParent,
        TP319** in_ptrToThisTP);
    ~TP319();
};
/*TP325: OMPTaskDirective*/
class TP325 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets326 : public darts::Codelet {
    public:
        TP325* myTP;
        TP325* inputsTPParent;
        _task325Inputs* taskInputs;
        _checkInCodelets326()
            : darts::Codelet()
        {
        }
        _checkInCodelets326(uint32_t dep, uint32_t res, TP325* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task325Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets335 : public darts::Codelet {
    public:
        TP325* myTP;
        TP325* inputsTPParent;
        _task325Inputs* taskInputs;
        _checkInCodelets335()
            : darts::Codelet()
        {
        }
        _checkInCodelets335(uint32_t dep, uint32_t res, TP325* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task325Inputs)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP325* controlTPParent;
    TP325* inputsTPParent;
    _task325Inputs* task325Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task325Inputs*>* parentTask325List;
    std::mutex* parentTask325ListLock;
    std::list<_task325Inputs*>::iterator parentTask325It;
    double* DST_darts325 /*VARIABLE*/;
    double* SRC_darts325 /*VARIABLE*/;
    size_t i_darts325 /*VARIABLE*/;
    size_t j_darts325 /*VARIABLE*/;
    size_t pos1_darts325 /*VARIABLE*/;
    size_t pos2_darts325 /*VARIABLE*/;
    _checkInCodelets326 checkInCodelets326;
    _checkInCodelets335 checkInCodelets335;
    TP325(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
        _task325Inputs* in_task325Inputs, std::list<_task325Inputs*>* in_parentTask325List,
        std::mutex* in_parentTask325ListLock,
        std::list<_task325Inputs*>::iterator in_parentTask325It);
    ~TP325();
};
/*TP371: OMPTaskDirective*/
class TP371 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets372 : public darts::Codelet {
    public:
        TP371* myTP;
        TP371* inputsTPParent;
        _task371Inputs* taskInputs;
        _checkInCodelets372()
            : darts::Codelet()
        {
        }
        _checkInCodelets372(uint32_t dep, uint32_t res, TP371* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task371Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets382 : public darts::Codelet {
    public:
        TP371* myTP;
        TP371* inputsTPParent;
        _task371Inputs* taskInputs;
        _checkInCodelets382()
            : darts::Codelet()
        {
        }
        _checkInCodelets382(uint32_t dep, uint32_t res, TP371* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task371Inputs)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP371* controlTPParent;
    TP371* inputsTPParent;
    _task371Inputs* task371Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task371Inputs*>* parentTask371List;
    std::mutex* parentTask371ListLock;
    std::list<_task371Inputs*>::iterator parentTask371It;
    double* DST_darts371 /*VARIABLE*/;
    double* SRC_darts371 /*VARIABLE*/;
    size_t i_darts371 /*VARIABLE*/;
    size_t j_darts371 /*VARIABLE*/;
    size_t pos1_darts371 /*VARIABLE*/;
    size_t pos2_darts371 /*VARIABLE*/;
    _checkInCodelets372 checkInCodelets372;
    _checkInCodelets382 checkInCodelets382;
    TP371(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
        _task371Inputs* in_task371Inputs, std::list<_task371Inputs*>* in_parentTask371List,
        std::mutex* in_parentTask371ListLock,
        std::list<_task371Inputs*>::iterator in_parentTask371It);
    ~TP371();
};
/*TP420: ForStmt*/
class TP420 : public ompTP {
public:
    class _checkInCodelets425 : public darts::Codelet {
    public:
        TP420* myTP;
        TP317* inputsTPParent;
        _checkInCodelets425()
            : darts::Codelet()
        {
        }
        _checkInCodelets425(uint32_t dep, uint32_t res, TP420* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP420* controlTPParent;
    TP317* inputsTPParent;
    TP420** ptrToThisTP;
    _checkInCodelets425 checkInCodelets425;
    TP420(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent, TP317* in_inputsTPParent,
        TP420** in_ptrToThisTP);
    ~TP420();
};
/*TP425: OMPTaskDirective*/
class TP425 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets426 : public darts::Codelet {
    public:
        TP425* myTP;
        TP425* inputsTPParent;
        _task425Inputs* taskInputs;
        _checkInCodelets426()
            : darts::Codelet()
        {
        }
        _checkInCodelets426(uint32_t dep, uint32_t res, TP425* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task425Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets440 : public darts::Codelet {
    public:
        TP425* myTP;
        TP425* inputsTPParent;
        _task425Inputs* taskInputs;
        _checkInCodelets440()
            : darts::Codelet()
        {
        }
        _checkInCodelets440(uint32_t dep, uint32_t res, TP425* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task425Inputs)
        {
        }
        void fire(void);
    };
    TP420* TPParent;
    TP425* controlTPParent;
    TP425* inputsTPParent;
    _task425Inputs* task425Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task425Inputs*>* parentTask425List;
    std::mutex* parentTask425ListLock;
    std::list<_task425Inputs*>::iterator parentTask425It;
    double* DST_darts425 /*VARIABLE*/;
    double* SRC_darts425 /*VARIABLE*/;
    size_t i_darts425 /*VARIABLE*/;
    size_t j_darts425 /*VARIABLE*/;
    size_t pos1_darts425 /*VARIABLE*/;
    size_t pos2_darts425 /*VARIABLE*/;
    _checkInCodelets426 checkInCodelets426;
    _checkInCodelets440 checkInCodelets440;
    TP425(int in_numThreads, int in_mainCodeletID, TP420* in_TPParent,
        _task425Inputs* in_task425Inputs, std::list<_task425Inputs*>* in_parentTask425List,
        std::mutex* in_parentTask425ListLock,
        std::list<_task425Inputs*>::iterator in_parentTask425It);
    ~TP425();
};
/*TP480: OMPTaskDirective*/
class TP480 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets481 : public darts::Codelet {
    public:
        TP480* myTP;
        TP480* inputsTPParent;
        _task480Inputs* taskInputs;
        _checkInCodelets481()
            : darts::Codelet()
        {
        }
        _checkInCodelets481(uint32_t dep, uint32_t res, TP480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task480Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets488 : public darts::Codelet {
    public:
        TP480* myTP;
        TP480* inputsTPParent;
        _task480Inputs* taskInputs;
        _checkInCodelets488()
            : darts::Codelet()
        {
        }
        _checkInCodelets488(uint32_t dep, uint32_t res, TP480* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task480Inputs)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP480* controlTPParent;
    TP480* inputsTPParent;
    _task480Inputs* task480Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task480Inputs*>* parentTask480List;
    std::mutex* parentTask480ListLock;
    std::list<_task480Inputs*>::iterator parentTask480It;
    double* DST_darts480 /*VARIABLE*/;
    double* SRC_darts480 /*VARIABLE*/;
    size_t i_darts480 /*VARIABLE*/;
    size_t j_darts480 /*VARIABLE*/;
    size_t pos1_darts480 /*VARIABLE*/;
    size_t pos2_darts480 /*VARIABLE*/;
    _checkInCodelets481 checkInCodelets481;
    _checkInCodelets488 checkInCodelets488;
    TP480(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
        _task480Inputs* in_task480Inputs, std::list<_task480Inputs*>* in_parentTask480List,
        std::mutex* in_parentTask480ListLock,
        std::list<_task480Inputs*>::iterator in_parentTask480It);
    ~TP480();
};
/*TP524: OMPTaskDirective*/
class TP524 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets525 : public darts::Codelet {
    public:
        TP524* myTP;
        TP524* inputsTPParent;
        _task524Inputs* taskInputs;
        _checkInCodelets525()
            : darts::Codelet()
        {
        }
        _checkInCodelets525(uint32_t dep, uint32_t res, TP524* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task524Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets533 : public darts::Codelet {
    public:
        TP524* myTP;
        TP524* inputsTPParent;
        _task524Inputs* taskInputs;
        _checkInCodelets533()
            : darts::Codelet()
        {
        }
        _checkInCodelets533(uint32_t dep, uint32_t res, TP524* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task524Inputs)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP524* controlTPParent;
    TP524* inputsTPParent;
    _task524Inputs* task524Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task524Inputs*>* parentTask524List;
    std::mutex* parentTask524ListLock;
    std::list<_task524Inputs*>::iterator parentTask524It;
    double* DST_darts524 /*VARIABLE*/;
    double* SRC_darts524 /*VARIABLE*/;
    size_t i_darts524 /*VARIABLE*/;
    size_t j_darts524 /*VARIABLE*/;
    size_t pos1_darts524 /*VARIABLE*/;
    size_t pos2_darts524 /*VARIABLE*/;
    _checkInCodelets525 checkInCodelets525;
    _checkInCodelets533 checkInCodelets533;
    TP524(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
        _task524Inputs* in_task524Inputs, std::list<_task524Inputs*>* in_parentTask524List,
        std::mutex* in_parentTask524ListLock,
        std::list<_task524Inputs*>::iterator in_parentTask524It);
    ~TP524();
};
/*TP571: ForStmt*/
class TP571 : public ompTP {
public:
    class _checkInCodelets576 : public darts::Codelet {
    public:
        TP571* myTP;
        TP317* inputsTPParent;
        _checkInCodelets576()
            : darts::Codelet()
        {
        }
        _checkInCodelets576(uint32_t dep, uint32_t res, TP571* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP319* TPParent;
    TP571* controlTPParent;
    TP317* inputsTPParent;
    TP571** ptrToThisTP;
    _checkInCodelets576 checkInCodelets576;
    TP571(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent, TP317* in_inputsTPParent,
        TP571** in_ptrToThisTP);
    ~TP571();
};
/*TP576: OMPTaskDirective*/
class TP576 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets577 : public darts::Codelet {
    public:
        TP576* myTP;
        TP576* inputsTPParent;
        _task576Inputs* taskInputs;
        _checkInCodelets577()
            : darts::Codelet()
        {
        }
        _checkInCodelets577(uint32_t dep, uint32_t res, TP576* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task576Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets588 : public darts::Codelet {
    public:
        TP576* myTP;
        TP576* inputsTPParent;
        _task576Inputs* taskInputs;
        _checkInCodelets588()
            : darts::Codelet()
        {
        }
        _checkInCodelets588(uint32_t dep, uint32_t res, TP576* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task576Inputs)
        {
        }
        void fire(void);
    };
    TP571* TPParent;
    TP576* controlTPParent;
    TP576* inputsTPParent;
    _task576Inputs* task576Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task576Inputs*>* parentTask576List;
    std::mutex* parentTask576ListLock;
    std::list<_task576Inputs*>::iterator parentTask576It;
    double* DST_darts576 /*VARIABLE*/;
    double* SRC_darts576 /*VARIABLE*/;
    size_t i_darts576 /*VARIABLE*/;
    size_t j_darts576 /*VARIABLE*/;
    size_t pos1_darts576 /*VARIABLE*/;
    size_t pos2_darts576 /*VARIABLE*/;
    _checkInCodelets577 checkInCodelets577;
    _checkInCodelets588 checkInCodelets588;
    TP576(int in_numThreads, int in_mainCodeletID, TP571* in_TPParent,
        _task576Inputs* in_task576Inputs, std::list<_task576Inputs*>* in_parentTask576List,
        std::mutex* in_parentTask576ListLock,
        std::list<_task576Inputs*>::iterator in_parentTask576It);
    ~TP576();
};
/*TP631: OMPTaskDirective*/
class TP631 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets632 : public darts::Codelet {
    public:
        TP631* myTP;
        TP631* inputsTPParent;
        _task631Inputs* taskInputs;
        _checkInCodelets632()
            : darts::Codelet()
        {
        }
        _checkInCodelets632(uint32_t dep, uint32_t res, TP631* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task631Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets641 : public darts::Codelet {
    public:
        TP631* myTP;
        TP631* inputsTPParent;
        _task631Inputs* taskInputs;
        _checkInCodelets641()
            : darts::Codelet()
        {
        }
        _checkInCodelets641(uint32_t dep, uint32_t res, TP631* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task631Inputs)
        {
        }
        void fire(void);
    };
    TP317* TPParent;
    TP631* controlTPParent;
    TP631* inputsTPParent;
    _task631Inputs* task631Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task631Inputs*>* parentTask631List;
    std::mutex* parentTask631ListLock;
    std::list<_task631Inputs*>::iterator parentTask631It;
    double* DST_darts631 /*VARIABLE*/;
    double* SRC_darts631 /*VARIABLE*/;
    size_t i_darts631 /*VARIABLE*/;
    size_t j_darts631 /*VARIABLE*/;
    size_t pos1_darts631 /*VARIABLE*/;
    size_t pos2_darts631 /*VARIABLE*/;
    _checkInCodelets632 checkInCodelets632;
    _checkInCodelets641 checkInCodelets641;
    TP631(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent,
        _task631Inputs* in_task631Inputs, std::list<_task631Inputs*>* in_parentTask631List,
        std::mutex* in_parentTask631ListLock,
        std::list<_task631Inputs*>::iterator in_parentTask631It);
    ~TP631();
};
/*TP677: OMPTaskDirective*/
class TP677 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets678 : public darts::Codelet {
    public:
        TP677* myTP;
        TP677* inputsTPParent;
        _task677Inputs* taskInputs;
        _checkInCodelets678()
            : darts::Codelet()
        {
        }
        _checkInCodelets678(uint32_t dep, uint32_t res, TP677* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task677Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets688 : public darts::Codelet {
    public:
        TP677* myTP;
        TP677* inputsTPParent;
        _task677Inputs* taskInputs;
        _checkInCodelets688()
            : darts::Codelet()
        {
        }
        _checkInCodelets688(uint32_t dep, uint32_t res, TP677* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task677Inputs)
        {
        }
        void fire(void);
    };
    TP317* TPParent;
    TP677* controlTPParent;
    TP677* inputsTPParent;
    _task677Inputs* task677Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task677Inputs*>* parentTask677List;
    std::mutex* parentTask677ListLock;
    std::list<_task677Inputs*>::iterator parentTask677It;
    double* DST_darts677 /*VARIABLE*/;
    double* SRC_darts677 /*VARIABLE*/;
    size_t i_darts677 /*VARIABLE*/;
    size_t j_darts677 /*VARIABLE*/;
    size_t pos1_darts677 /*VARIABLE*/;
    size_t pos2_darts677 /*VARIABLE*/;
    _checkInCodelets678 checkInCodelets678;
    _checkInCodelets688 checkInCodelets688;
    TP677(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent,
        _task677Inputs* in_task677Inputs, std::list<_task677Inputs*>* in_parentTask677List,
        std::mutex* in_parentTask677ListLock,
        std::list<_task677Inputs*>::iterator in_parentTask677It);
    ~TP677();
};
/*TP726: ForStmt*/
class TP726 : public ompTP {
public:
    class _checkInCodelets731 : public darts::Codelet {
    public:
        TP726* myTP;
        TP317* inputsTPParent;
        _checkInCodelets731()
            : darts::Codelet()
        {
        }
        _checkInCodelets731(uint32_t dep, uint32_t res, TP726* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
        {
        }
        void fire(void);
    };
    TP317* TPParent;
    TP726* controlTPParent;
    TP317* inputsTPParent;
    TP726** ptrToThisTP;
    _checkInCodelets731 checkInCodelets731;
    TP726(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent, TP317* in_inputsTPParent,
        TP726** in_ptrToThisTP);
    ~TP726();
};
/*TP731: OMPTaskDirective*/
class TP731 : public ompOMPTaskDirectiveTP {
public:
    class _checkInCodelets732 : public darts::Codelet {
    public:
        TP731* myTP;
        TP731* inputsTPParent;
        _task731Inputs* taskInputs;
        _checkInCodelets732()
            : darts::Codelet()
        {
        }
        _checkInCodelets732(uint32_t dep, uint32_t res, TP731* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task731Inputs)
        {
        }
        void fire(void);
    };
    class _checkInCodelets746 : public darts::Codelet {
    public:
        TP731* myTP;
        TP731* inputsTPParent;
        _task731Inputs* taskInputs;
        _checkInCodelets746()
            : darts::Codelet()
        {
        }
        _checkInCodelets746(uint32_t dep, uint32_t res, TP731* myTP, uint32_t id)
            : darts::Codelet(dep, res, myTP, LONGWAIT, id)
            , myTP(myTP)
            , inputsTPParent(myTP->inputsTPParent)
            , taskInputs(inputsTPParent->task731Inputs)
        {
        }
        void fire(void);
    };
    TP726* TPParent;
    TP731* controlTPParent;
    TP731* inputsTPParent;
    _task731Inputs* task731Inputs;
    /*Pointer to the list storing this task data*/
    std::list<_task731Inputs*>* parentTask731List;
    std::mutex* parentTask731ListLock;
    std::list<_task731Inputs*>::iterator parentTask731It;
    double* DST_darts731 /*VARIABLE*/;
    double* SRC_darts731 /*VARIABLE*/;
    size_t i_darts731 /*VARIABLE*/;
    size_t j_darts731 /*VARIABLE*/;
    size_t pos1_darts731 /*VARIABLE*/;
    size_t pos2_darts731 /*VARIABLE*/;
    _checkInCodelets732 checkInCodelets732;
    _checkInCodelets746 checkInCodelets746;
    TP731(int in_numThreads, int in_mainCodeletID, TP726* in_TPParent,
        _task731Inputs* in_task731Inputs, std::list<_task731Inputs*>* in_parentTask731List,
        std::mutex* in_parentTask731ListLock,
        std::list<_task731Inputs*>::iterator in_parentTask731It);
    ~TP731();
};
#endif
