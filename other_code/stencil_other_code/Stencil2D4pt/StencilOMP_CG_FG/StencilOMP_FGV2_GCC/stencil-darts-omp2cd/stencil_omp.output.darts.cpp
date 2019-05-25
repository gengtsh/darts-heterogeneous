#include "stencil_omp.output.darts.h"
using namespace darts;
using namespace std;
/*Function: stencil2D4pt, ID: 11*/
void stencil2D4pt(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps)
{
    /*stencil2D4pt:11*/
    /*CompoundStmt:89*/
    typedef double(*Array2D)[n_cols];
    Array2D DST = (Array2D)dst, SRC = (Array2D)src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
        for (size_t i = 1; i < n_rows - 1; ++i) {
            for (size_t j = 1; j < n_cols - 1; ++j) {
                DST[i][j] = (SRC[i - 1][j] + SRC[i + 1][j] + SRC[i][j - 1] + SRC[i][j + 1]) / 4;
            }
        }
        swap_ptr((void**)&DST, (void**)&SRC);
    }
}
/*Function: stencil2D4pt_omp, ID: 12*/
void stencil2D4pt_omp(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps)
{
    /*stencil2D4pt_omp:12*/
    /*CompoundStmt:135*/
    typedef double(*Array2D)[n_cols];
    Array2D DST = (Array2D)dst, SRC = (Array2D)src;
    for (size_t ts = 0; ts < n_tsteps; ++ts) {
        /*CompoundStmt:144*/
        if (affinMaskRes) {
            myDARTSRuntime->run(launch<TP145>(ompNumThreads * DARTS_CODELETS_MULT, 0,
                RuntimeFinalCodelet, 1, n_rows - 1, (double**)&((DST)), (double**)&((SRC)),
                (size_t*)&((n_cols)), (size_t*)&((n_rows)), (size_t*)&((n_tsteps))));
        }
        swap_ptr((void**)&DST, (void**)&SRC);
    }
}
/*Function: stencil2D4pt_omp_v2, ID: 13*/
void stencil2D4pt_omp_v2(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps)
{
    /*stencil2D4pt_omp_v2:13*/
    /*CompoundStmt:206*/
    typedef double(*Array2D)[n_cols];
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP208>(ompNumThreads * DARTS_CODELETS_MULT, 0,
            RuntimeFinalCodelet, (double* __restrict*)&((dst)), (size_t*)&((n_cols)),
            (size_t*)&((n_rows)), (size_t*)&((n_tsteps)), (double* __restrict*)&((src))));
    }
}
/*Function: stencil2D4pt_omp_v3, ID: 14*/
void stencil2D4pt_omp_v3(double* restrict dst, double* restrict src, const size_t n_rows,
    const size_t n_cols, const size_t n_tsteps)
{
    /*stencil2D4pt_omp_v3:14*/
    /*CompoundStmt:280*/
    typedef double(*Array2D)[n_cols];
    char* str_n_threads = getenv("OMP_NUM_THREADS");
    size_t nThrds = str_n_threads ? strtoul(str_n_threads, ((void*)0), 0) : 1L;
    size_t nThrdsMinus1 = nThrds - 1;
    size_t chunk = (n_rows - 2) / nThrds;
    int* dCalc[nThrds];
    int* dSwap[nThrds];
    for (size_t i = 0; i < nThrds; ++i) {
        dCalc[i] = smalloc(sizeof(int) * (n_tsteps + 1));
        dSwap[i] = smalloc(sizeof(int) * (n_tsteps + 1));
    }
    int ts;
    if (affinMaskRes) {
        myDARTSRuntime->run(launch<TP313>(ompNumThreads * DARTS_CODELETS_MULT, 0,
            RuntimeFinalCodelet, (size_t*)&((chunk)), (int**)((dCalc)), (nThrds), (int**)((dSwap)),
            (nThrds), (double* __restrict*)&((dst)), (size_t*)&((nThrds)),
            (size_t*)&((nThrdsMinus1)), (size_t*)&((n_cols)), (size_t*)&((n_rows)),
            (size_t*)&((n_tsteps)), (double* __restrict*)&((src)), (int*)&((ts))));
    }
    for (size_t i = 0; i < nThrds; ++i) {
        sfree(dCalc[i]);
        sfree(dSwap[i]);
    }
}
/*Function: stencil_run, ID: 15*/
void* stencil_run(void* arg)
{
    /*stencil_run:15*/
    /*CompoundStmt:795*/
    stencil_t* stencil = (stencil_t*)arg;
    (stencil->stencil)((stencil->arg)->dst, (stencil->arg)->src, (stencil->arg)->n_rows,
        (stencil->arg)->n_cols, (stencil->arg)->n_tsteps);
    return (void*)((void*)0);
}
/*TP145: OMPParallelForDirective*/
void TP145::_barrierCodelets145::fire(void)
{
    TP145* myTP = static_cast<TP145*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
bool TP145::requestNewRangeIterations145(size_t* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        size_t tempStartRange = rangePerCodelet145 * codeletID;
        size_t tempEndRange = rangePerCodelet145 * (codeletID + 1);
        if (remainderRange145 != 0) {
            if (codeletID < (uint32_t)remainderRange145) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange145;
                tempEndRange += remainderRange145;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration145;
        tempEndRange = tempEndRange * 1 + minIteration145;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration145 < lastIteration145) {
            (this->inputsTPParent->i_darts145[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts145[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration145;
        }
    }
    return isThereNewIteration;
}
void TP145::_checkInCodelets146::fire(void)
{
    /*Init the vars for this region*/

	typedef double (*Array2D)[this->inputsTPParent->n_cols_darts145[this->getLocalID()]];
	
    /*printing node 146: ForStmt*/
    /*var: DST*/
    /*var: SRC*/
    /*var: n_cols*/
    /*var: n_rows*/
    /*var: n_tsteps*/
    Array2D* DST = (Array2D*)(this->inputsTPParent->DST_darts145);
    (void)DST /*OMP_SHARED*/;
    Array2D* SRC = (Array2D*)(this->inputsTPParent->SRC_darts145);
    (void)SRC /*OMP_SHARED*/;
    size_t* n_cols = &(this->inputsTPParent->n_cols_darts145[this->getLocalID()]);
    (void)n_cols /*OMP_FIRSTPRIVATE*/;
    size_t* i = &(this->inputsTPParent->i_darts145[this->getLocalID()]);
    (void)i /*PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations145(
        (size_t*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->barrierCodelets145[0].decDep();
        return;
    }
    for (size_t i_darts_counter_temp145 = (*i); i_darts_counter_temp145 < endRange
         && i_darts_counter_temp145 < this->inputsTPParent->lastIteration145;
         ++i_darts_counter_temp145) {
        {
            {
                /*Loop's init*/
                size_t j = 1;
                for (; j < (*n_cols) - 1; ++j) {
                    (*(DST))[(i_darts_counter_temp145)][j]
                        = ((*(SRC))[(i_darts_counter_temp145)-1][j]
                              + (*(SRC))[(i_darts_counter_temp145) + 1][j]
                              + (*(SRC))[(i_darts_counter_temp145)][j - 1]
                              + (*(SRC))[(i_darts_counter_temp145)][j + 1])
                        / 4;
                }
            }
        }
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->barrierCodelets145[0].decDep();
}
TP145::TP145(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet,
    size_t in_initIteration, size_t in_lastIteration, double** in_DST, double** in_SRC,
    size_t* in_n_cols_outer145_ptr, size_t* in_n_rows_outer145_ptr,
    size_t* in_n_tsteps_outer145_ptr)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , DST_darts145(in_DST) /*OMP_SHARED - INPUT*/
    , SRC_darts145(in_SRC) /*OMP_SHARED - INPUT*/
    , n_cols_darts145(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_outer145_ptr(in_n_cols_outer145_ptr)
    , n_rows_darts145(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_outer145_ptr(in_n_rows_outer145_ptr)
    , n_tsteps_darts145(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_tsteps_outer145_ptr(in_n_tsteps_outer145_ptr)
    , i_darts145(new size_t[this->numThreads]) /*VARIABLE*/
    , initIteration145(in_initIteration)
    , lastIteration145(in_lastIteration)
    , barrierCodelets145(new _barrierCodelets145[1])
    , checkInCodelets146(new _checkInCodelets146[this->numThreads])
{
    /*Initialize the loop parameters*/
    range145 = abs(lastIteration145 - initIteration145) / 1;
    rangePerCodelet145 = range145 / numThreads;
    minIteration145 = min<size_t>(lastIteration145, initIteration145);
    remainderRange145 = range145 % numThreads;
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets145[0] = _barrierCodelets145(ompNumThreads, ompNumThreads, this, 0);
    _checkInCodelets146* checkInCodelets146Ptr = (this->checkInCodelets146);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        n_cols_darts145[codeletCounter] = *n_cols_outer145_ptr; /*OMP_FIRSTPRIVATE*/
        n_rows_darts145[codeletCounter] = *n_rows_outer145_ptr; /*OMP_FIRSTPRIVATE*/
        n_tsteps_darts145[codeletCounter] = *n_tsteps_outer145_ptr; /*OMP_FIRSTPRIVATE*/
        (*checkInCodelets146Ptr) = _checkInCodelets146(1, 1, this, codeletCounter);
        (*checkInCodelets146Ptr).decDep();
        checkInCodelets146Ptr++;
    }
}
TP145::~TP145()
{
    delete[] i_darts145;
    delete[] barrierCodelets145;
    delete[] checkInCodelets146;
}
/*TP208: OMPParallelDirective*/
void TP208::_barrierCodelets208::fire(void)
{
    TP208* myTP = static_cast<TP208*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP208::_checkInCodelets210::fire(void)
{
    /*Init the vars for this region*/

    /*printing node 210: DeclStmt*/
    this->inputsTPParent->DST_darts208[this->getID()]
        = (double*)(*(this->inputsTPParent->dst_darts208));
    this->inputsTPParent->SRC_darts208[this->getID()]
        = (double*)(*(this->inputsTPParent->src_darts208));

    /*printing node 213: DeclStmt*/
    this->inputsTPParent->n_ts_darts208[this->getID()]
        = (this->inputsTPParent->n_tsteps_darts208[this->getID()]);

    /*printing node 215: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->n_ts_darts208[this->getID()])-- > 0) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets214[this->getID()].decDep();
        return;
    }
    /*Signal the codelet after the loop from the end's condional node.*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->barrierCodelets208[0].decDep();
}
void TP208::_checkInCodelets214::fire(void)
{

    /*printing node 214: WhileStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP214_LoopCounter),
        myTP->controlTPParent->TP214_LoopCounterPerThread[this->getLocalID()],
        myTP->controlTPParent->TP214_LoopCounterPerThread[this->getLocalID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP214_LoopCounterPerThread[this->getLocalID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP214PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP214_LoopCounterPerThread[this->getLocalID()] += 1;
        invoke<TP214>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP214PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP214PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP214PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP214PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP214PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP214PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP214_LoopCounterPerThread[this->getLocalID()] += 1;
        }
    }
}
void TP208::_checkInCodelets817::fire(void)
{

    /*printing node 817: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->n_ts_darts208[this->getID()])-- > 0) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets214[this->getID()].decDep();
        return;
    }
    /*Signal the codelet after the loop from the condtional node.*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->barrierCodelets208[0].decDep();
}
TP208::TP208(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet,
    double* __restrict* in_dst, size_t* in_n_cols_outer208_ptr, size_t* in_n_rows_outer208_ptr,
    size_t* in_n_tsteps_outer208_ptr, double* __restrict* in_src)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , dst_darts208(in_dst) /*OMP_SHARED - INPUT*/
    , n_cols_darts208(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_outer208_ptr(in_n_cols_outer208_ptr)
    , n_rows_darts208(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_outer208_ptr(in_n_rows_outer208_ptr)
    , n_tsteps_darts208(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_tsteps_outer208_ptr(in_n_tsteps_outer208_ptr)
    , src_darts208(in_src) /*OMP_SHARED - INPUT*/
    , DST_darts208(new double*[this->numThreads]) /*VARIABLE*/
    , SRC_darts208(new double*[this->numThreads]) /*VARIABLE*/
    , n_ts_darts208(new size_t[this->numThreads]) /*VARIABLE*/
    , TP214_LoopCounter(0)
    , TP214_LoopCounterPerThread(new unsigned int[this->numThreads])
    , barrierCodelets208(new _barrierCodelets208[1])
    , checkInCodelets210(new _checkInCodelets210[this->numThreads])
    , checkInCodelets214(new _checkInCodelets214[this->numThreads])
    , checkInCodelets817(new _checkInCodelets817[this->numThreads])
{
    memset((void*)TP214_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets208[0] = _barrierCodelets208(ompNumThreads, ompNumThreads, this, 0);
    _checkInCodelets817* checkInCodelets817Ptr = (this->checkInCodelets817);
    _checkInCodelets214* checkInCodelets214Ptr = (this->checkInCodelets214);
    _checkInCodelets210* checkInCodelets210Ptr = (this->checkInCodelets210);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        n_cols_darts208[codeletCounter] = *n_cols_outer208_ptr; /*OMP_FIRSTPRIVATE*/
        n_rows_darts208[codeletCounter] = *n_rows_outer208_ptr; /*OMP_FIRSTPRIVATE*/
        n_tsteps_darts208[codeletCounter] = *n_tsteps_outer208_ptr; /*OMP_FIRSTPRIVATE*/
        (*checkInCodelets817Ptr) = _checkInCodelets817(1, 1, this, codeletCounter);
        checkInCodelets817Ptr++;
        (*checkInCodelets214Ptr) = _checkInCodelets214(1, 1, this, codeletCounter);
        checkInCodelets214Ptr++;
        (*checkInCodelets210Ptr) = _checkInCodelets210(1, 1, this, codeletCounter);
        (*checkInCodelets210Ptr).decDep();
        checkInCodelets210Ptr++;
    }
}
TP208::~TP208()
{
    delete[] TP214_LoopCounterPerThread;
    delete[] DST_darts208;
    delete[] SRC_darts208;
    delete[] n_ts_darts208;
    delete[] barrierCodelets208;
    delete[] checkInCodelets817;
    delete[] checkInCodelets214;
    delete[] checkInCodelets210;
}
/*TP214: WhileStmt*/
void TP214::_checkInCodelets218::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*region 218 0*/
    /*Determine the TP to which this codelet belongs and check if this codelet spawns the TP or if
     * it signals it using dispatchCodelet()*/
    size_t idx = this->getID() / myTP->codeletsPerTP218;
    if (idx < myTP->TPsToUse218) {
        if (!__sync_val_compare_and_swap(&(myTP->TP218_alreadyLaunched[idx]), 0, 1)) {
            size_t range = abs((this->inputsTPParent->n_rows_darts208[this->getID()]) - 1 - 1) / 1;
            size_t rangePerCodelet = range / myTP->TPsToUse218;
            size_t minIteration
                = min<size_t>((this->inputsTPParent->n_rows_darts208[this->getID()]) - 1, 1);
            size_t remainderRange = range % myTP->TPsToUse218;
            size_t initIteration = rangePerCodelet * idx;
            size_t lastIteration = rangePerCodelet * (idx + 1);
            if (remainderRange != 0) {
                if (idx < (uint32_t)remainderRange) {
                    initIteration += idx;
                    lastIteration += (idx + 1);
                } else {
                    initIteration += remainderRange;
                    lastIteration += remainderRange;
                }
            }
            initIteration = initIteration * 1 + minIteration;
            lastIteration = lastIteration * 1 + minIteration;
            if (1 < (this->inputsTPParent->n_rows_darts208[this->getID()]) - 1) {
                initIteration = min(initIteration, lastIteration);
                lastIteration = max(initIteration, lastIteration);
            } else {
                initIteration = max(initIteration, lastIteration);
                lastIteration = min(initIteration, lastIteration);
            }
            if (idx == myTP->TPsToUse218 - 1) {
                lastIteration = (this->inputsTPParent->n_rows_darts208[this->getID()]) - 1;
            }
#if USEINVOKE == 1
            invoke<TP218>(myTP, myTP->codeletsPerTP218 * DARTS_CODELETS_MULT, this->getID(), myTP,
                initIteration, lastIteration,
                &((this->inputsTPParent->DST_darts208[this->getID()])),
                &((this->inputsTPParent->SRC_darts208[this->getID()])),
                &((this->inputsTPParent->n_cols_darts208[this->getID()])),
                &((this->inputsTPParent->n_rows_darts208[this->getID()])), &(myTP->TP218Ptr[idx]));
#else
            place<TP218>(idx, myTP, myTP->codeletsPerTP218 * DARTS_CODELETS_MULT, this->getID(),
                myTP, initIteration, lastIteration,
                &((this->inputsTPParent->DST_darts208[this->getID()])),
                &((this->inputsTPParent->SRC_darts208[this->getID()])),
                &((this->inputsTPParent->n_cols_darts208[this->getID()])),
                &((this->inputsTPParent->n_rows_darts208[this->getID()])), &(myTP->TP218Ptr[idx]));
#endif
        } else {
            if (myTP->TP218Ptr[idx] != nullptr) {
                myTP->TP218Ptr[idx]->dispatchCodelet(this->getID());
            } else {
                this->resetCodelet();
                this->decDep();
            }
        }
    } else {
        /*Signaling next codelet region: 218 nextRegion: 274 */
        myTP->controlTPParent->checkInCodelets274[this->getID()].decDep();
    }
}
void TP214::_checkInCodelets274::fire(void)
{

    /*printing node 274: CallExpr*/
    swap_ptr((void**)&(this->inputsTPParent->DST_darts208[this->getID()]),
        (void**)&(this->inputsTPParent->SRC_darts208[this->getID()]));
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 274 nextRegion: 279 */
    myTP->controlTPParent->barrierCodelets279[0].decDep();
}
void TP214::_barrierCodelets279::fire(void)
{
    TP214* myTP = static_cast<TP214*>(myTP_);
    for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
        myTP->TPParent->checkInCodelets817[codeletsCounter].decDep();
    }
}
TP214::TP214(int in_numThreads, int in_mainCodeletID, TP208* in_TPParent, TP208* in_inputsTPParent,
    TP214** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP218Ptr(new TP218*[NUMTPS218])
    , TP218_alreadyLaunched(new size_t[NUMTPS218])
    , numTPsSet218(0)
    , numTPsReady218(0)
    , TPsToUse218(NUMTPS218)
    , codeletsPerTP218(this->numThreads / NUMTPS218)
    , totalCodelets218(this->TPsToUse218 * this->codeletsPerTP218)
    , checkInCodelets218(new _checkInCodelets218[this->numThreads])
    , checkInCodelets274(new _checkInCodelets274[this->numThreads])
    , barrierCodelets279(new _barrierCodelets279[1])
{
    /*Initialize Codelets*/
    barrierCodelets279[0] = _barrierCodelets279(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets274* checkInCodelets274Ptr = (this->checkInCodelets274);
    _checkInCodelets218* checkInCodelets218Ptr = (this->checkInCodelets218);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets218);
#endif
    for (int i = 0; i < NUMTPS218; i++) {
        TP218Ptr[i] = nullptr;
        TP218_alreadyLaunched[i] = 0;
    }
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        (*checkInCodelets274Ptr) = _checkInCodelets274(1, 1, this, codeletCounter);
        checkInCodelets274Ptr++;
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets218Ptr) = _checkInCodelets218(2, 1, this, codeletCounter);
#else
        (*checkInCodelets218Ptr) = _checkInCodelets218(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets218Ptr).decDep();
        checkInCodelets218Ptr++;
    }
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP214::~TP214()
{
    delete[] barrierCodelets279;
    delete[] checkInCodelets274;
    delete[] checkInCodelets218;
}
/*TP218: OMPForDirective*/
bool TP218::requestNewRangeIterations218(size_t* endRange, uint32_t codeletID)
{
    /*Scheduling Policy = Static */
    /*Chunk = 0*/
    bool isThereNewIteration = false;
    {
        /*Static Scheduling*/
        size_t tempStartRange = rangePerCodelet218 * codeletID;
        size_t tempEndRange = rangePerCodelet218 * (codeletID + 1);
        if (remainderRange218 != 0) {
            if (codeletID < (uint32_t)remainderRange218) {
                tempStartRange += codeletID;
                tempEndRange += (codeletID + 1);
            } else {
                tempStartRange += remainderRange218;
                tempEndRange += remainderRange218;
            }
        }
        tempStartRange = tempStartRange * 1 + minIteration218;
        tempEndRange = tempEndRange * 1 + minIteration218;
        if (tempStartRange != tempEndRange) {
            isThereNewIteration = true;
        }
        if (initIteration218 < lastIteration218) {
            (this->inputsTPParent->i_darts218[codeletID]) = min(tempStartRange, tempEndRange);
            *endRange = max(tempStartRange, tempEndRange);
        } else {
            (this->inputsTPParent->i_darts218[codeletID]) = max(tempStartRange, tempEndRange);
            *endRange = min(tempStartRange, tempEndRange);
        }
        if (codeletID == this->numThreads - 1) {
            *endRange = lastIteration218;
        }
    }
    return isThereNewIteration;
}
void TP218::_checkInCodelets219::fire(void)
{
#if USE_SPIN_CODELETS == 1
    /*Wait until the codelet with the same ID finishes in the previous TP*/
    if (myTP->availableCodelets[this->getLocalID()] == 0) {
        myTP->add(this);
        return;
    }
#endif
    /*Init the vars for this region*/
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->DST_darts218[this->getLocalID()]
        = (double**)&(myTP->TPParent->inputsTPParent->DST_darts208[this->getID()]);
    /*Get pointer from parent for variable
     with shared scope in this region but private
     in the enclosing one.*/
    this->inputsTPParent->SRC_darts218[this->getLocalID()]
        = (double**)&(myTP->TPParent->inputsTPParent->SRC_darts208[this->getID()]);

	typedef double (*Array2D)[this->inputsTPParent->n_cols_darts218[this->getLocalID()]];
		
    /*printing node 219: ForStmt*/
    /*var: DST*/
    /*var: SRC*/
    /*var: n_cols*/
    /*var: n_rows*/
    Array2D** DST = (Array2D**)&(this->inputsTPParent->DST_darts218[this->getLocalID()]);
    (void)DST /*OMP_SHARED_PRIVATE*/;
    Array2D** SRC = (Array2D**)&(this->inputsTPParent->SRC_darts218[this->getLocalID()]);
    (void)SRC /*OMP_SHARED_PRIVATE*/;
    size_t* n_cols = &(this->inputsTPParent->n_cols_darts218[this->getLocalID()]);
    (void)n_cols /*OMP_FIRSTPRIVATE*/;
    size_t* i = &(this->inputsTPParent->i_darts218[this->getLocalID()]);
    (void)i /*PRIVATE*/;
    bool isThereNewIteration = this->inputsTPParent->requestNewRangeIterations218(
        (size_t*)&(this->endRange), this->getLocalID());
    if (isThereNewIteration == false) {
        /*Find and signal the next codelet*/
        myTP->controlTPParent->TPParent->checkInCodelets274[this->getID()].decDep();
        return;
    }
    for (size_t i_darts_counter_temp218 = (*i); i_darts_counter_temp218 < endRange
         && i_darts_counter_temp218 < this->inputsTPParent->lastIteration218;
         ++i_darts_counter_temp218) {
        {
            {
                /*Loop's init*/
                size_t j = 1;
                for (; j < (*n_cols) - 1; ++j) {
                    (*(*DST))[(i_darts_counter_temp218)][j]
                        = ((*(*SRC))[(i_darts_counter_temp218)-1][j]
                              + (*(*SRC))[(i_darts_counter_temp218) + 1][j]
                              + (*(*SRC))[(i_darts_counter_temp218)][j - 1]
                              + (*(*SRC))[(i_darts_counter_temp218)][j + 1])
                        / 4;
                }
            }
        }
    }
    /*If this omp for has no barrier,
    check if all the codelets
    replicated from the same
    global ID has finished and
    signal the next codelet.
    Otherwise, return.*/
    uint32_t completedMultCodelet = __sync_fetch_and_add(
        &(myTP->signalNextReady[this->getLocalID() % myTP->baseNumThreads]), 1);
    if (completedMultCodelet < (uint32_t)(DARTS_CODELETS_MULT - 1))
        return;
    /*Signaling next codelet from last stmt in the codelet*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->checkInCodelets274[this->getID()].decDep();
}
TP218::TP218(int in_numThreads, int in_mainCodeletID, TP214* in_TPParent, size_t in_initIteration,
    size_t in_lastIteration, double** in_DST, double** in_SRC, size_t* in_n_cols_outer218_ptr,
    size_t* in_n_rows_outer218_ptr, TP218** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , DST_darts218(new double**[this->numThreads]) /*OMP_SHARED_PRIVATE - INPUT*/
    , SRC_darts218(new double**[this->numThreads]) /*OMP_SHARED_PRIVATE - INPUT*/
    , n_cols_darts218(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_outer218_ptr(in_n_cols_outer218_ptr)
    , n_rows_darts218(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_outer218_ptr(in_n_rows_outer218_ptr)
    , i_darts218(new size_t[this->numThreads]) /*VARIABLE*/
    , initIteration218(in_initIteration)
    , lastIteration218(in_lastIteration)
    , readyCodelets(0)
    , baseNumThreads(this->numThreads / DARTS_CODELETS_MULT)
    , signalNextReady(new int[baseNumThreads])
    , checkInCodelets219(new _checkInCodelets219[this->numThreads])
{
    /*Initialize the loop parameters*/
    range218 = abs(lastIteration218 - initIteration218) / 1;
    rangePerCodelet218 = range218 / numThreads;
    minIteration218 = min<size_t>(lastIteration218, initIteration218);
    remainderRange218 = range218 % numThreads;
    /*Initialize inputs and vars.*/
    this->DST_darts218
        = (double***)malloc(sizeof(double**) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    this->SRC_darts218
        = (double***)malloc(sizeof(double**) * this->numThreads) /*OMP_SHARED_PRIVATE*/;
    /*Initialize Codelets*/
    _checkInCodelets219* checkInCodelets219Ptr = (this->checkInCodelets219);
#if USE_SPIN_CODELETS == 0
    firstCodelet = (this->checkInCodelets219);
#endif
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->baseNumThreads;
         codeletCounter++) {
        this->signalNextReady[codeletCounter] = 0;
        n_cols_darts218[codeletCounter] = *n_cols_outer218_ptr; /*OMP_FIRSTPRIVATE*/
        for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
            n_cols_darts218[codeletCounter + this->baseNumThreads * i] = *n_cols_outer218_ptr;
        }
        n_rows_darts218[codeletCounter] = *n_rows_outer218_ptr; /*OMP_FIRSTPRIVATE*/
        for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
            n_rows_darts218[codeletCounter + this->baseNumThreads * i] = *n_rows_outer218_ptr;
        }
#if USE_SPIN_CODELETS == 0
        (*checkInCodelets219Ptr) = _checkInCodelets219(2, 1, this, codeletCounter);
#else
        (*checkInCodelets219Ptr) = _checkInCodelets219(1, 1, this, codeletCounter);
#endif
        (*checkInCodelets219Ptr).decDep();
        checkInCodelets219Ptr++;
    }
    this->dispatchCodelet(this->getID());
    *(in_ptrToThisTP) = this;
}
void TP218::dispatchCodelet(size_t codeletID)
{
    int idx = codeletID / this->baseNumThreads;
    int localID = codeletID - this->baseNumThreads * idx;
    this->checkInCodelets219[localID].setID(codeletID);
    this->checkInCodelets219[localID].setLocalID(localID);
    /*Check if we want to replicate codelets*/
    for (size_t i = 1; i < (size_t)DARTS_CODELETS_MULT; i++) {
#if USE_SPIN_CODELETS == 0
        this->checkInCodelets219[localID + this->baseNumThreads * i]
            = _checkInCodelets219(2, 1, this, localID + this->baseNumThreads * i);
#else
        this->checkInCodelets219[localID + this->baseNumThreads * i]
            = _checkInCodelets219(1, 1, this, localID + this->baseNumThreads * i);
#endif
        this->checkInCodelets219[localID + this->baseNumThreads * i].setID(codeletID);
        this->checkInCodelets219[localID + this->baseNumThreads * i].decDep();
#if USE_SPIN_CODELETS == 0
        this->firstCodelet[localID + this->baseNumThreads * i].decDep();
#else
        this->availableCodelets[localID + this->baseNumThreads * i] = 1;
#endif
    }
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[localID].decDep();
#else
    this->availableCodelets[localID] = 1;
#endif
}
TP218::~TP218()
{
    delete[] DST_darts218;
    delete[] SRC_darts218;
    delete[] i_darts218;
    delete[] checkInCodelets219;
}
/*TP313: OMPParallelDirective*/
void TP313::_barrierCodelets313::fire(void)
{
    TP313* myTP = static_cast<TP313*>(myTP_);
    myTP->controlTPParent->nextCodelet->decDep();
}
void TP313::_checkInCodelets315::fire(void)
{
    /*Select the thread executing OMPSingleDirective 315*/
    if (!__sync_val_compare_and_swap(&(myTP->TP315_alreadyLaunched), 0, 1)) {
        invoke<TP315>(myTP, 1, this->getID(), myTP,
            &((this->inputsTPParent->chunk_darts313[this->getID()])),
            ((this->inputsTPParent->dCalc_darts313)), (this->inputsTPParent->dCalc_outer313_size),
            ((this->inputsTPParent->dSwap_darts313)), (this->inputsTPParent->dSwap_outer313_size),
            &(*(this->inputsTPParent->dst_darts313)),
            &((this->inputsTPParent->nThrdsMinus1_darts313[this->getID()])),
            &((this->inputsTPParent->n_cols_darts313[this->getID()])),
            &((this->inputsTPParent->n_rows_darts313[this->getID()])),
            &((this->inputsTPParent->n_tsteps_darts313[this->getID()])),
            &(*(this->inputsTPParent->src_darts313)), &(*(this->inputsTPParent->ts_darts313)));
    } else {
        myTP->barrierCodelets315[0].decDep();
    }
}
void TP313::_barrierCodelets315::fire(void)
{
    TP313* myTP = static_cast<TP313*>(myTP_);
    for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
        myTP->TPParent->barrierCodelets313[0].decDep();
    }
}
TP313::TP313(int in_numThreads, int in_mainCodeletID, Codelet* in_nextCodelet,
    size_t* in_chunk_outer313_ptr, int** in_dCalc, int in_dCalc_outer313_size, int** in_dSwap,
    int in_dSwap_outer313_size, double* __restrict* in_dst, size_t* in_nThrds_outer313_ptr,
    size_t* in_nThrdsMinus1_outer313_ptr, size_t* in_n_cols_outer313_ptr,
    size_t* in_n_rows_outer313_ptr, size_t* in_n_tsteps_outer313_ptr, double* __restrict* in_src,
    int* in_ts)
    : ThreadedProcedure(in_numThreads, in_mainCodeletID)
    , nextCodelet(in_nextCodelet)
    , TPParent(this)
    , controlTPParent(this)
    , inputsTPParent(this)
    , chunk_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , chunk_outer313_ptr(in_chunk_outer313_ptr)
    , dCalc_darts313(in_dCalc) /*OMP_SHARED - INPUT*/
    , dCalc_outer313_size(in_dCalc_outer313_size)
    , dSwap_darts313(in_dSwap) /*OMP_SHARED - INPUT*/
    , dSwap_outer313_size(in_dSwap_outer313_size)
    , dst_darts313(in_dst) /*OMP_SHARED - INPUT*/
    , nThrds_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , nThrds_outer313_ptr(in_nThrds_outer313_ptr)
    , nThrdsMinus1_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , nThrdsMinus1_outer313_ptr(in_nThrdsMinus1_outer313_ptr)
    , n_cols_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_outer313_ptr(in_n_cols_outer313_ptr)
    , n_rows_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_outer313_ptr(in_n_rows_outer313_ptr)
    , n_tsteps_darts313(new size_t[this->numThreads]) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_tsteps_outer313_ptr(in_n_tsteps_outer313_ptr)
    , src_darts313(in_src) /*OMP_SHARED - INPUT*/
    , ts_darts313(in_ts) /*OMP_SHARED - INPUT*/
    , TP315Ptr(nullptr)
    , TP315_alreadyLaunched(0)
    , barrierCodelets313(new _barrierCodelets313[1])
    , checkInCodelets315(new _checkInCodelets315[this->numThreads])
    , barrierCodelets315(new _barrierCodelets315[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets313[0] = _barrierCodelets313(ompNumThreads, ompNumThreads, this, 0);
    barrierCodelets315[0] = _barrierCodelets315(this->numThreads, this->numThreads, this, 0);
    _checkInCodelets315* checkInCodelets315Ptr = (this->checkInCodelets315);
    for (size_t codeletCounter = 0; codeletCounter < (size_t)this->numThreads; codeletCounter++) {
        chunk_darts313[codeletCounter] = *chunk_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        nThrds_darts313[codeletCounter] = *nThrds_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        nThrdsMinus1_darts313[codeletCounter] = *nThrdsMinus1_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        n_cols_darts313[codeletCounter] = *n_cols_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        n_rows_darts313[codeletCounter] = *n_rows_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        n_tsteps_darts313[codeletCounter] = *n_tsteps_outer313_ptr; /*OMP_FIRSTPRIVATE*/
        (*checkInCodelets315Ptr) = _checkInCodelets315(1, 1, this, codeletCounter);
        (*checkInCodelets315Ptr).decDep();
        checkInCodelets315Ptr++;
    }
}
TP313::~TP313()
{
    delete[] barrierCodelets313;
    delete[] barrierCodelets315;
    delete[] checkInCodelets315;
}
/*TP315: OMPSingleDirective*/
void TP315::_checkInCodelets317::fire(void)
{
    invoke<TP317>(myTP, myTP->numThreads, this->getID(), myTP,
        &((this->inputsTPParent->chunk_darts315)), ((this->inputsTPParent->dCalc_darts315)),
        (this->inputsTPParent->dCalc_outer315_size), ((this->inputsTPParent->dSwap_darts315)),
        (this->inputsTPParent->dSwap_outer315_size), &(*(this->inputsTPParent->dst_darts315)),
        &((this->inputsTPParent->nThrdsMinus1_darts315)),
        &((this->inputsTPParent->n_cols_darts315)), &((this->inputsTPParent->n_rows_darts315)),
        &((this->inputsTPParent->n_tsteps_darts315)), &(*(this->inputsTPParent->src_darts315)));
}
void TP315::_barrierCodelets317::fire(void)
{	
    TP315* myTP = static_cast<TP315*>(myTP_);
    for (size_t codeletsCounter = 0; codeletsCounter < myTP->numThreads; codeletsCounter++) {
        myTP->TPParent->barrierCodelets315[0].decDep();
    }
}
TP315::TP315(int in_numThreads, int in_mainCodeletID, TP313* in_TPParent,
    size_t* in_chunk_outer315_ptr, int** in_dCalc, int in_dCalc_outer315_size, int** in_dSwap,
    int in_dSwap_outer315_size, double* __restrict* in_dst, size_t* in_nThrdsMinus1_outer315_ptr,
    size_t* in_n_cols_outer315_ptr, size_t* in_n_rows_outer315_ptr,
    size_t* in_n_tsteps_outer315_ptr, double* __restrict* in_src, int* in_ts)
    : ompOMPSingleDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , chunk_darts315(*in_chunk_outer315_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , dCalc_darts315(in_dCalc) /*OMP_SHARED - INPUT*/
    , dCalc_outer315_size(in_dCalc_outer315_size)
    , dSwap_darts315(in_dSwap) /*OMP_SHARED - INPUT*/
    , dSwap_outer315_size(in_dSwap_outer315_size)
    , dst_darts315(in_dst) /*OMP_SHARED - INPUT*/
    , nThrdsMinus1_darts315(*in_nThrdsMinus1_outer315_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_darts315(*in_n_cols_outer315_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_darts315(*in_n_rows_outer315_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_tsteps_darts315(*in_n_tsteps_outer315_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , src_darts315(in_src) /*OMP_SHARED - INPUT*/
    , ts_darts315(in_ts) /*OMP_SHARED - INPUT*/
    , checkInCodelets317(1, 1, this, this->mainCodeletID)
    , barrierCodelets317(new _barrierCodelets317[1])
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    barrierCodelets317[0] = _barrierCodelets317(this->numThreads, this->numThreads, this, 0);
    checkInCodelets317.setLocalID(0);
    checkInCodelets317.decDep();
}
TP315::~TP315() { delete[] barrierCodelets317; }
/*TP317: OMPTaskgroupDirective*/
void TP317::_checkInCodelets320::fire(void)
{
	
    /*Init the vars for this region*/

    /*printing node 320: BinaryOperator*/
    (this->inputsTPParent->ts_darts317) = (this->inputsTPParent->n_tsteps_darts317);

    /*printing node 321: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->ts_darts317) > 1) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets319.decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 322 nextRegion: 628 */
        myTP->controlTPParent->checkInCodelets628.decDep();
        return;
    }
}
void TP317::_checkInCodelets319::fire(void)
{

    /*printing node 319: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP319_LoopCounter),
        myTP->controlTPParent->TP319_LoopCounterPerThread[this->getLocalID()],
        myTP->controlTPParent->TP319_LoopCounterPerThread[this->getLocalID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP319_LoopCounterPerThread[this->getLocalID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP319PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP319_LoopCounterPerThread[this->getLocalID()] += 1;
        invoke<TP319>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP319PtrVec.back()));
    } else {

        if (myTP->controlTPParent->TP319PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP319PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP319PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP319PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP319PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP319_LoopCounterPerThread[this->getLocalID()] += 1;
        }
    }
}
void TP317::_checkInCodelets322::fire(void)
{

    /*printing node 322: BinaryOperator*/
    (this->inputsTPParent->ts_darts317) = (this->inputsTPParent->ts_darts317) - 2;

    /*printing node 815: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->ts_darts317) > 1) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets319.decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 322 nextRegion: 628 */
        myTP->controlTPParent->checkInCodelets628.decDep();
        return;
    }
}
void TP317::_checkInCodelets628::fire(void)
{
	
    /*Printing conditional branch node 628: inlining: 0*/
    if ((this->inputsTPParent->ts_darts317) == 1) {
		
        myTP->checkInCodelets631.decDep();
    } else {
        /*Signaling the region after the if stmt*/
        /*Find and signal the next codelet*/	
        myTP->controlTPParent->TPParent->barrierCodelets317[0].decDep();
    }
}
void TP317::_checkInCodelets631::fire(void)
{

    /*printing node 631: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 631*/
    _task631Inputs* task631Inputs = new _task631Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->n_cols_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task631Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[0][(this->inputsTPParent->ts_darts317) + 2]));
    task631Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[1][(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task631Inputs->taskOutVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[0][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task631Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task631ListLock.lock();
    list<_task631Inputs*>::iterator task631It
        = task631List[this->getID()].insert(task631List[this->getID()].end(), task631Inputs);
    task631ListLock.unlock();
    invoke<TP631>(myTP, 1, this->getID(), myTP, task631Inputs, &(task631List[this->getID()]),
        &(task631ListLock), task631It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 631 nextRegion: 677 */
	myTP->controlTPParent->checkInCodelets677.decDep();
}
void TP317::_checkInCodelets677::fire(void)
{
    /*printing node 677: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 677*/
    _task677Inputs* task677Inputs = new _task677Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->nThrdsMinus1_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &((this->inputsTPParent->n_rows_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task677Inputs->taskInVarDependencies.push_back(&(
        ((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->nThrdsMinus1_darts317) - 1]
                                                [(this->inputsTPParent->ts_darts317) + 2]));
    task677Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task677Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task677Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task677ListLock.lock();
    list<_task677Inputs*>::iterator task677It
        = task677List[this->getID()].insert(task677List[this->getID()].end(), task677Inputs);
    task677ListLock.unlock();
    invoke<TP677>(myTP, 1, this->getID(), myTP, task677Inputs, &(task677List[this->getID()]),
        &(task677ListLock), task677It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 677 nextRegion: 727 */
    myTP->controlTPParent->checkInCodelets727.decDep();
}
void TP317::_checkInCodelets727::fire(void)
{
    /*printing node 727: DeclStmt*/
    this->inputsTPParent->k_darts317 = 1;

    /*printing node 728: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets726.decDep();
        return;
    }
    /*Signal the codelet after the loop from the end's condional node.*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->barrierCodelets317[0].decDep();
}
void TP317::_checkInCodelets726::fire(void)
{

    /*printing node 726: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP726_LoopCounter),
        myTP->controlTPParent->TP726_LoopCounterPerThread[this->getLocalID()],
        myTP->controlTPParent->TP726_LoopCounterPerThread[this->getLocalID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP726_LoopCounterPerThread[this->getLocalID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP726PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP726_LoopCounterPerThread[this->getLocalID()] += 1;
        invoke<TP726>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP726PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP726PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP726PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP726PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP726PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP726PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP726_LoopCounterPerThread[this->getLocalID()] += 1;
        }
    }
}
void TP317::_checkInCodelets729::fire(void)
{
    /*printing node 729: UnaryOperator*/
    ++(this->inputsTPParent->k_darts317);

    /*printing node 816: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets726.decDep();
        return;
    }
    /*Signal the codelet after the loop from the condtional node.*/
    /*Find and signal the next codelet*/
    myTP->controlTPParent->TPParent->barrierCodelets317[0].decDep();
}
TP317::TP317(int in_numThreads, int in_mainCodeletID, TP315* in_TPParent,
    size_t* in_chunk_outer317_ptr, int** in_dCalc, int in_dCalc_outer317_size, int** in_dSwap,
    int in_dSwap_outer317_size, double* __restrict* in_dst, size_t* in_nThrdsMinus1_outer317_ptr,
    size_t* in_n_cols_outer317_ptr, size_t* in_n_rows_outer317_ptr,
    size_t* in_n_tsteps_outer317_ptr, double* __restrict* in_src)
    : ompOMPTaskgroupDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , chunk_darts317(*in_chunk_outer317_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , dCalc_darts317(in_dCalc) /*OMP_SHARED - INPUT*/
    , dCalc_outer317_size(in_dCalc_outer317_size)
    , dSwap_darts317(in_dSwap) /*OMP_SHARED - INPUT*/
    , dSwap_outer317_size(in_dSwap_outer317_size)
    , dst_darts317(in_dst) /*OMP_SHARED - INPUT*/
    , nThrdsMinus1_darts317(*in_nThrdsMinus1_outer317_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_cols_darts317(*in_n_cols_outer317_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_rows_darts317(*in_n_rows_outer317_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , n_tsteps_darts317(*in_n_tsteps_outer317_ptr) /*OMP_FIRSTPRIVATE - INPUT*/
    , src_darts317(in_src) /*OMP_SHARED - INPUT*/
    , TP319_LoopCounter(0)
    , TP319_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP726_LoopCounter(0)
    , TP726_LoopCounterPerThread(new unsigned int[this->numThreads])
    , checkInCodelets320(1, 1, this, this->mainCodeletID)
    , checkInCodelets319(1, 1, this, this->mainCodeletID)
    , checkInCodelets322(1, 1, this, this->mainCodeletID)
    , checkInCodelets628(1, 1, this, this->mainCodeletID)
    , checkInCodelets631(1, 1, this, this->mainCodeletID)
    , checkInCodelets677(1, 1, this, this->mainCodeletID)
    , checkInCodelets727(1, 1, this, this->mainCodeletID)
    , checkInCodelets726(1, 1, this, this->mainCodeletID)
    , checkInCodelets729(1, 1, this, this->mainCodeletID)
{
    memset((void*)TP319_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    memset((void*)TP726_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets729.setLocalID(0);
    checkInCodelets726.setLocalID(0);
    checkInCodelets727.setLocalID(0);
    checkInCodelets677.setLocalID(0);
    checkInCodelets631.setLocalID(0);
    checkInCodelets628.setLocalID(0);
    checkInCodelets322.setLocalID(0);
    checkInCodelets319.setLocalID(0);
    checkInCodelets320.setLocalID(0);
    checkInCodelets320.decDep();
}
TP317::~TP317()
{
    delete[] TP319_LoopCounterPerThread;
    delete[] TP726_LoopCounterPerThread;
}
/*TP319: ForStmt*/
void TP319::_checkInCodelets325::fire(void)
{
    /*printing node 325: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 325*/
    _task325Inputs* task325Inputs = new _task325Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->n_cols_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task325Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[0][(this->inputsTPParent->ts_darts317) + 2]));
    task325Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[1][(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task325Inputs->taskOutVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[0][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task325Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task325ListLock.lock();
    list<_task325Inputs*>::iterator task325It
        = task325List[this->getID()].insert(task325List[this->getID()].end(), task325Inputs);
    task325ListLock.unlock();
    invoke<TP325>(myTP, 1, this->getID(), myTP, task325Inputs, &(task325List[this->getID()]),
        &(task325ListLock), task325It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 325 nextRegion: 371 */
    myTP->controlTPParent->checkInCodelets371.decDep();
}
void TP319::_checkInCodelets371::fire(void)
{
    /*printing node 371: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 371*/
    _task371Inputs* task371Inputs = new _task371Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->nThrdsMinus1_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &((this->inputsTPParent->n_rows_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task371Inputs->taskInVarDependencies.push_back(&(
        ((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->nThrdsMinus1_darts317) - 1]
                                                [(this->inputsTPParent->ts_darts317) + 2]));
    task371Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task371Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task371Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task371ListLock.lock();
    list<_task371Inputs*>::iterator task371It
        = task371List[this->getID()].insert(task371List[this->getID()].end(), task371Inputs);
    task371ListLock.unlock();
    invoke<TP371>(myTP, 1, this->getID(), myTP, task371Inputs, &(task371List[this->getID()]),
        &(task371ListLock), task371It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 371 nextRegion: 421 */
    myTP->controlTPParent->checkInCodelets421.decDep();
}
void TP319::_checkInCodelets421::fire(void)
{
    /*printing node 421: DeclStmt*/
    this->inputsTPParent->k_darts317 = 1;

    /*printing node 422: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets420.decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the end condional node.*/
        /*Signaling next codelet region: 423 nextRegion: 480 */
        myTP->controlTPParent->checkInCodelets480.decDep();
        return;
    }
}
void TP319::_checkInCodelets420::fire(void)
{

    /*printing node 420: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP420_LoopCounter),
        myTP->controlTPParent->TP420_LoopCounterPerThread[this->getLocalID()],
        myTP->controlTPParent->TP420_LoopCounterPerThread[this->getLocalID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP420_LoopCounterPerThread[this->getLocalID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP420PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP420_LoopCounterPerThread[this->getLocalID()] += 1;
        invoke<TP420>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP420PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP420PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP420PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP420PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP420PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP420PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP420_LoopCounterPerThread[this->getLocalID()] += 1;
        }
    }
}
void TP319::_checkInCodelets423::fire(void)
{

    /*printing node 423: UnaryOperator*/
    ++(this->inputsTPParent->k_darts317);

    /*printing node 813: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets420.decDep();
        return;
    } else {
        /*Signal the codelet after the loop from the condtional node.*/
        /*Signaling next codelet region: 423 nextRegion: 480 */
        myTP->controlTPParent->checkInCodelets480.decDep();
        return;
    }
}
void TP319::_checkInCodelets480::fire(void)
{
    /*printing node 480: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 480*/
    _task480Inputs* task480Inputs = new _task480Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->n_cols_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task480Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[0][(this->inputsTPParent->ts_darts317)]));
    task480Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[1][(this->inputsTPParent->ts_darts317)]));
    /*depend(out) vars*/
    task480Inputs->taskOutVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[0][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task480Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task480ListLock.lock();
    list<_task480Inputs*>::iterator task480It
        = task480List[this->getID()].insert(task480List[this->getID()].end(), task480Inputs);
    task480ListLock.unlock();
    invoke<TP480>(myTP, 1, this->getID(), myTP, task480Inputs, &(task480List[this->getID()]),
        &(task480ListLock), task480It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 480 nextRegion: 524 */
    myTP->controlTPParent->checkInCodelets524.decDep();
}
void TP319::_checkInCodelets524::fire(void)
{
    /*printing node 524: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 524*/
    _task524Inputs* task524Inputs = new _task524Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->nThrdsMinus1_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &((this->inputsTPParent->n_rows_darts317)),
        &(*(this->inputsTPParent->src_darts317)), &((this->inputsTPParent->ts_darts317)),
        nextSyncCodelet);
    /*depend(in) vars*/
    task524Inputs->taskInVarDependencies.push_back(&(
        ((this->inputsTPParent->dCalc_darts317))[(this->inputsTPParent->nThrdsMinus1_darts317) - 1]
                                                [(this->inputsTPParent->ts_darts317)]));
    task524Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*depend(out) vars*/
    task524Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->nThrdsMinus1_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task524Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task524ListLock.lock();
    list<_task524Inputs*>::iterator task524It
        = task524List[this->getID()].insert(task524List[this->getID()].end(), task524Inputs);
    task524ListLock.unlock();
    invoke<TP524>(myTP, 1, this->getID(), myTP, task524Inputs, &(task524List[this->getID()]),
        &(task524ListLock), task524It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 524 nextRegion: 572 */
    myTP->controlTPParent->checkInCodelets572.decDep();
}
void TP319::_checkInCodelets572::fire(void)
{
    /*printing node 572: DeclStmt*/
    this->inputsTPParent->k_darts317 = 1;

    /*printing node 573: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets571.decDep();
        return;
    }
    /*Signal the codelet after the loop from the end's condional node.*/
    /*The node is the last one in a complex loop, so signal the inc node*/
    myTP->controlTPParent->TPParent->checkInCodelets322.decDep();
}
void TP319::_checkInCodelets571::fire(void)
{

    /*printing node 571: ForStmt*/
    bool haveToLaunch = __sync_bool_compare_and_swap(&(myTP->controlTPParent->TP571_LoopCounter),
        myTP->controlTPParent->TP571_LoopCounterPerThread[this->getLocalID()],
        myTP->controlTPParent->TP571_LoopCounterPerThread[this->getLocalID()] + 1);
    unsigned int iterIdx = myTP->controlTPParent->TP571_LoopCounterPerThread[this->getLocalID()];
    if (haveToLaunch) {
        this->resetCodelet();
        myTP->controlTPParent->TP571PtrVec.push_back(nullptr);
        myTP->controlTPParent->TP571_LoopCounterPerThread[this->getLocalID()] += 1;
        invoke<TP571>(myTP, myTP->numThreads, this->getID(), myTP, myTP->inputsTPParent,
            &(myTP->controlTPParent->TP571PtrVec.back()));
    } else {
        if (myTP->controlTPParent->TP571PtrVec.size() == 0) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP571PtrVec.size() < (iterIdx + 1)) {
            this->resetCodelet();
            this->decDep();
            return;
        } else if (myTP->controlTPParent->TP571PtrVec[iterIdx] == nullptr) {
            this->resetCodelet();
            this->decDep();
            return;
        } else {
            this->resetCodelet();
#if USE_SPIN_CODELETS == 0
            myTP->controlTPParent->TP571PtrVec[iterIdx]->firstCodelet[this->getID()].decDep();
#else
            myTP->controlTPParent->TP571PtrVec[iterIdx]->availableCodelets[this->getID()] = 1;
#endif
            myTP->controlTPParent->TP571_LoopCounterPerThread[this->getLocalID()] += 1;
        }
    }
}
void TP319::_checkInCodelets574::fire(void)
{

    /*printing node 574: UnaryOperator*/
    ++(this->inputsTPParent->k_darts317);

    /*printing node 814: BinaryOperator*/
    /*Print the code for a condition node in a complex loop stmt */
    if ((this->inputsTPParent->k_darts317) < (this->inputsTPParent->nThrdsMinus1_darts317)) {
        this->resetCodelet();
        /*Signal the first codelet in the loop*/
        myTP->checkInCodelets571.decDep();
        return;
    }
    /*Signal the codelet after the loop from the condtional node.*/
    /*The node is the last one in a complex loop, so signal the inc node*/
    myTP->controlTPParent->TPParent->checkInCodelets322.decDep();
}
TP319::TP319(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent, TP317* in_inputsTPParent,
    TP319** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , TP420_LoopCounter(0)
    , TP420_LoopCounterPerThread(new unsigned int[this->numThreads])
    , TP571_LoopCounter(0)
    , TP571_LoopCounterPerThread(new unsigned int[this->numThreads])
    , checkInCodelets325(1, 1, this, this->mainCodeletID)
    , checkInCodelets371(1, 1, this, this->mainCodeletID)
    , checkInCodelets421(1, 1, this, this->mainCodeletID)
    , checkInCodelets420(1, 1, this, this->mainCodeletID)
    , checkInCodelets423(1, 1, this, this->mainCodeletID)
    , checkInCodelets480(1, 1, this, this->mainCodeletID)
    , checkInCodelets524(1, 1, this, this->mainCodeletID)
    , checkInCodelets572(1, 1, this, this->mainCodeletID)
    , checkInCodelets571(1, 1, this, this->mainCodeletID)
    , checkInCodelets574(1, 1, this, this->mainCodeletID)
{
    memset((void*)TP420_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    memset((void*)TP571_LoopCounterPerThread, 0, this->numThreads * sizeof(unsigned int));
    /*Initialize Codelets*/
    checkInCodelets574.setLocalID(0);
    checkInCodelets571.setLocalID(0);
    checkInCodelets572.setLocalID(0);
    checkInCodelets524.setLocalID(0);
    checkInCodelets480.setLocalID(0);
    checkInCodelets423.setLocalID(0);
    checkInCodelets420.setLocalID(0);
    checkInCodelets421.setLocalID(0);
    checkInCodelets371.setLocalID(0);
    checkInCodelets325.setLocalID(0);
    checkInCodelets325.decDep();
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP319::~TP319()
{
    delete[] TP420_LoopCounterPerThread;
    delete[] TP571_LoopCounterPerThread;
}
/*TP325: OMPTaskDirective*/
void TP325::_checkInCodelets326::fire(void)
{
    if (myTP->task325Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 326: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts325))[0][(*(this->taskInputs->ts_darts325))];

    /*printing node 328: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts325))[0][(*(this->taskInputs->ts_darts325)) + 2];

    /*printing node 331: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts325))[1][(*(this->taskInputs->ts_darts325)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 326 nextRegion: 335 */
    myTP->controlTPParent->checkInCodelets335.decDep();
}
void TP325::_checkInCodelets335::fire(void)
{
    /*printing node 335: DeclStmt*/
    this->taskInputs->DST_darts325 = (double*)(*(this->taskInputs->dst_darts325));
    this->taskInputs->SRC_darts325 = (double*)(*(this->taskInputs->src_darts325));

    /*printing node 338: DeclStmt*/
    this->taskInputs->pos1_darts325 = 1;

    /*printing node 339: DeclStmt*/
    this->taskInputs->pos2_darts325 = 1 + (this->taskInputs->chunk_darts325);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts325];
	
    /*printing node 341: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts325);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts325);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts325);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts325);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts325);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts325);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts325 = (this->taskInputs->pos1_darts325);
        size_t i_darts_counter_temp325 = (this->taskInputs->i_darts325);
        for (; (i_darts_counter_temp325) < (*pos2); ++(i_darts_counter_temp325)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp325 = (*j);
                for (; j_darts_counter_temp325 < (*n_cols) - 1; ++j_darts_counter_temp325) {
                    (*DST)[(i_darts_counter_temp325)][j_darts_counter_temp325]
                        = ((*SRC)[(i_darts_counter_temp325)-1][j_darts_counter_temp325]
                              + (*SRC)[(i_darts_counter_temp325) + 1][j_darts_counter_temp325]
                              + (*SRC)[(i_darts_counter_temp325)][j_darts_counter_temp325 - 1]
                              + (*SRC)[(i_darts_counter_temp325)][j_darts_counter_temp325 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp325;
            }
        }
        (this->taskInputs->i_darts325) = i_darts_counter_temp325;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task325Inputs->taskCompleted = true;
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask325ListLock->lock();
    myTP->controlTPParent->parentTask325List->erase(myTP->controlTPParent->parentTask325It);
    myTP->controlTPParent->parentTask325ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task325Inputs->nextSyncCodelet->decDep();
}
TP325::TP325(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
    _task325Inputs* in_task325Inputs, std::list<_task325Inputs*>* in_parentTask325List,
    std::mutex* in_parentTask325ListLock, std::list<_task325Inputs*>::iterator in_parentTask325It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task325Inputs(in_task325Inputs)
    , parentTask325List(in_parentTask325List)
    , parentTask325ListLock(in_parentTask325ListLock)
    , parentTask325It(in_parentTask325It)
    , checkInCodelets326(1, 1, this, this->mainCodeletID)
    , checkInCodelets335(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets335.setLocalID(0);
    checkInCodelets326.setLocalID(0);
    checkInCodelets326.decDep();
}
TP325::~TP325()
{
	task325ToDeleteVector.push_back(task325Inputs);
}
/*TP371: OMPTaskDirective*/
void TP371::_checkInCodelets372::fire(void)
{
    if (myTP->task371Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 372: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts371))[(this->taskInputs->nThrdsMinus1_darts371)]
                                        [(*(this->taskInputs->ts_darts371))];

    /*printing node 374: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts371))[(this->taskInputs->nThrdsMinus1_darts371) - 1]
                                        [(*(this->taskInputs->ts_darts371)) + 2];

    /*printing node 378: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts371))[(this->taskInputs->nThrdsMinus1_darts371)]
                                        [(*(this->taskInputs->ts_darts371)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 372 nextRegion: 382 */
    myTP->controlTPParent->checkInCodelets382.decDep();
}
void TP371::_checkInCodelets382::fire(void)
{

    /*printing node 382: DeclStmt*/
    this->taskInputs->DST_darts371 = (double*)(*(this->taskInputs->dst_darts371));
    this->taskInputs->SRC_darts371 = (double*)(*(this->taskInputs->src_darts371));

    /*printing node 385: DeclStmt*/
    this->taskInputs->pos1_darts371
        = 1 + (this->taskInputs->chunk_darts371) * (this->taskInputs->nThrdsMinus1_darts371);

    /*printing node 388: DeclStmt*/
    this->taskInputs->pos2_darts371 = (this->taskInputs->n_rows_darts371) - 1;

	typedef double (*Array2D)[this->taskInputs->n_cols_darts371];
	
    /*printing node 390: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts371);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts371);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts371);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts371);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts371);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts371);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts371 = (this->taskInputs->pos1_darts371);
        size_t i_darts_counter_temp371 = (this->taskInputs->i_darts371);
        for (; (i_darts_counter_temp371) < (*pos2); ++(i_darts_counter_temp371)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp371 = (*j);
                for (; j_darts_counter_temp371 < (*n_cols) - 1; ++j_darts_counter_temp371) {
                    (*DST)[(i_darts_counter_temp371)][j_darts_counter_temp371]
                        = ((*SRC)[(i_darts_counter_temp371)-1][j_darts_counter_temp371]
                              + (*SRC)[(i_darts_counter_temp371) + 1][j_darts_counter_temp371]
                              + (*SRC)[(i_darts_counter_temp371)][j_darts_counter_temp371 - 1]
                              + (*SRC)[(i_darts_counter_temp371)][j_darts_counter_temp371 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp371;
            }
        }
        (this->taskInputs->i_darts371) = i_darts_counter_temp371;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task371Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask371ListLock->lock();
    myTP->controlTPParent->parentTask371List->erase(myTP->controlTPParent->parentTask371It);
    myTP->controlTPParent->parentTask371ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task371Inputs->nextSyncCodelet->decDep();
}
TP371::TP371(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
    _task371Inputs* in_task371Inputs, std::list<_task371Inputs*>* in_parentTask371List,
    std::mutex* in_parentTask371ListLock, std::list<_task371Inputs*>::iterator in_parentTask371It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task371Inputs(in_task371Inputs)
    , parentTask371List(in_parentTask371List)
    , parentTask371ListLock(in_parentTask371ListLock)
    , parentTask371It(in_parentTask371It)
    , checkInCodelets372(1, 1, this, this->mainCodeletID)
    , checkInCodelets382(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets382.setLocalID(0);
    checkInCodelets372.setLocalID(0);
    checkInCodelets372.decDep();
}
TP371::~TP371()
{
	task371ToDeleteVector.push_back(task371Inputs);
}
/*TP420: ForStmt*/
void TP420::_checkInCodelets425::fire(void)
{

    /*printing node 425: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet
        = &(myTP->controlTPParent->TPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 425*/
    _task425Inputs* task425Inputs = new _task425Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->k_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &(*(this->inputsTPParent->src_darts317)),
        &((this->inputsTPParent->ts_darts317)), nextSyncCodelet);
    /*depend(in) vars*/
    task425Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->k_darts317) - 1]
                                                  [(this->inputsTPParent->ts_darts317) + 2]));
    task425Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317) + 2]));
    task425Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->k_darts317) + 1]
                                                  [(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task425Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task425Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task425ListLock.lock();
    list<_task425Inputs*>::iterator task425It
        = task425List[this->getID()].insert(task425List[this->getID()].end(), task425Inputs);
    task425ListLock.unlock();
    invoke<TP425>(myTP, 1, this->getID(), myTP, task425Inputs, &(task425List[this->getID()]),
        &(task425ListLock), task425It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*The node is the last one in a complex loop, so signal the inc node*/
    myTP->controlTPParent->TPParent->checkInCodelets423.decDep();
}
TP420::TP420(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent, TP317* in_inputsTPParent,
    TP420** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , checkInCodelets425(1, 1, this, this->mainCodeletID)
{
    /*Initialize Codelets*/
    checkInCodelets425.setLocalID(0);
    checkInCodelets425.decDep();
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP420::~TP420() {}
/*TP425: OMPTaskDirective*/
void TP425::_checkInCodelets426::fire(void)
{
    if (myTP->task425Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 426: ArraySubscriptExpr*/
    ((this->taskInputs
            ->dCalc_darts425))[(this->taskInputs->k_darts425)][(*(this->taskInputs->ts_darts425))];

    /*printing node 428: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts425))[(this->taskInputs->k_darts425) - 1]
                                        [(*(this->taskInputs->ts_darts425)) + 2];

    /*printing node 432: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts425))[(this->taskInputs->k_darts425)]
                                        [(*(this->taskInputs->ts_darts425)) + 2];

    /*printing node 435: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts425))[(this->taskInputs->k_darts425) + 1]
                                        [(*(this->taskInputs->ts_darts425)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 426 nextRegion: 440 */
    myTP->controlTPParent->checkInCodelets440.decDep();
}
void TP425::_checkInCodelets440::fire(void)
{
    /*printing node 440: DeclStmt*/
    this->taskInputs->DST_darts425 = (double*)(*(this->taskInputs->dst_darts425));
    this->taskInputs->SRC_darts425 = (double*)(*(this->taskInputs->src_darts425));

    /*printing node 443: DeclStmt*/
    this->taskInputs->pos1_darts425
        = 1 + (this->taskInputs->chunk_darts425) * (this->taskInputs->k_darts425);

    /*printing node 446: DeclStmt*/
    this->taskInputs->pos2_darts425
        = 1 + (this->taskInputs->chunk_darts425) * ((this->taskInputs->k_darts425) + 1);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts425];
		
    /*printing node 450: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts425);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts425);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts425);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts425);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts425);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts425);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts425 = (this->taskInputs->pos1_darts425);
        size_t i_darts_counter_temp425 = (this->taskInputs->i_darts425);
        for (; (i_darts_counter_temp425) < (*pos2); ++(i_darts_counter_temp425)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp425 = (*j);
                for (; j_darts_counter_temp425 < (*n_cols) - 1; ++j_darts_counter_temp425) {
                    (*DST)[(i_darts_counter_temp425)][j_darts_counter_temp425]
                        = ((*SRC)[(i_darts_counter_temp425)-1][j_darts_counter_temp425]
                              + (*SRC)[(i_darts_counter_temp425) + 1][j_darts_counter_temp425]
                              + (*SRC)[(i_darts_counter_temp425)][j_darts_counter_temp425 - 1]
                              + (*SRC)[(i_darts_counter_temp425)][j_darts_counter_temp425 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp425;
            }
        }
        (this->taskInputs->i_darts425) = i_darts_counter_temp425;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task425Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask425ListLock->lock();
    myTP->controlTPParent->parentTask425List->erase(myTP->controlTPParent->parentTask425It);
    myTP->controlTPParent->parentTask425ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task425Inputs->nextSyncCodelet->decDep();
}
TP425::TP425(int in_numThreads, int in_mainCodeletID, TP420* in_TPParent,
    _task425Inputs* in_task425Inputs, std::list<_task425Inputs*>* in_parentTask425List,
    std::mutex* in_parentTask425ListLock, std::list<_task425Inputs*>::iterator in_parentTask425It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task425Inputs(in_task425Inputs)
    , parentTask425List(in_parentTask425List)
    , parentTask425ListLock(in_parentTask425ListLock)
    , parentTask425It(in_parentTask425It)
    , checkInCodelets426(1, 1, this, this->mainCodeletID)
    , checkInCodelets440(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets440.setLocalID(0);
    checkInCodelets426.setLocalID(0);
    checkInCodelets426.decDep();
}
TP425::~TP425()
{
	task425ToDeleteVector.push_back(task425Inputs);
}
/*TP480: OMPTaskDirective*/
void TP480::_checkInCodelets481::fire(void)
{
    if (myTP->task480Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 481: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts480))[0][(*(this->taskInputs->ts_darts480))];

    /*printing node 483: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts480))[0][(*(this->taskInputs->ts_darts480))];

    /*printing node 485: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts480))[1][(*(this->taskInputs->ts_darts480))];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 481 nextRegion: 488 */
    myTP->controlTPParent->checkInCodelets488.decDep();
}
void TP480::_checkInCodelets488::fire(void)
{
    /*printing node 488: DeclStmt*/
    this->taskInputs->DST_darts480 = (double*)(*(this->taskInputs->dst_darts480));
    this->taskInputs->SRC_darts480 = (double*)(*(this->taskInputs->src_darts480));

    /*printing node 491: DeclStmt*/
    this->taskInputs->pos1_darts480 = 1;

    /*printing node 492: DeclStmt*/
    this->taskInputs->pos2_darts480 = 1 + (this->taskInputs->chunk_darts480);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts480];
	
    /*printing node 494: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts480);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts480);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts480);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts480);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts480);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts480);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts480 = (this->taskInputs->pos1_darts480);
        size_t i_darts_counter_temp480 = (this->taskInputs->i_darts480);
        for (; (i_darts_counter_temp480) < (*pos2); ++(i_darts_counter_temp480)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp480 = (*j);
                for (; j_darts_counter_temp480 < (*n_cols) - 1; ++j_darts_counter_temp480) {
                    (*SRC)[(i_darts_counter_temp480)][j_darts_counter_temp480]
                        = ((*DST)[(i_darts_counter_temp480)-1][j_darts_counter_temp480]
                              + (*DST)[(i_darts_counter_temp480) + 1][j_darts_counter_temp480]
                              + (*DST)[(i_darts_counter_temp480)][j_darts_counter_temp480 - 1]
                              + (*DST)[(i_darts_counter_temp480)][j_darts_counter_temp480 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp480;
            }
        }
        (this->taskInputs->i_darts480) = i_darts_counter_temp480;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task480Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask480ListLock->lock();
    myTP->controlTPParent->parentTask480List->erase(myTP->controlTPParent->parentTask480It);
    myTP->controlTPParent->parentTask480ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task480Inputs->nextSyncCodelet->decDep();
}
TP480::TP480(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
    _task480Inputs* in_task480Inputs, std::list<_task480Inputs*>* in_parentTask480List,
    std::mutex* in_parentTask480ListLock, std::list<_task480Inputs*>::iterator in_parentTask480It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task480Inputs(in_task480Inputs)
    , parentTask480List(in_parentTask480List)
    , parentTask480ListLock(in_parentTask480ListLock)
    , parentTask480It(in_parentTask480It)
    , checkInCodelets481(1, 1, this, this->mainCodeletID)
    , checkInCodelets488(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets488.setLocalID(0);
    checkInCodelets481.setLocalID(0);
    checkInCodelets481.decDep();
}
TP480::~TP480()
{
	task480ToDeleteVector.push_back(task480Inputs);
}
/*TP524: OMPTaskDirective*/
void TP524::_checkInCodelets525::fire(void)
{
    if (myTP->task524Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 525: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts524))[(this->taskInputs->nThrdsMinus1_darts524)]
                                        [(*(this->taskInputs->ts_darts524))];

    /*printing node 527: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts524))[(this->taskInputs->nThrdsMinus1_darts524) - 1]
                                        [(*(this->taskInputs->ts_darts524))];

    /*printing node 530: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts524))[(this->taskInputs->nThrdsMinus1_darts524)]
                                        [(*(this->taskInputs->ts_darts524))];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 525 nextRegion: 533 */
    myTP->controlTPParent->checkInCodelets533.decDep();
}
void TP524::_checkInCodelets533::fire(void)
{
    /*printing node 533: DeclStmt*/
    this->taskInputs->DST_darts524 = (double*)(*(this->taskInputs->dst_darts524));
    this->taskInputs->SRC_darts524 = (double*)(*(this->taskInputs->src_darts524));

    /*printing node 536: DeclStmt*/
    this->taskInputs->pos1_darts524
        = 1 + (this->taskInputs->chunk_darts524) * (this->taskInputs->nThrdsMinus1_darts524);

    /*printing node 539: DeclStmt*/
    this->taskInputs->pos2_darts524 = (this->taskInputs->n_rows_darts524) - 1;

	typedef double (*Array2D)[this->taskInputs->n_cols_darts524];
	
    /*printing node 541: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts524);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts524);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts524);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts524);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts524);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts524);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts524 = (this->taskInputs->pos1_darts524);
        size_t i_darts_counter_temp524 = (this->taskInputs->i_darts524);
        for (; (i_darts_counter_temp524) < (*pos2); ++(i_darts_counter_temp524)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp524 = (*j);
                for (; j_darts_counter_temp524 < (*n_cols) - 1; ++j_darts_counter_temp524) {
                    (*SRC)[(i_darts_counter_temp524)][j_darts_counter_temp524]
                        = ((*DST)[(i_darts_counter_temp524)-1][j_darts_counter_temp524]
                              + (*DST)[(i_darts_counter_temp524) + 1][j_darts_counter_temp524]
                              + (*DST)[(i_darts_counter_temp524)][j_darts_counter_temp524 - 1]
                              + (*DST)[(i_darts_counter_temp524)][j_darts_counter_temp524 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp524;
            }
        }
        (this->taskInputs->i_darts524) = i_darts_counter_temp524;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task524Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask524ListLock->lock();
    myTP->controlTPParent->parentTask524List->erase(myTP->controlTPParent->parentTask524It);
    myTP->controlTPParent->parentTask524ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task524Inputs->nextSyncCodelet->decDep();
}
TP524::TP524(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent,
    _task524Inputs* in_task524Inputs, std::list<_task524Inputs*>* in_parentTask524List,
    std::mutex* in_parentTask524ListLock, std::list<_task524Inputs*>::iterator in_parentTask524It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task524Inputs(in_task524Inputs)
    , parentTask524List(in_parentTask524List)
    , parentTask524ListLock(in_parentTask524ListLock)
    , parentTask524It(in_parentTask524It)
    , checkInCodelets525(1, 1, this, this->mainCodeletID)
    , checkInCodelets533(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets533.setLocalID(0);
    checkInCodelets525.setLocalID(0);
    checkInCodelets525.decDep();
}
TP524::~TP524()
{
	task524ToDeleteVector.push_back(task524Inputs);
}
/*TP571: ForStmt*/
void TP571::_checkInCodelets576::fire(void)
{

    /*printing node 576: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet
        = &(myTP->controlTPParent->TPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 576*/
    _task576Inputs* task576Inputs = new _task576Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->k_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &(*(this->inputsTPParent->src_darts317)),
        &((this->inputsTPParent->ts_darts317)), nextSyncCodelet);
    /*depend(in) vars*/
    task576Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[(this->inputsTPParent->k_darts317) - 1]
                                                  [(this->inputsTPParent->ts_darts317)]));
    task576Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317)]));
    task576Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dCalc_darts317))[(this->inputsTPParent->k_darts317) + 1]
                                                  [(this->inputsTPParent->ts_darts317)]));
    /*depend(out) vars*/
    task576Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task576Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task576ListLock.lock();
    list<_task576Inputs*>::iterator task576It
        = task576List[this->getID()].insert(task576List[this->getID()].end(), task576Inputs);
    task576ListLock.unlock();
    invoke<TP576>(myTP, 1, this->getID(), myTP, task576Inputs, &(task576List[this->getID()]),
        &(task576ListLock), task576It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*The node is the last one in a complex loop, so signal the inc node*/
    myTP->controlTPParent->TPParent->checkInCodelets574.decDep();
}
TP571::TP571(int in_numThreads, int in_mainCodeletID, TP319* in_TPParent, TP317* in_inputsTPParent,
    TP571** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , checkInCodelets576(1, 1, this, this->mainCodeletID)
{
    /*Initialize Codelets*/
    checkInCodelets576.setLocalID(0);
    checkInCodelets576.decDep();
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP571::~TP571() {}
/*TP576: OMPTaskDirective*/
void TP576::_checkInCodelets577::fire(void)
{
    if (myTP->task576Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 577: ArraySubscriptExpr*/
    ((this->taskInputs
            ->dSwap_darts576))[(this->taskInputs->k_darts576)][(*(this->taskInputs->ts_darts576))];

    /*printing node 579: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts576))[(this->taskInputs->k_darts576) - 1]
                                        [(*(this->taskInputs->ts_darts576))];

    /*printing node 582: ArraySubscriptExpr*/
    ((this->taskInputs
            ->dCalc_darts576))[(this->taskInputs->k_darts576)][(*(this->taskInputs->ts_darts576))];

    /*printing node 584: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts576))[(this->taskInputs->k_darts576) + 1]
                                        [(*(this->taskInputs->ts_darts576))];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 577 nextRegion: 588 */
    myTP->controlTPParent->checkInCodelets588.decDep();
}
void TP576::_checkInCodelets588::fire(void)
{
    /*printing node 588: DeclStmt*/
    this->taskInputs->DST_darts576 = (double*)(*(this->taskInputs->dst_darts576));
    this->taskInputs->SRC_darts576 = (double*)(*(this->taskInputs->src_darts576));

    /*printing node 591: DeclStmt*/
    this->taskInputs->pos1_darts576
        = 1 + (this->taskInputs->chunk_darts576) * (this->taskInputs->k_darts576);

    /*printing node 594: DeclStmt*/
    this->taskInputs->pos2_darts576
        = 1 + (this->taskInputs->chunk_darts576) * ((this->taskInputs->k_darts576) + 1);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts576];
		
    /*printing node 598: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts576);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts576);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts576);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts576);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts576);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts576);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts576 = (this->taskInputs->pos1_darts576);
        size_t i_darts_counter_temp576 = (this->taskInputs->i_darts576);
        for (; (i_darts_counter_temp576) < (*pos2); ++(i_darts_counter_temp576)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp576 = (*j);
                for (; j_darts_counter_temp576 < (*n_cols) - 1; ++j_darts_counter_temp576) {
                    (*SRC)[(i_darts_counter_temp576)][j_darts_counter_temp576]
                        = ((*DST)[(i_darts_counter_temp576)-1][j_darts_counter_temp576]
                              + (*DST)[(i_darts_counter_temp576) + 1][j_darts_counter_temp576]
                              + (*DST)[(i_darts_counter_temp576)][j_darts_counter_temp576 - 1]
                              + (*DST)[(i_darts_counter_temp576)][j_darts_counter_temp576 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp576;
            }
        }
        (this->taskInputs->i_darts576) = i_darts_counter_temp576;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task576Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask576ListLock->lock();
    myTP->controlTPParent->parentTask576List->erase(myTP->controlTPParent->parentTask576It);
    myTP->controlTPParent->parentTask576ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task576Inputs->nextSyncCodelet->decDep();
}
TP576::TP576(int in_numThreads, int in_mainCodeletID, TP571* in_TPParent,
    _task576Inputs* in_task576Inputs, std::list<_task576Inputs*>* in_parentTask576List,
    std::mutex* in_parentTask576ListLock, std::list<_task576Inputs*>::iterator in_parentTask576It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task576Inputs(in_task576Inputs)
    , parentTask576List(in_parentTask576List)
    , parentTask576ListLock(in_parentTask576ListLock)
    , parentTask576It(in_parentTask576It)
    , checkInCodelets577(1, 1, this, this->mainCodeletID)
    , checkInCodelets588(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets588.setLocalID(0);
    checkInCodelets577.setLocalID(0);
    checkInCodelets577.decDep();
}
TP576::~TP576()
{
	task576ToDeleteVector.push_back(task576Inputs);
}
/*TP631: OMPTaskDirective*/
void TP631::_checkInCodelets632::fire(void)
{
    if (myTP->task631Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 632: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts631))[0][(*(this->taskInputs->ts_darts631))];

    /*printing node 634: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts631))[0][(*(this->taskInputs->ts_darts631)) + 2];

    /*printing node 637: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts631))[1][(*(this->taskInputs->ts_darts631)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 632 nextRegion: 641 */
    myTP->controlTPParent->checkInCodelets641.decDep();
}
void TP631::_checkInCodelets641::fire(void)
{

    /*printing node 641: DeclStmt*/
    this->taskInputs->DST_darts631 = (double*)(*(this->taskInputs->dst_darts631));
    this->taskInputs->SRC_darts631 = (double*)(*(this->taskInputs->src_darts631));

    /*printing node 644: DeclStmt*/
    this->taskInputs->pos1_darts631 = 1;

    /*printing node 645: DeclStmt*/
    this->taskInputs->pos2_darts631 = 1 + (this->taskInputs->chunk_darts631);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts631];
	
    /*printing node 647: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts631);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts631);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts631);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts631);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts631);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts631);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts631 = (this->taskInputs->pos1_darts631);
        size_t i_darts_counter_temp631 = (this->taskInputs->i_darts631);
        for (; (i_darts_counter_temp631) < (*pos2); ++(i_darts_counter_temp631)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp631 = (*j);
                for (; j_darts_counter_temp631 < (*n_cols) - 1; ++j_darts_counter_temp631) {
                    (*DST)[(i_darts_counter_temp631)][j_darts_counter_temp631]
                        = ((*SRC)[(i_darts_counter_temp631)-1][j_darts_counter_temp631]
                              + (*SRC)[(i_darts_counter_temp631) + 1][j_darts_counter_temp631]
                              + (*SRC)[(i_darts_counter_temp631)][j_darts_counter_temp631 - 1]
                              + (*SRC)[(i_darts_counter_temp631)][j_darts_counter_temp631 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp631;
            }
        }
        (this->taskInputs->i_darts631) = i_darts_counter_temp631;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task631Inputs->taskCompleted = true;
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask631ListLock->lock();
    myTP->controlTPParent->parentTask631List->erase(myTP->controlTPParent->parentTask631It);
    myTP->controlTPParent->parentTask631ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task631Inputs->nextSyncCodelet->decDep();
}
TP631::TP631(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent,
    _task631Inputs* in_task631Inputs, std::list<_task631Inputs*>* in_parentTask631List,
    std::mutex* in_parentTask631ListLock, std::list<_task631Inputs*>::iterator in_parentTask631It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task631Inputs(in_task631Inputs)
    , parentTask631List(in_parentTask631List)
    , parentTask631ListLock(in_parentTask631ListLock)
    , parentTask631It(in_parentTask631It)
    , checkInCodelets632(1, 1, this, this->mainCodeletID)
    , checkInCodelets641(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets641.setLocalID(0);
    checkInCodelets632.setLocalID(0);
    checkInCodelets632.decDep();
}
TP631::~TP631()
{ 
	task631ToDeleteVector.push_back(task631Inputs);
}
/*TP677: OMPTaskDirective*/
void TP677::_checkInCodelets678::fire(void)
{
    if (myTP->task677Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 678: ArraySubscriptExpr*/
    ((this->taskInputs->dCalc_darts677))[(this->taskInputs->nThrdsMinus1_darts677)]
                                        [(*(this->taskInputs->ts_darts677))];

    /*printing node 680: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts677))[(this->taskInputs->nThrdsMinus1_darts677) - 1]
                                        [(*(this->taskInputs->ts_darts677)) + 2];

    /*printing node 684: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts677))[(this->taskInputs->nThrdsMinus1_darts677)]
                                        [(*(this->taskInputs->ts_darts677)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 678 nextRegion: 688 */
    myTP->controlTPParent->checkInCodelets688.decDep();
}
void TP677::_checkInCodelets688::fire(void)
{
    /*printing node 688: DeclStmt*/
    this->taskInputs->DST_darts677 = (double*)(*(this->taskInputs->dst_darts677));
    this->taskInputs->SRC_darts677 = (double*)(*(this->taskInputs->src_darts677));

    /*printing node 691: DeclStmt*/
    this->taskInputs->pos1_darts677
        = 1 + (this->taskInputs->chunk_darts677) * (this->taskInputs->nThrdsMinus1_darts677);

    /*printing node 694: DeclStmt*/
    this->taskInputs->pos2_darts677 = (this->taskInputs->n_rows_darts677) - 1;

	typedef double (*Array2D)[this->taskInputs->n_cols_darts677];
	
    /*printing node 696: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts677);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts677);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts677);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts677);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts677);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts677);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts677 = (this->taskInputs->pos1_darts677);
        size_t i_darts_counter_temp677 = (this->taskInputs->i_darts677);
        for (; (i_darts_counter_temp677) < (*pos2); ++(i_darts_counter_temp677)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp677 = (*j);
                for (; j_darts_counter_temp677 < (*n_cols) - 1; ++j_darts_counter_temp677) {
                    (*DST)[(i_darts_counter_temp677)][j_darts_counter_temp677]
                        = ((*SRC)[(i_darts_counter_temp677)-1][j_darts_counter_temp677]
                              + (*SRC)[(i_darts_counter_temp677) + 1][j_darts_counter_temp677]
                              + (*SRC)[(i_darts_counter_temp677)][j_darts_counter_temp677 - 1]
                              + (*SRC)[(i_darts_counter_temp677)][j_darts_counter_temp677 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp677;
            }
        }
        (this->taskInputs->i_darts677) = i_darts_counter_temp677;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task677Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask677ListLock->lock();
    myTP->controlTPParent->parentTask677List->erase(myTP->controlTPParent->parentTask677It);
    myTP->controlTPParent->parentTask677ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task677Inputs->nextSyncCodelet->decDep();
}
TP677::TP677(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent,
    _task677Inputs* in_task677Inputs, std::list<_task677Inputs*>* in_parentTask677List,
    std::mutex* in_parentTask677ListLock, std::list<_task677Inputs*>::iterator in_parentTask677It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task677Inputs(in_task677Inputs)
    , parentTask677List(in_parentTask677List)
    , parentTask677ListLock(in_parentTask677ListLock)
    , parentTask677It(in_parentTask677It)
    , checkInCodelets678(1, 1, this, this->mainCodeletID)
    , checkInCodelets688(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets688.setLocalID(0);
    checkInCodelets678.setLocalID(0);
    checkInCodelets678.decDep();
}
TP677::~TP677()
{
	task677ToDeleteVector.push_back(task677Inputs);
}
/*TP726: ForStmt*/
void TP726::_checkInCodelets731::fire(void)
{

    /*printing node 731: OMPTaskDirective*/
    /*syncNode: OMPTaskgroupDirective 317*/
    Codelet* nextSyncCodelet = &(myTP->controlTPParent->TPParent->TPParent->barrierCodelets317[0]);
    /*Increment sync point's dependency to account for the task to be launched*/
    nextSyncCodelet->incDep();
    /*Encapsulating data for task 731*/
    _task731Inputs* task731Inputs = new _task731Inputs(&((this->inputsTPParent->chunk_darts317)),
        ((this->inputsTPParent->dCalc_darts317)), (this->inputsTPParent->dCalc_outer317_size),
        ((this->inputsTPParent->dSwap_darts317)), (this->inputsTPParent->dSwap_outer317_size),
        &(*(this->inputsTPParent->dst_darts317)), &((this->inputsTPParent->k_darts317)),
        &((this->inputsTPParent->n_cols_darts317)), &(*(this->inputsTPParent->src_darts317)),
        &((this->inputsTPParent->ts_darts317)), nextSyncCodelet);
    /*depend(in) vars*/
    task731Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->k_darts317) - 1]
                                                  [(this->inputsTPParent->ts_darts317) + 2]));
    task731Inputs->taskInVarDependencies.push_back(&(((this->inputsTPParent->dSwap_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317) + 2]));
    task731Inputs->taskInVarDependencies.push_back(
        &(((this->inputsTPParent->dSwap_darts317))[(this->inputsTPParent->k_darts317) + 1]
                                                  [(this->inputsTPParent->ts_darts317) + 2]));
    /*depend(out) vars*/
    task731Inputs->taskOutVarDependencies.push_back(&(((this->inputsTPParent->dCalc_darts317))[(
        this->inputsTPParent->k_darts317)][(this->inputsTPParent->ts_darts317)]));
    /*Dependencies wrt task325*/
    task325ListLock.lock();
    for (_task325Inputs* tempTaskInput : task325List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task325ListLock.unlock();
    /*Dependencies wrt task371*/
    task371ListLock.lock();
    for (_task371Inputs* tempTaskInput : task371List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task371ListLock.unlock();
    /*Dependencies wrt task425*/
    task425ListLock.lock();
    for (_task425Inputs* tempTaskInput : task425List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task425ListLock.unlock();
    /*Dependencies wrt task480*/
    task480ListLock.lock();
    for (_task480Inputs* tempTaskInput : task480List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task480ListLock.unlock();
    /*Dependencies wrt task524*/
    task524ListLock.lock();
    for (_task524Inputs* tempTaskInput : task524List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task524ListLock.unlock();
    /*Dependencies wrt task576*/
    task576ListLock.lock();
    for (_task576Inputs* tempTaskInput : task576List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task576ListLock.unlock();
    /*Dependencies wrt task631*/
    task631ListLock.lock();
    for (_task631Inputs* tempTaskInput : task631List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task631ListLock.unlock();
    /*Dependencies wrt task677*/
    task677ListLock.lock();
    for (_task677Inputs* tempTaskInput : task677List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task677ListLock.unlock();
    /*Dependencies wrt task731*/
    task731ListLock.lock();
    for (_task731Inputs* tempTaskInput : task731List[this->getID()])
        task731Inputs->setDeps(tempTaskInput);
    task731ListLock.unlock();
    /*Save the pointer to the recently created task's data*/
    task731ListLock.lock();
    list<_task731Inputs*>::iterator task731It
        = task731List[this->getID()].insert(task731List[this->getID()].end(), task731Inputs);
    task731ListLock.unlock();
    invoke<TP731>(myTP, 1, this->getID(), myTP, task731Inputs, &(task731List[this->getID()]),
        &(task731ListLock), task731It);
    /*Signaling next codelet from last stmt in the codelet*/
    /*The node is the last one in a complex loop, so signal the inc node*/
    myTP->controlTPParent->TPParent->checkInCodelets729.decDep();
}
TP726::TP726(int in_numThreads, int in_mainCodeletID, TP317* in_TPParent, TP317* in_inputsTPParent,
    TP726** in_ptrToThisTP)
    : ompTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(in_inputsTPParent)
    , ptrToThisTP(in_ptrToThisTP)
    , checkInCodelets731(1, 1, this, this->mainCodeletID)
{
    /*Initialize Codelets*/
    checkInCodelets731.setLocalID(0);
    checkInCodelets731.decDep();
    *(this->ptrToThisTP) = this;
#if USE_SPIN_CODELETS == 0
    this->firstCodelet[this->getID()].decDep();
#else
    this->availableCodelets[this->getID()] = 1;
#endif
}
TP726::~TP726() {}
/*TP731: OMPTaskDirective*/
void TP731::_checkInCodelets732::fire(void)
{
    if (myTP->task731Inputs->checkForDeps() == false) {
        myTP->add(this);
        return;
    }
    /*Init the vars for this region*/

    /*printing node 732: ArraySubscriptExpr*/
    ((this->taskInputs
            ->dCalc_darts731))[(this->taskInputs->k_darts731)][(*(this->taskInputs->ts_darts731))];

    /*printing node 734: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts731))[(this->taskInputs->k_darts731) - 1]
                                        [(*(this->taskInputs->ts_darts731)) + 2];

    /*printing node 738: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts731))[(this->taskInputs->k_darts731)]
                                        [(*(this->taskInputs->ts_darts731)) + 2];

    /*printing node 741: ArraySubscriptExpr*/
    ((this->taskInputs->dSwap_darts731))[(this->taskInputs->k_darts731) + 1]
                                        [(*(this->taskInputs->ts_darts731)) + 2];
    /*Signaling next codelet from last stmt in the codelet*/
    /*Signaling next codelet region: 732 nextRegion: 746 */
    myTP->controlTPParent->checkInCodelets746.decDep();
}
void TP731::_checkInCodelets746::fire(void)
{
    /*printing node 746: DeclStmt*/
    this->taskInputs->DST_darts731 = (double*)(*(this->taskInputs->dst_darts731));
    this->taskInputs->SRC_darts731 = (double*)(*(this->taskInputs->src_darts731));

    /*printing node 749: DeclStmt*/
    this->taskInputs->pos1_darts731
        = 1 + (this->taskInputs->chunk_darts731) * (this->taskInputs->k_darts731);

    /*printing node 752: DeclStmt*/
    this->taskInputs->pos2_darts731
        = 1 + (this->taskInputs->chunk_darts731) * ((this->taskInputs->k_darts731) + 1);

	typedef double (*Array2D)[this->taskInputs->n_cols_darts731];
		
    /*printing node 756: ForStmt*/
    {
        size_t* n_cols = &(this->taskInputs->n_cols_darts731);
        (void)n_cols /*OMP_FIRSTPRIVATE*/;
        Array2D* DST = (Array2D*)&(this->taskInputs->DST_darts731);
        (void)DST /*PRIVATE*/;
        Array2D* SRC = (Array2D*)&(this->taskInputs->SRC_darts731);
        (void)SRC /*PRIVATE*/;
        size_t* i = &(this->taskInputs->i_darts731);
        (void)i /*PRIVATE*/;
        size_t* j = &(this->taskInputs->j_darts731);
        (void)j /*PRIVATE*/;
        size_t* pos2 = &(this->taskInputs->pos2_darts731);
        (void)pos2 /*PRIVATE*/;
        /*Loop's init*/
        this->taskInputs->i_darts731 = (this->taskInputs->pos1_darts731);
        size_t i_darts_counter_temp731 = (this->taskInputs->i_darts731);
        for (; (i_darts_counter_temp731) < (*pos2); ++(i_darts_counter_temp731)) {
            {
                /*Loop's init*/
                *j = 1;
                size_t j_darts_counter_temp731 = (*j);
                for (; j_darts_counter_temp731 < (*n_cols) - 1; ++j_darts_counter_temp731) {
                    (*DST)[(i_darts_counter_temp731)][j_darts_counter_temp731]
                        = ((*SRC)[(i_darts_counter_temp731)-1][j_darts_counter_temp731]
                              + (*SRC)[(i_darts_counter_temp731) + 1][j_darts_counter_temp731]
                              + (*SRC)[(i_darts_counter_temp731)][j_darts_counter_temp731 - 1]
                              + (*SRC)[(i_darts_counter_temp731)][j_darts_counter_temp731 + 1])
                        / 4;
                }
                (*j) = j_darts_counter_temp731;
            }
        }
        (this->taskInputs->i_darts731) = i_darts_counter_temp731;
    }
    /*Signaling next codelet from last stmt in the codelet*/
    /*Mark this task as completed*/
    myTP->controlTPParent->task731Inputs->taskCompleted = true;	
    /*Remove this task from the pool*/
    myTP->controlTPParent->parentTask731ListLock->lock();
    myTP->controlTPParent->parentTask731List->erase(myTP->controlTPParent->parentTask731It);
    myTP->controlTPParent->parentTask731ListLock->unlock();
    /*Signal the task's synchronization point*/
    myTP->controlTPParent->task731Inputs->nextSyncCodelet->decDep();
}
TP731::TP731(int in_numThreads, int in_mainCodeletID, TP726* in_TPParent,
    _task731Inputs* in_task731Inputs, std::list<_task731Inputs*>* in_parentTask731List,
    std::mutex* in_parentTask731ListLock, std::list<_task731Inputs*>::iterator in_parentTask731It)
    : ompOMPTaskDirectiveTP(in_numThreads, in_mainCodeletID)
    , TPParent(in_TPParent)
    , controlTPParent(this)
    , inputsTPParent(this)
    , task731Inputs(in_task731Inputs)
    , parentTask731List(in_parentTask731List)
    , parentTask731ListLock(in_parentTask731ListLock)
    , parentTask731It(in_parentTask731It)
    , checkInCodelets732(1, 1, this, this->mainCodeletID)
    , checkInCodelets746(1, 1, this, this->mainCodeletID)
{
    /*Initialize inputs and vars.*/
    /*Initialize Codelets*/
    checkInCodelets746.setLocalID(0);
    checkInCodelets732.setLocalID(0);
    checkInCodelets732.decDep();
}
TP731::~TP731()
{
	task731ToDeleteVector.push_back(task731Inputs);
}
