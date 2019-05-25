#include "timer.output.darts.h"
using namespace darts;
using namespace std;
/*Function: timer, ID: 17*/
void* timer(void* arg)
{
    /*timer:17*/
    /*CompoundStmt:37*/
    loc_timer_t* t = (loc_timer_t*)arg;
    cpu_set_t cpuset;
    do
        __builtin_memset(&cpuset, '\x00', sizeof(cpu_set_t));
    while (0);
    __extension__({
        size_t __cpu = (t->cpu);
        __cpu < 8 * (sizeof(cpu_set_t))
            ? (((__cpu_mask*)((&cpuset)->__bits))[((__cpu) / (8 * sizeof(__cpu_mask)))]
                  |= ((__cpu_mask)1 << ((__cpu) % (8 * sizeof(__cpu_mask)))))
            : 0;
    });
    spthread_setaffinity_np(pthread_self(), 1UL, &cpuset);
    get_timespec(&t->start);
    for (size_t i = 0; i < t->n_reps; ++i)
        t->code(t->data);
    get_timespec(&t->stop);
    return ((void*)0);
}
