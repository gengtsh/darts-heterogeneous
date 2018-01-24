#include "timer.h"
#include <pthread.h>
//#include <sched.h>

void*
timer(void* arg)
{
    loc_timer_t*  t = (loc_timer_t*) arg;
    cpu_set_t cpuset;

    CPU_ZERO ( &cpuset );
    CPU_SET  ( t->cpu, &cpuset );

    spthread_setaffinity_np(pthread_self(), 1UL, &cpuset);

    get_timespec(&t->start);
    for (size_t i = 0; i < 20; ++i)
        t->code(t->data);
    get_timespec(&t->stop);

    return NULL;
}


