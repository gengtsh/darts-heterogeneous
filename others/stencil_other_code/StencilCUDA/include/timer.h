#ifndef TIMER_H_GUARD
#define TIMER_H_GUARD

#include "util.h"
#include <time.h>
#include <stdint.h>

typedef struct timespec timespec_t;
typedef void* (*code_t)(void*);
typedef void*   data_t;

typedef struct loc_timer_s {
        timespec_t start, stop;
        long       cpu;
        code_t     code;
        data_t     data;
} loc_timer_t;

static inline void sclock_gettime(clockid_t clk_id, struct timespec *tp) {
    int res = clock_gettime ( clk_id, tp );
    if ( res ) fatal(errno, "clock_gettime");
}

static inline void get_timespec(timespec_t* t) { sclock_gettime(CLOCK_REALTIME, t); }

static inline uint64_t get_time ( void ) {
    timespec_t time;
    clock_gettime( CLOCK_REALTIME, &time );
    return (uint64_t)( time.tv_sec*1000000000UL + time.tv_nsec );
}

static inline void
timespec_diff ( timespec_t* end, timespec_t* start, timespec_t* result ) {
    if ( end->tv_nsec-start->tv_nsec < 0 ) {
        /* Transfer one second from seconds to nanoseconds  */
        result->tv_sec  =         -1 + end->tv_sec  - start->tv_sec;  
        result->tv_nsec = 1000000000 + end->tv_nsec - start->tv_nsec;  
    } else {
        result->tv_sec  =              end->tv_sec  - start->tv_sec;
        result->tv_nsec =              end->tv_nsec - start->tv_nsec;
    }
}

void* timer(void* arg);


inline uint64_t getTime(void)
{
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);

/*TODO: Determine which one of these is the best.*/
#if 0
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
#endif
#if 0
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time);
#endif
    return ((uint64_t)time.tv_sec) * 1000000000 + time.tv_nsec;
}


#endif // TIMER_H_GUARD
