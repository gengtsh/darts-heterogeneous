#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#define _GNU_SOURCE	       /* See feature_test_macros(7) */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <pthread.h>

static inline void fatal(int err, const char* msg) {
    errno = err;
    perror (  msg  );
    exit   ( errno );
}

static inline void spthread_create(pthread_t* thread, pthread_attr_t* attr, 
                                   void* (*start_routine)(void*), void* arg) {
        int err = pthread_create(thread,attr,start_routine,arg);
        if (0 != err) fatal(err, "pthread_create");
}

static inline void spthread_setaffinity_np(pthread_t thread, size_t cpusetsize, 
                                          const cpu_set_t *cpuset) {
    int err = pthread_setaffinity_np(thread, cpusetsize, cpuset);
    if (0 != err) fatal(err, "pthread_setaffinity_np");
}
 
static inline void* smalloc( const size_t sz ) {
    void* p = malloc ( sz );
    if ( !p ) fatal(errno,"malloc");
    return p;
}

static inline void sfree ( void* p ) {
    *(char*)p = 0;
    free(p);
}

#endif // UTIL_H_GUARD
