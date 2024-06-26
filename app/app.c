#include <stdlib.h>
#include <string.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "app.h"

int stdl_app_set_debug_log_level() {
    char* end = NULL;
    char* loglevel = getenv("STDL_LOG_LEVEL");
    if(loglevel != NULL) {
        long loglevel_d = strtol(loglevel, &end, 10);
        if(loglevel != end)
            stdl_set_log_level((int) loglevel_d);
    }

    char* debuglevel = getenv("STDL_DEBUG_LEVEL");
    if(debuglevel != NULL) {
        long debuglevel_d = strtol(debuglevel, &end, 10);
        if(debuglevel != end)
            stdl_set_debug_level((int) debuglevel_d);
    }

    return STDL_ERR_OK;
}

int stdl_app_log_env() {
    stdl_log_msg(0, "STDL_LOG_LEVEL=%d\n", stdl_get_log_level());
    stdl_log_msg(0, "STDL_DEBUG_LEVEL=%d\n", stdl_get_debug_level());

#ifdef USE_OPENMP
    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
            stdl_log_msg(0, "OMP_NUM_THREADS=%d\n", omp_get_num_threads());
    }
#endif

    return STDL_ERR_OK;
}


void stdl_timer_start(struct timespec* t) {
    clock_gettime(CLOCK_MONOTONIC, t);
}

double stdl_timer_stop(struct timespec* start) {
    struct timespec stop;
    clock_gettime(CLOCK_MONOTONIC, &stop);
    return (double) (stop.tv_sec - start->tv_sec) + (double) (stop.tv_nsec - start->tv_nsec) * 1e-9;
}
