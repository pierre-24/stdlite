
#ifndef STDLITE_TIMER_H
#define STDLITE_TIMER_H

#include <time.h>

/** Create a timer
 */
void timer_start(struct timespec* t) {
    clock_gettime(CLOCK_MONOTONIC, t);
}

/** Get the elapsed time (in second)
 */
double timer_stop(struct timespec* start) {
    struct timespec stop;
    clock_gettime(CLOCK_MONOTONIC, &stop);
    return (double) (stop.tv_sec - start->tv_sec) + (double) (stop.tv_nsec - start->tv_nsec) * 1e-9;
}

#endif //STDLITE_TIMER_H
