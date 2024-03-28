#ifndef STDLITE_APP_H
#define STDLITE_APP_H

#include <time.h>

#include <stdlite/context.h>

#include "user_input_handler.h"

#ifndef APP_NAME
#define APP_NAME "stdlite_app"
#endif

/**
 * Set the debug and log level, using two environment variables: `STDL_LOG_LEVEL` and `STDL_DEBUG_LEVEL`
 * @return error code
 * @ingroup app_main
 */
int stdl_app_set_debug_log_level();

/**
 * Log the environment variables
 * @return error code
 * @ingroup app_main
 */
int stdl_app_log_env();

typedef struct timespec tmspec;

/** Create a timer
 * @param t a valid `timespec` structure
 * @ingroup app_main
 */
void stdl_timer_start(tmspec * t);

/** Get the elapsed time (in second)
 * @param start a valid `timespec` structure initialized with `stdl_timer_start`
 * @return the number of second since `start`
 * @ingroup app_main
 */
double stdl_timer_stop(tmspec * start);

#endif //STDLITE_APP_H
