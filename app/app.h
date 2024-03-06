#ifndef STDLITE_APP_H
#define STDLITE_APP_H

#include <stdlite/context.h>

#include "user_input.h"

#ifndef APP_NAME
#define APP_NAME "stdlite_app"
#define APP_VERSION "unknown"
#endif

/**
 * Create user input from program input
 *
 * @param argc number of arguments
 * @param argv arguments
 * @param[out] inp
 * @return error code
 * @ingroup app
 */
int stdl_app_user_input(int argc, char* argv[], stdl_user_input** inp);

/**
 * Set the debug and log level
 * @return error code
 * @ingroup app
 */
int stdl_app_set_debug_log_level();

/**
 * Create context
 * @param inp a valid user input
 * @param[out] ctx_ptr context to be created
 * @return error code
 * @ingroup app
 */
int stdl_app_context(stdl_user_input* inp, stdl_context **ctx_ptr);

#endif //STDLITE_APP_H
