#ifndef STDLITE_APP_H
#define STDLITE_APP_H

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

#endif //STDLITE_APP_H
