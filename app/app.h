#ifndef STDLITE_APP_H
#define STDLITE_APP_H

#include <stdlite/context.h>

#include "user_input.h"

#ifndef APP_NAME
#define APP_NAME "stdlite_app"
#define APP_VERSION "unknown"
#endif

/**
 * Set the debug and log level, using two environment variables: `LOGLEVEL` and `DEBUGLEVEL`
 * @return error code
 * @ingroup app
 */
int stdl_app_set_debug_log_level();

#endif //STDLITE_APP_H
