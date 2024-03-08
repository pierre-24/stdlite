#include <stdlib.h>
#include <string.h>

#include "app.h"

int stdl_app_set_debug_log_level() {
    char* end = NULL;
    char* loglevel = getenv("LOGLEVEL");
    if(loglevel != NULL) {
        long loglevel_d = strtol(loglevel, &end, 10);
        if(loglevel != end)
            stdl_set_log_level((int) loglevel_d);
    }

    char* debuglevel = getenv("DEBUGLEVEL");
    if(debuglevel != NULL) {
        long debuglevel_d = strtol(debuglevel, &end, 10);
        if(debuglevel != end)
            stdl_set_debug_level((int) debuglevel_d);
    }

    return STDL_ERR_OK;
}
