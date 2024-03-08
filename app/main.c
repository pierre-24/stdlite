#include <stdlib.h>
#include <string.h>

#include <stdlite/context.h>

#include "app.h"
#include "timer.h"

#define NDASHES 80

void title(char* title) {
    int sz = strlen(title);
    stdl_log_msg(0, "\n--- %s -", title);
    for (int i = 0; i < NDASHES - sz - 6; ++i) {
        stdl_log_msg(0, "-");
    }
    stdl_log_msg(0, "\n");
}

int main(int argc, char* argv[]) {

    // record time
    struct timespec elapsed_prog;
    stdl_timer_start(&elapsed_prog);

    int err;
    stdl_user_input* input = NULL;
    stdl_context* ctx = NULL;

    // Environment variables
    err = stdl_app_set_debug_log_level();
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // read input
    err = stdl_app_user_input(argc, argv, &input);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // scream aloud
    stdl_log_msg(0,
                 "==========================================\n"
                 "Running %s\n"
                 "using %s (version %s)\n"
                 "==========================================\n",
                 APP_NAME, stdl_library_name(), stdl_library_version()
                 );

    title("User input");
    err = stdl_user_input_log(input);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // create context
    title("Create context");
    err = stdl_app_context(input, &ctx);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // the end
    _end:
    if(input != NULL)
        stdl_user_input_delete(input);

    if(ctx != NULL)
        stdl_context_delete(ctx);

    if(err < STDL_ERR_LAST) {
        title("End");
        stdl_log_msg(0, "Elapsed time in %s: %.2f secs\n", APP_NAME, stdl_timer_stop(&elapsed_prog));
    }

    if(err == STDL_ERR_OK || err > STDL_ERR_LAST) {
        return EXIT_SUCCESS;
    } else {
        stdl_log_msg(0, "\n\n** Something went bad :(\n");
        return err;
    }
}
