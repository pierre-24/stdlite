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
    struct timespec elapsed_time_prog;
    stdl_timer_start(&elapsed_time_prog);

    int err;
    stdl_user_input_handler* input = NULL;
    stdl_responses_handler * rh = NULL;
    stdl_context* ctx = NULL;

    // Environment variables
    err = stdl_app_set_debug_log_level();
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // read input
    err = stdl_user_input_handler_new_from_args(argc, argv, &input);
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
    err = stdl_user_input_handler_log(input);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // create context
    struct timespec elapsed_time_ctx;
    stdl_timer_start(&elapsed_time_ctx);

    title("Create context");
    err = stdl_user_input_handler_make_context(input, &ctx);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    stdl_log_msg(0, "Elapsed time in context: %.2f secs\n", stdl_timer_stop(&elapsed_time_ctx));

    // Compute responses
    struct timespec elapsed_time_response;
    stdl_timer_start(&elapsed_time_response);

    title("Compute responses");
    if(input->res_resreqs != NULL) {
        err = stdl_user_input_handler_prepare_responses(input, ctx, &rh);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        err = stdl_responses_handler_compute(rh, ctx);
        STDL_ERROR_CODE_HANDLE(err, goto _end);
    } else {
        stdl_log_msg(0, "No response requested, exiting\n");
    }

    stdl_log_msg(0, "Elapsed time in responses: %.2f secs\n", stdl_timer_stop(&elapsed_time_response));

    // the end
    _end:
    if(err < STDL_ERR_LAST) {
        title("End");
    }

    if(input != NULL)
        stdl_user_input_handler_delete(input);

    if(ctx != NULL)
        stdl_context_delete(ctx);

    if(rh != NULL)
        stdl_responses_handler_delete(rh);

    if(err < STDL_ERR_LAST)
        stdl_log_msg(0, "Elapsed time in %s: %.2f secs\n", APP_NAME, stdl_timer_stop(&elapsed_time_prog));

    if(err == STDL_ERR_OK || err > STDL_ERR_LAST) {
        return EXIT_SUCCESS;
    } else {
        stdl_log_msg(0, "\n\n** Something went bad :(\n");
        return err;
    }
}
