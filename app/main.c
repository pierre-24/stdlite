#include <stdlib.h>
#include <string.h>

#include <stdlite/context.h>
#include <stdlite/helpers.h>

#include "app.h"

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
    struct timespec elapsed_time_app;
    stdl_timer_start(&elapsed_time_app);

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
                 "=========================================================\n"
                 "Program: %s\n"
                 "Library: %s (version %s)\n"
                 "Build date: %s\n"
                 "Build commit: %s\n"
                 "=========================================================\n",
                 APP_NAME, stdl_library_name(), stdl_library_version(), stdl_library_build_date(), stdl_library_build_commit()
                 );

    // Show what was understood
    title("User input");
    err = stdl_user_input_handler_log(input);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    title("Environment");
    err = stdl_app_log_env();
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // create context
    struct timespec elapsed_time_ctx;
    stdl_timer_start(&elapsed_time_ctx);

    title("Create context");
    err = stdl_user_input_handler_make_context(input, &ctx);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    stdl_log_msg(0, "Elapsed time in context: %.2f secs\n", stdl_timer_stop(&elapsed_time_ctx));

    // Compute responses
    /*struct timespec elapsed_time_response;
    stdl_timer_start(&elapsed_time_response);

    title("Compute responses");
    if(input->res_resreqs != NULL) {
        err = stdl_user_input_handler_prepare_responses(input, ctx, &rh);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        err = stdl_responses_handler_compute(rh, ctx);
        STDL_ERROR_CODE_HANDLE(err, goto _end);
    } else {
        stdl_log_msg(0, "No response requested.\n");
    }

    stdl_log_msg(0, "Elapsed time in responses: %.2f secs\n", stdl_timer_stop(&elapsed_time_response));*/


    // compute properties
    /*struct timespec elapsed_time_properties;
    stdl_timer_start(&elapsed_time_properties);

    title("Properties");
    if(input->res_resreqs != NULL) {
        err = stdl_user_input_handler_compute_properties(input, ctx, rh);
        STDL_ERROR_CODE_HANDLE(err, goto _end);
    } else {
        stdl_log_msg(0, "No properties.\n");
    }

    stdl_log_msg(0, "Elapsed time in properties: %.2f secs\n", stdl_timer_stop(&elapsed_time_properties));*/

    // the end
    _end:
    if(err < STDL_ERR_LAST) {
        title("End");
    }

    if(err == STDL_ERR_OK) {
        // report memory usage
        size_t user_input_sz = 0, ctx_sz = 0, res_sz = 0;

        stdl_log_msg(0, "Note: temporary memory allocated during run\n (including by linear algebra libraries)\n is not reported in the table below.\n RSS will be larger.\n");
        stdl_log_msg(0, "** Approximate memory usage ---\n");

        if(input != NULL) {
            size_t  resreq_sz;
            stdl_user_input_handler_approximate_size(input, &user_input_sz, &resreq_sz);

            double user_input_asz, resreq_asz;
            char* user_input_usz, *resreq_usz;
            stdl_convert_size(user_input_sz - resreq_sz, &user_input_asz, &user_input_usz);
            stdl_convert_size(resreq_sz, &resreq_asz, &resreq_usz);
            stdl_log_msg(0, "User input      REQ %8.1f%s\n", resreq_asz, resreq_usz);
            stdl_log_msg(0, "                OTH %8.1f%s\n", user_input_asz, user_input_usz);
        }

        if (ctx != NULL) {
            size_t wf_sz, bs_sz;
            stdl_context_approximate_size(ctx, &ctx_sz, &bs_sz, &wf_sz);

            double ctx_asz, wf_asz, bs_asz;
            char *ctx_usz, *wf_usz, *bs_usz;

            stdl_convert_size(ctx_sz - wf_sz - bs_sz, &ctx_asz, &ctx_usz);
            stdl_convert_size(wf_sz, &wf_asz, &wf_usz);
            stdl_convert_size(bs_sz, &bs_asz, &bs_usz);

            stdl_log_msg(0, "--------------  --- -----------\n");
            stdl_log_msg(0, "Context         WF  %8.1f%s\n", wf_asz, wf_usz);
            stdl_log_msg(0, "                BS  %8.1f%s\n", bs_asz, bs_usz);
            stdl_log_msg(0, "                OTH %8.1f%s\n", ctx_asz, ctx_usz);
        }

        if(rh != NULL) {
            size_t ev_sz, lrv_sz, amp_sz;
            stdl_responses_handler_approximate_size(rh, ctx->nmo, ctx->ncsfs, &res_sz, &ev_sz, &lrv_sz, &amp_sz);

            double res_asz, ev_asz, lrv_asz, amp_asz;
            char *res_usz, *ev_usz, *lrv_usz, *amp_usz;

            stdl_convert_size(res_sz - ev_sz - lrv_sz - amp_sz, &res_asz, &res_usz);
            stdl_convert_size(ev_sz, &ev_asz, &ev_usz);
            stdl_convert_size(lrv_sz, &lrv_asz, &lrv_usz);
            stdl_convert_size(amp_sz, &amp_asz, &amp_usz);

            stdl_log_msg(0, "--------------  --- -----------\n");
            stdl_log_msg(0, "Responses       INT %8.1f%s\n", ev_asz, ev_usz);
            stdl_log_msg(0, "                LRV %8.1f%s\n", lrv_asz, lrv_usz);
            stdl_log_msg(0, "                AMP %8.1f%s\n", amp_asz, amp_usz);
            stdl_log_msg(0, "                OTH %8.1f%s\n", res_asz, res_usz);
        }

        double tot_asz;
        char* tot_usz;
        stdl_convert_size(user_input_sz + ctx_sz + res_sz, &tot_asz, &tot_usz);
        stdl_log_msg(0, "--------------  --- -----------\n");
        stdl_log_msg(0, "Total               %8.1f%s\n", tot_asz, tot_usz);

    }

    if(input != NULL)
        stdl_user_input_handler_delete(input);

    if(ctx != NULL)
        stdl_context_delete(ctx);

    if(rh != NULL)
        stdl_responses_handler_delete(rh);

    if(err < STDL_ERR_LAST)
        stdl_log_msg(0, "Elapsed time in %s: %.2f secs\n", APP_NAME, stdl_timer_stop(&elapsed_time_app));

    if(err == STDL_ERR_OK || err > STDL_ERR_LAST) {
        return EXIT_SUCCESS;
    } else {
        stdl_log_msg(0, "\n\n** Something went bad :(\n");
        return err;
    }
}
