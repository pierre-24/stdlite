#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <assert.h>
#include <string.h>

#include "responses_handler.h"


int stdl_responses_handler_new(size_t nops, size_t nlrvreqs, size_t nexci, stdl_context *ctx, stdl_responses_handler **rh_ptr) {
    assert(nops > 0 && ctx != NULL && rh_ptr != NULL);

    *rh_ptr = malloc(sizeof(stdl_responses_handler));
    STDL_ERROR_HANDLE_AND_REPORT(*rh_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*rh_ptr)->nops = nops;
    (*rh_ptr)->nlrvreqs = nlrvreqs;
    (*rh_ptr)->nexci = nexci;

    (*rh_ptr)->eexci = NULL;
    (*rh_ptr)->Xamp = NULL;
    (*rh_ptr)->Yamp = NULL;

    (*rh_ptr)->ops = malloc(nops * sizeof(stdl_operator));
    (*rh_ptr)->lrvreqs = calloc(nops, sizeof(stdl_lrv_request*));

    STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->ops == NULL || (*rh_ptr)->lrvreqs == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");

    if(nexci > 0) {
        (*rh_ptr)->eexci = malloc(nexci * sizeof(float ));
        (*rh_ptr)->Xamp = malloc(ctx->ncsfs * nexci * sizeof(float ));

        STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->eexci == NULL || (*rh_ptr)->Xamp == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");

        if(ctx->B != NULL) {
            (*rh_ptr)->Yamp = malloc(ctx->ncsfs * nexci * sizeof(float));
            STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->Yamp == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");
        }
    }

    return STDL_ERR_OK;
}


int stdl_responses_handler_delete(stdl_responses_handler* rh) {
    assert(rh != NULL);

    for (size_t i = 0; i < rh->nlrvreqs; ++i) {
        if(rh->lrvreqs[i] != NULL)
            stdl_lrv_request_delete(rh->lrvreqs[i]);
    }

    STDL_FREE_ALL(rh->ops, rh->lrvreqs, rh->eexci, rh->Xamp, rh->Yamp, rh);

    return STDL_ERR_OK;
}
