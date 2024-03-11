#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <string.h>

#include "requests.h"

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroot, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL && resp_order > 0);

    size_t nw = resp_order - res_order;

    *req_ptr = malloc(sizeof(stdl_response_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create response_request %p", *req_ptr);

    (*req_ptr)->resp_order = resp_order;
    (*req_ptr)->res_order = res_order;
    (*req_ptr)->nroot = nroot;

    (*req_ptr)->wpos = 0;
    (*req_ptr)->next = NULL;

    (*req_ptr)->ops = malloc((nw + 1) * sizeof(stdl_operator));
    (*req_ptr)->w = malloc(nw * sizeof(float ));
    (*req_ptr)->requests = malloc(nw * sizeof (stdl_lrv_request*));
    (*req_ptr)->wpos = malloc(nw * sizeof(size_t ));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->ops == NULL || (*req_ptr)->w == NULL || (*req_ptr)->requests == NULL || (*req_ptr)->wpos == NULL, stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

    memcpy((*req_ptr)->ops, ops, (nw + 1) * sizeof(stdl_operator));
    memcpy((*req_ptr)->w, w, nw * sizeof(float ));

    return STDL_ERR_OK;
}

int stdl_response_request_delete(stdl_response_request* req) {
    assert(req != NULL);

    STDL_DEBUG("delete response_request %p", req);

    if(req->next != NULL)
        stdl_response_request_delete(req->next);

    STDL_FREE_ALL(req->ops, req->w, req->requests, req->wpos, req);

    return STDL_ERR_OK;

}

int stdl_lrv_request_new(stdl_operator op, size_t nw, stdl_lrv_request** req_ptr) {
    assert(req_ptr != NULL);

    *req_ptr = malloc(sizeof(stdl_lrv_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create lrv_request %p", *req_ptr);

    (*req_ptr)->op = op;
    (*req_ptr)->nw = nw;

    (*req_ptr)->X = NULL;
    (*req_ptr)->Y = NULL;

    (*req_ptr)->w = malloc(nw * sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->w == NULL, stdl_lrv_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_lrv_request_delete(stdl_lrv_request* req) {
    assert(req != NULL);

    STDL_DEBUG("delete lrv_request %p", req);

    STDL_FREE_ALL(req->w, req->X, req->Y, req);

    return STDL_ERR_OK;
}
