#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <string.h>

#include "response_requests.h"

int stdl_operator_dim(stdl_operator op, size_t* dim) {
    assert(op < STDL_OP_COUNT && dim != NULL);

    switch (op) {
        case STDL_OP_DIPL:
            *dim = 3;
            break;
        default:
            *dim = 1;
    }

    return STDL_ERR_OK;
}

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroot, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL && resp_order > 0);

    *req_ptr = malloc(sizeof(stdl_response_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create response_request %p", *req_ptr);

    (*req_ptr)->resp_order = resp_order;
    (*req_ptr)->res_order = res_order;
    (*req_ptr)->nroot = nroot;
    (*req_ptr)->next = NULL;

    (*req_ptr)->w = NULL;
    (*req_ptr)->requests = NULL;
    (*req_ptr)->wpos = NULL;

    size_t nops = resp_order-res_order+1;
    (*req_ptr)->ops = malloc(nops * sizeof(stdl_operator));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->ops == NULL, stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");
    memcpy((*req_ptr)->ops, ops, nops * sizeof(stdl_operator));

    size_t nw = (resp_order == res_order)? 0: resp_order-res_order+1;
    if(nw > 0) {
        (*req_ptr)->w = malloc(nw * sizeof(float));
        (*req_ptr)->requests = malloc(nw * sizeof(stdl_lrv_request*));
        (*req_ptr)->wpos = malloc(nw * sizeof(size_t));
        STDL_ERROR_HANDLE_AND_REPORT(
                (*req_ptr)->ops == NULL || (*req_ptr)->w == NULL || (*req_ptr)->requests == NULL || (*req_ptr)->wpos == NULL,
                stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

        memcpy((*req_ptr)->w, w, nw * sizeof(float ));
    }

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

int stdl_lrv_request_new(stdl_operator op, size_t nw, size_t ncsfs, stdl_lrv_request **req_ptr) {
    assert(req_ptr != NULL);

    *req_ptr = malloc(sizeof(stdl_lrv_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create lrv_request %p", *req_ptr);

    (*req_ptr)->op = op;
    (*req_ptr)->nw = nw;

    size_t dim;
    int err = stdl_operator_dim(op, &dim);
    STDL_ERROR_HANDLE(err, stdl_lrv_request_delete(*req_ptr); return err);

    (*req_ptr)->w = malloc(nw * sizeof(float ));
    (*req_ptr)->X = malloc(nw * ncsfs * dim * sizeof(float ));
    (*req_ptr)->Y = malloc(nw * ncsfs * dim * sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->w == NULL || (*req_ptr)->X == NULL || (*req_ptr)->Y == NULL, stdl_lrv_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_lrv_request_delete(stdl_lrv_request* req) {
    assert(req != NULL);

    STDL_DEBUG("delete lrv_request %p", req);

    STDL_FREE_ALL(req->w, req->X, req->Y, req);

    return STDL_ERR_OK;
}
