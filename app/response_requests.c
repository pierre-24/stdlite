#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/matrix.h>
#include <string.h>

#include "response_requests.h"

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, size_t *iw, int nroots, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL && resp_order > 0);

    *req_ptr = malloc(sizeof(stdl_response_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create response_request %p", *req_ptr);

    (*req_ptr)->resp_order = resp_order;
    (*req_ptr)->res_order = res_order;
    (*req_ptr)->nroots = nroots;
    (*req_ptr)->nops = resp_order-res_order+1;
    (*req_ptr)->nlrvs = (resp_order == res_order) ? 0 : (*req_ptr)->nops;
    (*req_ptr)->iw = NULL;
    (*req_ptr)->next = NULL;

    (*req_ptr)->ops = malloc((*req_ptr)->nops * sizeof(stdl_operator));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->ops == NULL, stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

    memcpy((*req_ptr)->ops, ops, (*req_ptr)->nops * sizeof(stdl_operator));

    if((*req_ptr)->nlrvs > 0) {
        (*req_ptr)->iw = malloc((*req_ptr)->nlrvs * sizeof(size_t));
        STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->iw == NULL, stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

        memcpy((*req_ptr)->iw, iw, (*req_ptr)->nlrvs * sizeof(size_t ));
    }

    return STDL_ERR_OK;
}

int stdl_response_request_delete(stdl_response_request* req) {
    assert(req != NULL);

    STDL_DEBUG("delete response_request %p", req);

    if(req->next != NULL)
        stdl_response_request_delete(req->next);

    STDL_FREE_ALL(req->ops, req->iw, req);

    return STDL_ERR_OK;

}

int stdl_response_request_approximate_size(stdl_response_request *req, size_t *sz) {
    assert(req != NULL && sz != NULL);

    size_t next_sz = 0;

    if(req->next != NULL)
        stdl_response_request_approximate_size(req->next, &next_sz);

    *sz = sizeof(stdl_response_request)
            + req->nops * sizeof(size_t) // ops
            + req->nlrvs * (sizeof(size_t) + sizeof(stdl_lrv)) // iw + lrvs
            + next_sz;

    return STDL_ERR_OK;
}
