#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>

#include "requests.h"

int stdl_lrv_request_new(stdl_operator op, size_t nw, stdl_lrv_request** req_prt) {
    assert(req_prt != NULL);

    return STDL_ERR_OK;
}

int stdl_lrv_request_delete(stdl_lrv_request* req) {
    assert(req != NULL);

    STDL_FREE_ALL(req);

    return STDL_ERR_OK;
}

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroot, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL);

    return STDL_ERR_OK;
}

int stdl_response_request_delete(stdl_response_request* req) {
    assert(req != NULL);

    STDL_FREE_ALL(req);

    return STDL_ERR_OK;

}