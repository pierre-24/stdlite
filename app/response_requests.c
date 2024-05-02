#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/matrix.h>
#include <string.h>

#include "response_requests.h"

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, size_t *iw, int nroots, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL && resp_order > 0 && res_order <= resp_order);

    *req_ptr = malloc(sizeof(stdl_response_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create response_request %p", *req_ptr);

    (*req_ptr)->resp_order = resp_order;
    (*req_ptr)->res_order = res_order;
    (*req_ptr)->nroots = nroots;
    (*req_ptr)->nops = resp_order+1;
    (*req_ptr)->nlrvs = (resp_order == res_order) ? 0 : (*req_ptr)->nops;
    (*req_ptr)->iw = NULL;
    (*req_ptr)->property_tensor = NULL;
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

    STDL_FREE_ALL(req->ops, req->iw, req->property_tensor, req);

    return STDL_ERR_OK;
}

int stdl_response_request_approximate_size(stdl_response_request *req, size_t nexci, size_t *sz) {
    assert(req != NULL && sz != NULL);

    size_t next_sz = 0, tensor_sz = 0;

    if(req->next != NULL)
        stdl_response_request_approximate_size(req->next, nexci, &next_sz);

    if(req->property_tensor != NULL) {
        if(req->res_order == 0) {
            for (size_t iop = 0; iop < req->nops; ++iop) {
                tensor_sz *= STDL_OPERATOR_DIM[iop];
            }
        } else {
            for (size_t iop = 0; iop < req->nops; ++iop) {
                tensor_sz += STDL_OPERATOR_DIM[iop];
            }

            tensor_sz *= (req->nroots < 0 ? nexci : (size_t) req->nroots);
        }

        tensor_sz *= sizeof(float);
    }

    *sz = sizeof(stdl_response_request)
            + req->nops * sizeof(size_t) // ops
            + req->nlrvs * (sizeof(size_t) + sizeof(stdl_lrv)) // iw + lrvs
            + tensor_sz
            + next_sz;

    return STDL_ERR_OK;
}

int stdl_response_request_dump_h5(stdl_response_request *req, size_t maxnexci, hid_t group_id) {
    assert(req != NULL && group_id != H5I_INVALID_HID);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Saving property >");
    stdl_log_msg(1, "\n  | Saving data ");

    herr_t status = H5LTmake_dataset(group_id, "info", 1, (hsize_t[]) {5}, H5T_NATIVE_ULONG, (size_t[]) {req->resp_order, req->res_order,  req->nops, req->nlrvs, req->res_order > 0 ? maxnexci : 0});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", group_id);

    status = H5LTmake_dataset(group_id, "ops", 1, (hsize_t[]) {(hsize_t) req->nops}, H5T_NATIVE_INT, req->ops);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", group_id);

    if(req->nlrvs > 0) {
        status = H5LTmake_dataset(group_id, "w", 1, (hsize_t[]) {(hsize_t) req->nlrvs}, H5T_NATIVE_ULONG, req->iw);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", group_id);
    }

    if(req->property_tensor != NULL) {

        // get its shape
        hsize_t dims[3] = {0};
        int rank = 0;

        if (req->resp_order == 1 && req->res_order == 0) {
            dims[0] = (hsize_t) STDL_OPERATOR_DIM[req->ops[0]];
            dims[1] = (hsize_t) STDL_OPERATOR_DIM[req->ops[1]];
            rank = 2;
        } else if (req->resp_order == 2 && req->res_order == 0) {
            dims[0] = (hsize_t) STDL_OPERATOR_DIM[req->ops[0]];
            dims[1] = (hsize_t) STDL_OPERATOR_DIM[req->ops[1]];
            dims[2] = (hsize_t) STDL_OPERATOR_DIM[req->ops[2]];
            rank = 3;
        } else if (req->resp_order == 1 && req->res_order == 1) {
            dims[0] = (hsize_t) (STDL_OPERATOR_DIM[req->ops[0]] + STDL_OPERATOR_DIM[req->ops[1]]);
            dims[1] = (hsize_t) (req->nroots < 0? maxnexci : req->nroots);
            rank = 2;
        } else if (req->resp_order == 2 && req->res_order == 2) {
            dims[0] = (hsize_t) (STDL_OPERATOR_DIM[req->ops[0]] + STDL_OPERATOR_DIM[req->ops[1]] + STDL_OPERATOR_DIM[req->ops[2]]);
            dims[1] = (hsize_t) STDL_MATRIX_SP_SIZE(req->nroots < 0? maxnexci : req->nroots);
            rank = 2;
        }

        // write it
        status = H5LTmake_dataset(group_id, "property_tensor", rank, dims, H5T_NATIVE_FLOAT, req->property_tensor);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", group_id);
    }

    stdl_log_msg(0, "< done\n");
    return STDL_ERR_OK;
}
