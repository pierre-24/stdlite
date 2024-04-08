#include <assert.h>
#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <string.h>

#include "response_requests.h"

int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroots, stdl_response_request** req_ptr) {
    assert(req_ptr != NULL && resp_order > 0);

    *req_ptr = malloc(sizeof(stdl_response_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create response_request %p", *req_ptr);

    (*req_ptr)->resp_order = resp_order;
    (*req_ptr)->res_order = res_order;
    (*req_ptr)->nroots = nroots;
    (*req_ptr)->next = NULL;

    (*req_ptr)->w = NULL;
    (*req_ptr)->lrvreqs = NULL;
    (*req_ptr)->wpos = NULL;

    size_t nops = resp_order-res_order+1;
    (*req_ptr)->ops = malloc(nops * sizeof(stdl_operator));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->ops == NULL, stdl_response_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");
    memcpy((*req_ptr)->ops, ops, nops * sizeof(stdl_operator));

    size_t nw = (resp_order == res_order)? 0: resp_order-res_order+1;
    if(nw > 0) {
        (*req_ptr)->w = malloc(nw * sizeof(float));
        (*req_ptr)->lrvreqs = malloc(nw * sizeof(stdl_lrv_request*));
        (*req_ptr)->wpos = malloc(nw * sizeof(size_t));
        STDL_ERROR_HANDLE_AND_REPORT(
                (*req_ptr)->ops == NULL || (*req_ptr)->w == NULL || (*req_ptr)->lrvreqs == NULL || (*req_ptr)->wpos == NULL,
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

    STDL_FREE_ALL(req->ops, req->w, req->lrvreqs, req->wpos, req);

    return STDL_ERR_OK;

}

int stdl_lrv_request_new(stdl_operator op, size_t nw, size_t ncsfs, stdl_lrv_request **req_ptr) {
    assert(req_ptr != NULL);

    *req_ptr = malloc(sizeof(stdl_lrv_request));
    STDL_ERROR_HANDLE_AND_REPORT(*req_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create lrv_request %p", *req_ptr);

    (*req_ptr)->op = op;
    (*req_ptr)->nw = nw;

    (*req_ptr)->w = malloc(nw * sizeof(float ));
    (*req_ptr)->egrad = malloc(STDL_OPERATOR_DIM[op] * ncsfs * sizeof(float ));
    (*req_ptr)->X = malloc(nw * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
    (*req_ptr)->Y = malloc(nw * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT((*req_ptr)->w == NULL || (*req_ptr)->egrad == NULL || (*req_ptr)->X == NULL || (*req_ptr)->Y == NULL, stdl_lrv_request_delete(*req_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_lrv_request_delete(stdl_lrv_request* req) {
    assert(req != NULL);

    STDL_DEBUG("delete lrv_request %p", req);

    STDL_FREE_ALL(req->w, req->egrad, req->X, req->Y, req);

    return STDL_ERR_OK;
}

int stdl_lrv_request_compute(stdl_lrv_request *lrvreq, stdl_context *ctx) {
    assert(lrvreq != NULL && ctx != NULL && lrvreq->op_integrals != NULL);

    // get perturbed gradient
    int err = stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[lrvreq->op], STDL_OPERATOR_HERMITIAN[lrvreq->op], lrvreq->op_integrals, lrvreq->egrad);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // compute response vectors
    if(ctx->B == NULL)
        err = stdl_response_TDA_linear(ctx, lrvreq->nw, lrvreq->w, STDL_OPERATOR_DIM[lrvreq->op], lrvreq->egrad, lrvreq->X, lrvreq->Y);
    else
        err = stdl_response_TD_linear(ctx, lrvreq->nw, lrvreq->w, STDL_OPERATOR_DIM[lrvreq->op], lrvreq->egrad, lrvreq->X, lrvreq->Y);


    return err;
}

int stdl_lrv_request_dump_h5(stdl_lrv_request *req, stdl_context *ctx, hid_t group_id) {
    assert(req != NULL && group_id != H5I_INVALID_HID);

    stdl_log_msg(1, "+ ");

    herr_t status;

    char* gname;
    switch (req->op) {
        case STDL_OP_DIPL:
            gname = "dipl";
            break;
        default:
            gname = "unk";
            break;
    }

    stdl_log_msg(0, "Saving LRV >");
    stdl_log_msg(1, "\n  | Create group `responses/%s` ", gname);

    hid_t lrv_group_id = H5Gcreate(group_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(lrv_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Store ints, egrad, and LRV ");

    status = H5LTmake_dataset(lrv_group_id, "info", 1, (hsize_t[]) {4}, H5T_NATIVE_ULONG, (size_t[]) {req->op, STDL_OPERATOR_DIM[req->op], (size_t) STDL_OPERATOR_HERMITIAN[req->op],req->nw});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "w", 1, (hsize_t[]) {req->nw}, H5T_NATIVE_FLOAT, req->w);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "integrals", 1, (hsize_t[]) {STDL_MATRIX_SP_SIZE(ctx->nmo)}, H5T_NATIVE_DOUBLE, req->op_integrals);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "egrad", 2, (hsize_t[]) {ctx->ncsfs, STDL_OPERATOR_DIM[req->op]}, H5T_NATIVE_FLOAT, req->egrad);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "X", 3, (hsize_t[]) {req->nw, ctx->ncsfs, STDL_OPERATOR_DIM[req->op]}, H5T_NATIVE_FLOAT, req->X);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "Y", 3, (hsize_t[]) {req->nw, ctx->ncsfs, STDL_OPERATOR_DIM[req->op]}, H5T_NATIVE_FLOAT, req->Y);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    stdl_log_msg(0, "< done\n");

    H5Gclose(lrv_group_id);
    return STDL_ERR_OK;
}

int stdl_lrv_request_approximate_size(stdl_lrv_request *req, size_t ncsfs, size_t *sz) {
    assert(req != NULL && sz != NULL);

    *sz = sizeof(stdl_lrv_request)
            + (req->nw * (1 + 3 * ncsfs * STDL_OPERATOR_DIM[req->op])) * sizeof(float );

    return STDL_ERR_OK;
}

int stdl_response_request_approximate_size(stdl_response_request *req, size_t *sz) {
    assert(req != NULL && sz != NULL);

    size_t nops = req->resp_order-req->res_order+1, nw = (req->resp_order == req->res_order)? 0: req->resp_order-req->res_order+1, next_sz = 0;

    if(req->next != NULL)
        stdl_response_request_approximate_size(req->next, &next_sz);

    *sz = sizeof(stdl_response_request)
            + nops * sizeof(size_t)
            + nw * (sizeof(float ) + sizeof(size_t) + sizeof(stdl_lrv_request*))
            + next_sz;

    return STDL_ERR_OK;
}
