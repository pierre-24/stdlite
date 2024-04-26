#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/response.h>
#include <assert.h>
#include <string.h>

#include "responses_handler.h"
#include "log_property.h"

int stdl_op_data_new(stdl_operator op, size_t nmo, size_t ncsfs, size_t nlrvs, float *w, size_t *iw, stdl_op_data **data_ptr) {
    assert(data_ptr != NULL);

    *data_ptr = malloc(sizeof(stdl_op_data));
    STDL_ERROR_HANDLE_AND_REPORT(*data_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create op_data %p", *data_ptr);

    (*data_ptr)->op = op;
    (*data_ptr)->nlrvs = nlrvs;

    (*data_ptr)->op_ints_MO = malloc(STDL_OPERATOR_DIM[op] * STDL_MATRIX_SP_SIZE(nmo) * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT((*data_ptr)->op_ints_MO == NULL, stdl_op_data_delete(*data_ptr); return STDL_ERR_MALLOC, "malloc");

    (*data_ptr)->w = NULL;
    (*data_ptr)->iw = NULL;
    (*data_ptr)->egrad = NULL;
    (*data_ptr)->XpY = NULL;
    (*data_ptr)->XmY = NULL;
    (*data_ptr)->lrvs = NULL;

    if(nlrvs > 0) {
        (*data_ptr)->w = malloc(nlrvs * sizeof(float ));
        (*data_ptr)->iw = malloc(nlrvs * sizeof(size_t ));
        (*data_ptr)->egrad = malloc(STDL_OPERATOR_DIM[op] * ncsfs * sizeof(float ));
        (*data_ptr)->XpY = malloc(nlrvs * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
        (*data_ptr)->XmY = malloc(nlrvs * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
        (*data_ptr)->lrvs = malloc(nlrvs * sizeof(stdl_lrv));

        STDL_ERROR_HANDLE_AND_REPORT(
                (*data_ptr)->op_ints_MO == NULL
                    || (*data_ptr)->w == NULL
                    || (*data_ptr)->w == NULL
                    || (*data_ptr)->egrad == NULL
                    || (*data_ptr)->XpY == NULL
                    || (*data_ptr)->XmY == NULL
                    || (*data_ptr)->lrvs == NULL,
                stdl_op_data_delete(*data_ptr); return STDL_ERR_MALLOC, "malloc");

        memcpy((*data_ptr)->w, w, nlrvs * sizeof(float ));
        memcpy((*data_ptr)->iw, iw, nlrvs * sizeof(size_t ));
    }

    return STDL_ERR_OK;
}

int stdl_op_data_delete(stdl_op_data* data) {

    assert(data != NULL);

    STDL_DEBUG("delete op_data %p", data);

    STDL_FREE_ALL(data->w, data->iw, data->op_ints_MO, data->egrad, data->XpY, data->XmY, data->lrvs, data);

    return STDL_ERR_OK;
}

int stdl_op_data_approximate_size(stdl_op_data *data, size_t nmo, size_t ncsfs, size_t *sz) {
    assert(data != NULL && sz != NULL);

    *sz = sizeof(stdl_op_data)
            + STDL_MATRIX_SP_SIZE(nmo) * STDL_OPERATOR_DIM[data->op] * sizeof(double) // op_ints_MO
            + (data->nlrvs * (1 + 3 * ncsfs * STDL_OPERATOR_DIM[data->op])) * sizeof(float ); // egrad + X + Y

    return STDL_ERR_OK;
}

int stdl_op_data_dump_h5(stdl_op_data *data, stdl_context *ctx, hid_t group_id) {
    assert(data != NULL && group_id != H5I_INVALID_HID);

    stdl_log_msg(1, "+ ");

    herr_t status;

    char* gname = STDL_OPERATOR_NAME[data->op];

    stdl_log_msg(0, "Saving data >");
    stdl_log_msg(1, "\n  | Create group `responses/%s` ", gname);

    hid_t op_group_id = H5Gcreate(group_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(op_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Store integrals ");

    status = H5LTmake_dataset(op_group_id, "info", 1, (hsize_t[]) {4}, H5T_NATIVE_ULONG, (size_t[]) {data->op, STDL_OPERATOR_DIM[data->op], (size_t) STDL_OPERATOR_ISSYM[data->op], data->nlrvs});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);

    status = H5LTmake_dataset(op_group_id, "integrals", 1, (hsize_t[]) {STDL_MATRIX_SP_SIZE(ctx->nmo)}, H5T_NATIVE_DOUBLE, data->op_ints_MO);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);

    if(data->nlrvs > 0) {
        stdl_log_msg(0, "-");
        stdl_log_msg(1, "\n  | Store egrad, and LRV ");

        status = H5LTmake_dataset(op_group_id, "w", 1, (hsize_t[]) {data->nlrvs}, H5T_NATIVE_FLOAT, data->w);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);

        status = H5LTmake_dataset(op_group_id, "egrad", 2, (hsize_t[]) {ctx->ncsfs, STDL_OPERATOR_DIM[data->op]}, H5T_NATIVE_FLOAT, data->egrad);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);

        status = H5LTmake_dataset(op_group_id, "X+Y", 3, (hsize_t[]) {data->nlrvs, ctx->ncsfs, STDL_OPERATOR_DIM[data->op]}, H5T_NATIVE_FLOAT, data->XpY);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);

        status = H5LTmake_dataset(op_group_id, "X-Y", 3, (hsize_t[]) {data->nlrvs, ctx->ncsfs, STDL_OPERATOR_DIM[data->op]},H5T_NATIVE_FLOAT, data->XmY);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", op_group_id);
    }

    stdl_log_msg(0, "< done\n");

    H5Gclose(op_group_id);
    return STDL_ERR_OK;
}

int stdl_op_data_compute_lrvs(stdl_op_data *data, stdl_context *ctx) {
    assert(data != NULL && ctx != NULL);

    // get perturbed gradient
    int err = stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[data->op], STDL_OPERATOR_ISSYM[data->op], data->op_ints_MO, data->egrad);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // compute response vectors
    if(ctx->AmB == NULL)
        err = stdl_response_TDA_linear(ctx, data->nlrvs, data->w, STDL_OPERATOR_DIM[data->op], STDL_OPERATOR_HERMITIAN[data->op], data->egrad, data->XpY, data->XmY);
    else
        err = stdl_response_TD_linear(ctx, data->nlrvs, data->w, STDL_OPERATOR_DIM[data->op], STDL_OPERATOR_HERMITIAN[data->op], data->egrad, data->XpY, data->XmY);

    // distribute over `data->lrvs`
    for (size_t ilrv = 0; ilrv < data->nlrvs; ++ilrv) {
        data->lrvs[ilrv].op = data->op;
        data->lrvs[ilrv].op_ints_MO = data->op_ints_MO;
        data->lrvs[ilrv].w = data->w[ilrv];

        data->lrvs[ilrv].XpYw = data->XpY + ilrv * STDL_OPERATOR_DIM[data->op] * ctx->ncsfs;
        data->lrvs[ilrv].XmYw = data->XmY + ilrv * STDL_OPERATOR_DIM[data->op] * ctx->ncsfs;
    }

    return err;
}

int stdl_responses_handler_new(stdl_context *ctx, size_t nexci, stdl_responses_handler **rh_ptr) {
    assert(ctx != NULL && rh_ptr != NULL);

    *rh_ptr = malloc(sizeof(stdl_responses_handler));
    STDL_ERROR_HANDLE_AND_REPORT(*rh_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create responses_handler %p", *rh_ptr);

    (*rh_ptr)->nexci = nexci;

    for (int iop = 0; iop < STDL_OP_COUNT; ++iop) {
        (*rh_ptr)->lrvs_data[iop] = NULL;
    }

    (*rh_ptr)->eexci = NULL;
    (*rh_ptr)->XpYamp = NULL;
    (*rh_ptr)->XmYamp = NULL;

    if(nexci > 0) {
        (*rh_ptr)->eexci = malloc(nexci * sizeof(float ));
        (*rh_ptr)->XpYamp = malloc(ctx->ncsfs * nexci * sizeof(float ));

        STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->eexci == NULL || (*rh_ptr)->XpYamp == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");

        if(ctx->AmB != NULL) {
            (*rh_ptr)->XmYamp = malloc(ctx->ncsfs * nexci * sizeof(float));
            STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->XmYamp == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");
        }
    }

    return STDL_ERR_OK;
}

int stdl_responses_handler_new_from_input(stdl_user_input_handler* inp, stdl_context* ctx, stdl_responses_handler **rh_ptr) {
    assert(inp != NULL && inp->res_resreqs != NULL && ctx != NULL && rh_ptr != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Preparing responses >");
    stdl_log_msg(1, "\n  | Count requests ");

    int ops_needed[STDL_OP_COUNT] = {0};
    size_t res_nexci = 0, res_nops = 0, res_nlrvs = 0;

    int* ops_w = calloc(STDL_OP_COUNT * inp->res_nw, sizeof(int));

    STDL_ERROR_HANDLE_AND_REPORT(ops_w == NULL, return STDL_ERR_MALLOC, "malloc");

    stdl_response_request* req = inp->res_resreqs;
    while (req != NULL) {
        // check out if it contains a new operator
        for (size_t iop = 0; iop < req->nops; ++iop) {
            stdl_operator op = req->ops[iop];
            if(!ops_needed[op])
                res_nops += 1;

            ops_needed[op] = 1;
        }

        // check out if it contains new frequencies
        for (size_t ilrv = 0; ilrv < req->nlrvs; ++ilrv) {
            ops_w[req->ops[ilrv] * inp->res_nw + req->iw[ilrv]] = 1;
        }

        // check out if it requires amplitudes
        if(req->res_order > 0) {
            if(req->nroots < 0)
                res_nexci = ctx->ncsfs;
            else if((size_t) req->nroots > res_nexci) {
                if((size_t) req->nroots > ctx->ncsfs) {
                    STDL_WARN("%ld excited states requested, which is more than the number of CSFs", req->nroots);
                    req->nroots = (int) ctx->ncsfs;
                }
                res_nexci = (size_t) req->nroots;
            }
        }

        req = req->next;
    }

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Setup handler and data ");

    int err = stdl_responses_handler_new(ctx, res_nexci, rh_ptr);
    STDL_ERROR_CODE_HANDLE(err, return err);

    for (size_t iop = 0; iop < STDL_OP_COUNT; ++iop) {
        size_t nlrvs = 0;

        if(ops_needed[iop]) {
            float* op_w = malloc(inp->res_nw * sizeof(float ));
            size_t* op_iw = malloc(inp->res_nw * sizeof(size_t ));
            STDL_ERROR_HANDLE_AND_REPORT(op_w == NULL || op_iw == NULL, STDL_FREE_IF_USED(op_w); STDL_FREE_IF_USED(op_iw); return STDL_ERR_MALLOC, "malloc");

            STDL_DEBUG("%s selected", STDL_OPERATOR_NAME[iop]);
            for (size_t iw = 0; iw < inp->res_nw; ++iw) {
                if(ops_w[iop *  inp->res_nw + iw]) {
                    op_iw[nlrvs] = iw;
                    op_w[nlrvs] = inp->res_w[iw];
                    nlrvs += 1;
                }
            }

            err = stdl_op_data_new(iop, ctx->nmo, ctx->ncsfs, nlrvs, op_w, op_iw, &((*rh_ptr)->lrvs_data[iop]));
            STDL_ERROR_CODE_HANDLE(err,  STDL_FREE_IF_USED(op_iw); STDL_FREE_IF_USED(op_w); STDL_FREE_IF_USED(ops_w); return err);

            STDL_FREE_ALL(op_w, op_iw);

            res_nlrvs += nlrvs;
        }
    }

    stdl_log_msg(0, "< done\n");
    stdl_log_msg(0, "Will compute %ld matrix(ces) of MO integrals, %ld response vector(s), and %ld amplitude vector(s)\n", res_nops, res_nlrvs, res_nexci);

    STDL_FREE_ALL(ops_w);

    return STDL_ERR_OK;
}


int stdl_responses_handler_delete(stdl_responses_handler* rh) {
    assert(rh != NULL);

    STDL_DEBUG("delete responses_handler %p", rh);

    for (int iop = 0; iop < STDL_OP_COUNT; ++iop) {
        if(rh->lrvs_data[iop] != NULL)
            stdl_op_data_delete(rh->lrvs_data[iop]);
    }

    STDL_FREE_ALL(rh->eexci, rh->XpYamp, rh->XmYamp, rh);

    return STDL_ERR_OK;
}

int _dump_amplitudes_h5(stdl_responses_handler *rh, stdl_context *ctx, hid_t group_id) {

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Saving amplitudes >");
    stdl_log_msg(1, "\n  | Create group `responses/amplitudes` ");

    hid_t amp_group_id = H5Gcreate(group_id, "amplitudes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(amp_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Store amplitudes ");

    herr_t status = H5LTmake_dataset(amp_group_id, "info", 1, (hsize_t[]) {1}, H5T_NATIVE_ULONG, (size_t[]) {rh->nexci});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);

    status = H5LTmake_dataset(amp_group_id, "eexci", 1, (hsize_t[]) {rh->nexci}, H5T_NATIVE_FLOAT, rh->eexci);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);

    status = H5LTmake_dataset(amp_group_id, rh->XmYamp != NULL ? "X+Y" : "X", 2, (hsize_t[]) {rh->nexci, ctx->ncsfs}, H5T_NATIVE_FLOAT, rh->XpYamp);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);

    if(rh->XmYamp != NULL) {
        status = H5LTmake_dataset(amp_group_id, "X-Y", 2, (hsize_t[]) {rh->nexci, ctx->ncsfs}, H5T_NATIVE_FLOAT, rh->XmYamp);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);
    }

    H5Gclose(amp_group_id);

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}


int stdl_responses_handler_compute(stdl_responses_handler *rh, stdl_user_input_handler *inp, stdl_context *ctx) {
    assert(rh != NULL && inp != NULL && ctx != NULL);

    int err;

    // open H5 file
    hid_t file_id = H5Fopen(inp->data_output, H5F_ACC_RDWR, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", inp->data_output);

    // create group
    hid_t res_group_id = H5Gcreate(file_id, "responses", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(res_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    // amplitudes
    if(rh->nexci > 0) {
        stdl_log_msg(0, "~~ Compute amplitude vectors\n");

        if(ctx->AmB == NULL)
            stdl_response_TDA_casida(ctx, rh->nexci, rh->eexci, rh->XpYamp);
        else
            stdl_response_TD_casida(ctx, rh->nexci, rh->eexci, rh->XpYamp, rh->XmYamp);

        stdl_log_property_amplitude_contributions(rh, ctx, 5e-3f);
        _dump_amplitudes_h5(rh, ctx, res_group_id);
    }

    // integrals & LRVs (if any)
    for (size_t iop = 0; iop < STDL_OP_COUNT; ++iop) {

        stdl_op_data* op_data = rh->lrvs_data[iop];
        if(op_data == NULL)
            continue;

        stdl_log_msg(0, "~~ Compute integrals%s for `%s`\n", op_data->nlrvs > 0? " and linear response vectors": "", STDL_OPERATOR_NAME[op_data->op]);

        size_t ints_AO_sz = STDL_OPERATOR_DIM[op_data->op] * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao) * sizeof(double );
        double ints_AO_asz;
        char* ints_AO_usz;
        stdl_convert_size(ints_AO_sz, &ints_AO_asz, &ints_AO_usz);
        stdl_log_msg(0, "Extra memory required (for AO integrals): %.1f%s\n", ints_AO_asz, ints_AO_usz);

        double* op_ints_AO = malloc(ints_AO_sz);
        STDL_ERROR_HANDLE_AND_REPORT(op_ints_AO == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

        err = stdl_operator_int1e_dsp(ctx->bs, op_data->op, (op_data->op == STDL_OP_DIPL? -1. :1.), op_ints_AO);
        STDL_ERROR_CODE_HANDLE(err, free(op_ints_AO); goto _end);

        for (size_t cpt = 0; cpt < STDL_OPERATOR_DIM[op_data->op]; ++cpt) {
            char buff[128];
            sprintf(buff, "Component %ld", cpt);
            stdl_matrix_dsp_print(3, ctx->original_wf->nao,
                                  op_ints_AO + cpt * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao),
                                  buff);
        }

        stdl_log_msg(1, "+ ");
        stdl_log_msg(0, "Converting integrals from AO to MO basis (dim=%ld) >", STDL_OPERATOR_DIM[op_data->op]);

        for (size_t cpt = 0; cpt < STDL_OPERATOR_DIM[op_data->op]; ++cpt) {
            stdl_log_msg(1, "\n  | component %ld ", cpt);
            err = stdl_wavefunction_dsp_ao_to_dsp_mo(
                    ctx->original_wf->nao,
                    ctx->nmo,
                    STDL_OPERATOR_ISSYM[op_data->op],
                    ctx->C_ptr,
                    op_ints_AO + cpt * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao),
                    op_data->op_ints_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo));
            STDL_ERROR_CODE_HANDLE(err, free(op_ints_AO); goto _end);
            stdl_log_msg(0, "-");
        }

        stdl_log_msg(0, "< done\n");

        for (size_t cpt = 0; cpt < STDL_OPERATOR_DIM[op_data->op]; ++cpt) {
            char buff[128];
            sprintf(buff, "Component %ld", cpt);
            stdl_matrix_dsp_print(2, ctx->nmo, op_data->op_ints_MO + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo), buff);
        }

        STDL_FREE_IF_USED(op_ints_AO);

        if(op_data->nlrvs > 0) {
            err = stdl_op_data_compute_lrvs(op_data, ctx);
            STDL_ERROR_CODE_HANDLE(err, goto _end);
        }

        // dump
        err = stdl_op_data_dump_h5(op_data, ctx, res_group_id);
        STDL_ERROR_CODE_HANDLE(err, goto _end);
    }

    _end:
    H5Gclose(res_group_id);
    H5Fclose(file_id);

    return err;
}

// find the lrv data corresponding to a given request
int _find_lrvs(stdl_responses_handler* rh, stdl_user_input_handler* inp, stdl_response_request* req, stdl_lrv** lrvs) {
    assert(rh != NULL && inp != NULL && req != NULL && lrvs != NULL);

    for (size_t iop = 0; iop < req->nops; ++iop) {
        stdl_operator op = req->ops[iop];
        stdl_op_data* op_data = rh->lrvs_data[op];
        STDL_ERROR_HANDLE_AND_REPORT(op_data == NULL, return STDL_ERR_INPUT, "missing data for %s?!?", STDL_OPERATOR_NAME[op]);

        lrvs[iop] = NULL;

        for (size_t ilrv = 0; ilrv < op_data->nlrvs; ++ilrv) {
            if(op_data->iw[ilrv] == req->iw[iop]) {
                lrvs[iop] = &(op_data->lrvs[ilrv]);
                break;
            }
        }

        STDL_ERROR_HANDLE_AND_REPORT(lrvs[iop] == NULL, return STDL_ERR_INPUT, "cannot find lrv for %ld", iop);
    }

    return STDL_ERR_OK;
}

int stdl_response_handler_compute_properties(stdl_responses_handler* rh, stdl_user_input_handler* inp, stdl_context* ctx) {
    assert(rh != NULL && inp != NULL && ctx != NULL);

    int err = STDL_ERR_OK;
    hid_t prop_group_id = H5I_INVALID_HID;

    // open H5 file
    hid_t file_id = H5Fopen(inp->data_output, H5F_ACC_RDWR, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, err = STDL_ERR_OPEN; goto _end, "cannot open %s", inp->data_output);

    // create group
    prop_group_id = H5Gcreate(file_id, "properties", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(prop_group_id == H5I_INVALID_HID, err = STDL_ERR_WRITE; goto _end, "cannot create group");

    // add info
    herr_t status = H5LTmake_dataset(prop_group_id, "w", 1, (hsize_t[]) {inp->res_nw}, H5T_NATIVE_FLOAT, inp->res_w);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_WRITE; goto _end, "cannot create dataset in group %d", prop_group_id);

    size_t ireq = 0;
    stdl_response_request* req = inp->res_resreqs;
    char buffattr[256];
    while (req != NULL) {
        if(req->resp_order == 1 && req->res_order == 0) { // linear
            stdl_lrv* lrvs[] = {NULL, NULL};

            err = _find_lrvs(rh, inp, req, lrvs);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            req->property_tensor = malloc(STDL_OPERATOR_DIM[lrvs[0]->op] * STDL_OPERATOR_DIM[lrvs[1]->op] * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(req->property_tensor == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

            err = stdl_property_tensor_linear(ctx, lrvs, req->property_tensor);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            if(req->ops[0] == STDL_OP_DIPL && req->ops[1] == STDL_OP_DIPL)
                stdl_log_property_polarizability(req, req->property_tensor, inp->res_w[req->iw[1]]);
            else
                stdl_log_property_linear_tensor(req, req->property_tensor, inp->res_w[req->iw[1]]);

            sprintf(buffattr, "<<%s;%s>> @ wB=%f",
                    STDL_OPERATOR_NAME[req->ops[0]],
                    STDL_OPERATOR_NAME[req->ops[1]],
                    inp->res_w[req->iw[1]]);

        } else if(req->resp_order == 2 && req->res_order == 0) { // quadratic
            stdl_lrv* lrvs[] = {NULL, NULL, NULL};

            err = _find_lrvs(rh, inp, req, lrvs);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            sprintf(buffattr, "<<%s;%s,%s>> @ wB=%f wC=%f",
                    STDL_OPERATOR_NAME[req->ops[0]],
                    STDL_OPERATOR_NAME[req->ops[1]],
                    STDL_OPERATOR_NAME[req->ops[2]],
                    inp->res_w[req->iw[1]],
                    inp->res_w[req->iw[2]]
                    );

            req->property_tensor = malloc(STDL_OPERATOR_DIM[lrvs[0]->op] * STDL_OPERATOR_DIM[lrvs[1]->op] * STDL_OPERATOR_DIM[lrvs[2]->op] * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(req->property_tensor == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

            err = stdl_property_tensor_quadratic(ctx, lrvs, req->property_tensor);
            STDL_ERROR_CODE_HANDLE(err, goto _end);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            if(req->ops[0] == STDL_OP_DIPL && req->ops[1] == STDL_OP_DIPL && req->ops[2] == STDL_OP_DIPL)
                stdl_log_property_first_hyperpolarizability(req, req->property_tensor, inp->res_w[req->iw[1]], inp->res_w[req->iw[2]]);
            else
                stdl_log_property_quadratic_tensor(req, req->property_tensor, inp->res_w[req->iw[1]], inp->res_w[req->iw[2]]);

        } else if(req->resp_order == 1 && req->res_order == 1) { // linear SR
            req->property_tensor = malloc(rh->nexci * (STDL_OPERATOR_DIM[req->ops[0]] +  STDL_OPERATOR_DIM[req->ops[1]]) * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(req->property_tensor == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

            err = stdl_property_tensor_g2e_moments(
                    ctx,
                    req->ops,
                    (double *[]) {rh->lrvs_data[req->ops[0]]->op_ints_MO, rh->lrvs_data[req->ops[1]]->op_ints_MO},
                    req->nroots < 0 ? rh->nexci : req->nroots,
                    rh->XpYamp, rh->XmYamp,
                    req->property_tensor);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            stdl_log_property_g2e_moments(rh, ctx, req->ops, req->nroots < 0 ? rh->nexci : req->nroots, req->property_tensor);

            sprintf(buffattr, "lim_w <<%s;%s>>",
                    STDL_OPERATOR_NAME[req->ops[0]],
                    STDL_OPERATOR_NAME[req->ops[1]]
                    );
        }

        // write in H5
        char buff[256];
        sprintf(buff, "prop%ld", ireq);
        hid_t propx_group_id = H5Gcreate(prop_group_id, buff, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        STDL_ERROR_HANDLE_AND_REPORT(prop_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

        err = stdl_response_request_dump_h5(req, rh->nexci, propx_group_id);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        H5Gclose(propx_group_id);

        status = H5LTset_attribute_string(prop_group_id, buff, "name", buffattr);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, err = STDL_ERR_WRITE; goto _end, "cannot create attribute to group %d", propx_group_id);

        req = req->next;
        ireq++;
    }

    // write info in H5
    status = H5LTmake_dataset(prop_group_id, "info", 1, (hsize_t[]) {1}, H5T_NATIVE_ULONG, (size_t[]) {ireq});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", prop_group_id);

    _end:
    H5Gclose(prop_group_id);
    H5Fclose(file_id);

    return err;
}

int stdl_responses_handler_approximate_size(stdl_responses_handler *rh, size_t nmo, size_t ncsfs, size_t *sz, size_t *ops_sz, size_t *amp_sz) {
    assert(rh != NULL && sz != NULL && ops_sz != NULL && amp_sz != NULL);

    *ops_sz = 0;
    for (size_t iop = 0; iop < STDL_OP_COUNT; ++iop) {
        if(rh->lrvs_data[iop] != NULL) {
            size_t lsz = 0;
            stdl_op_data_approximate_size(rh->lrvs_data[iop], nmo, ncsfs, &lsz);
            *ops_sz += lsz;
        }
    }

    *amp_sz = rh->nexci * (1 + (rh->XmYamp == NULL ? 1 : 2) * ncsfs) * sizeof(float);

    *sz = sizeof(stdl_responses_handler)
            + *ops_sz
            + *amp_sz;

    return STDL_ERR_OK;
}
