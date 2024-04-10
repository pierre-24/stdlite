#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/response.h>
#include <assert.h>
#include <string.h>

#include "responses_handler.h"

int stdl_op_data_new(stdl_operator op, size_t nmo, size_t ncsfs, size_t nw, stdl_op_data **data_ptr) {
    assert(data_ptr != NULL);

    *data_ptr = malloc(sizeof(stdl_op_data));
    STDL_ERROR_HANDLE_AND_REPORT(*data_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create op_data %p", *data_ptr);

    (*data_ptr)->op = op;
    (*data_ptr)->nw = nw;

    (*data_ptr)->op_ints_MO = malloc(STDL_OPERATOR_DIM[op] * STDL_MATRIX_SP_SIZE(nmo) * sizeof(double));
    (*data_ptr)->w = malloc(nw * sizeof(float ));
    (*data_ptr)->egrad = malloc(STDL_OPERATOR_DIM[op] * ncsfs * sizeof(float ));
    (*data_ptr)->X = malloc(nw * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
    (*data_ptr)->Y = malloc(nw * ncsfs * STDL_OPERATOR_DIM[op] * sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT((*data_ptr)->op_ints_MO == NULL || (*data_ptr)->w == NULL || (*data_ptr)->egrad == NULL || (*data_ptr)->X == NULL || (*data_ptr)->Y == NULL, stdl_op_data_delete(*data_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_op_data_delete(stdl_op_data* data) {

    assert(data != NULL);

    STDL_DEBUG("delete op_data %p", data);

    STDL_FREE_ALL(data->w, data->egrad, data->X, data->Y, data);

    return STDL_ERR_OK;
}

int stdl_op_data_approximate_size(stdl_op_data *data, size_t nmo, size_t ncsfs, size_t *sz) {
    assert(data != NULL && sz != NULL);

    *sz = sizeof(stdl_op_data)
            + STDL_MATRIX_SP_SIZE(nmo) * STDL_OPERATOR_DIM[data->op] * sizeof(double) // op_ints_MO
            + (data->nw * (1 + 3 * ncsfs * STDL_OPERATOR_DIM[data->op])) * sizeof(float ); // egrad + X + Y

    return STDL_ERR_OK;
}

int stdl_op_data_dump_h5(stdl_op_data *data, stdl_context *ctx, hid_t group_id) {
    assert(data != NULL && group_id != H5I_INVALID_HID);

    stdl_log_msg(1, "+ ");

    herr_t status;

    char* gname = STDL_OPERATOR_NAME[data->op];

    stdl_log_msg(0, "Saving LRV >");
    stdl_log_msg(1, "\n  | Create group `responses/%s` ", gname);

    hid_t lrv_group_id = H5Gcreate(group_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(lrv_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Store ints, egrad, and LRV ");

    status = H5LTmake_dataset(lrv_group_id, "info", 1, (hsize_t[]) {4}, H5T_NATIVE_ULONG, (size_t[]) {data->op, STDL_OPERATOR_DIM[data->op], (size_t) STDL_OPERATOR_HERMITIAN[data->op], data->nw});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "w", 1, (hsize_t[]) {data->nw}, H5T_NATIVE_FLOAT, data->w);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "integrals", 1, (hsize_t[]) {STDL_MATRIX_SP_SIZE(ctx->nmo)}, H5T_NATIVE_DOUBLE, data->op_ints_MO);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "egrad", 2, (hsize_t[]) {ctx->ncsfs, STDL_OPERATOR_DIM[data->op]}, H5T_NATIVE_FLOAT, data->egrad);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "X", 3, (hsize_t[]) {data->nw, ctx->ncsfs, STDL_OPERATOR_DIM[data->op]}, H5T_NATIVE_FLOAT, data->X);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    status = H5LTmake_dataset(lrv_group_id, "Y", 3, (hsize_t[]) {data->nw, ctx->ncsfs, STDL_OPERATOR_DIM[data->op]}, H5T_NATIVE_FLOAT, data->Y);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", lrv_group_id);

    stdl_log_msg(0, "< done\n");

    H5Gclose(lrv_group_id);
    return STDL_ERR_OK;
}

int stdl_op_data_compute(stdl_op_data *data, stdl_context *ctx) {
    assert(data != NULL && ctx != NULL);

    // get perturbed gradient
    int err = stdl_response_perturbed_gradient(ctx, STDL_OPERATOR_DIM[data->op], STDL_OPERATOR_HERMITIAN[data->op], data->op_ints_MO, data->egrad);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // compute response vectors
    if(ctx->B == NULL)
        err = stdl_response_TDA_linear(ctx, data->nw, data->w, STDL_OPERATOR_DIM[data->op], 1, data->egrad, data->X, data->Y);
    else
        err = stdl_response_TD_linear(ctx, data->nw, data->w, STDL_OPERATOR_DIM[data->op], 1, data->egrad, data->X, data->Y);

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
    (*rh_ptr)->Xamp = NULL;
    (*rh_ptr)->Yamp = NULL;

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

    STDL_DEBUG("delete responses_handler %p", rh);

    STDL_FREE_ALL(rh->eexci, rh->Xamp, rh->Yamp, rh);

    return STDL_ERR_OK;
}

int _dump_amplitudes_h5(hid_t group_id, stdl_context* ctx, stdl_responses_handler* rh) {

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Saving amplitudes >");
    stdl_log_msg(1, "\n  | Create group `responses/amplitudes` ");

    hid_t amp_group_id = H5Gcreate(group_id, "amplitudes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(amp_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | Store amplitudes ");

    herr_t status = H5LTmake_dataset(amp_group_id, "info", 1, (hsize_t[]) {1}, H5T_NATIVE_ULONG, (size_t[]) {rh->nexci});
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);

    status = H5LTmake_dataset(amp_group_id, "X", 2, (hsize_t[]) {rh->nexci, ctx->ncsfs}, H5T_NATIVE_FLOAT, rh->Xamp);
    STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);

    if(rh->Yamp != NULL) {
        status = H5LTmake_dataset(amp_group_id, "Y", 2, (hsize_t[]) {rh->nexci, ctx->ncsfs}, H5T_NATIVE_FLOAT, rh->Yamp);
        STDL_ERROR_HANDLE_AND_REPORT(status < 0, return STDL_ERR_WRITE, "cannot create dataset in group %d", amp_group_id);
    }

    H5Gclose(amp_group_id);

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}

/*
int stdl_responses_handler_compute(stdl_responses_handler* rh, stdl_context* ctx) {
    assert(rh != NULL && rh->nops > 0 && ctx != NULL);

    int err;

    // open H5 file
    hid_t file_id = H5Fopen(rh->data_output, H5F_ACC_RDWR, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", rh->data_output);

    // create group
    hid_t res_group_id = H5Gcreate(file_id, "responses", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(res_group_id == H5I_INVALID_HID, return STDL_ERR_WRITE, "cannot create group");

    stdl_log_msg(0, "~~ Compute ops_integrals in MO basis\n");

    for (size_t iop = 0; iop < rh->nops; ++iop) {
        size_t dim = STDL_OPERATOR_DIM[rh->ops[iop]];

        size_t ints_AO_sz = dim * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao) * sizeof(double );
        double ints_AO_asz;
        char* ints_AO_usz;
        stdl_convert_size(ints_AO_sz, &ints_AO_asz, &ints_AO_usz);
        stdl_log_msg(0, "Extra memory required (for AO integrals): %.1f%s\n", ints_AO_asz, ints_AO_usz);

        double* op_ints_AO = malloc(ints_AO_sz);
        STDL_ERROR_HANDLE_AND_REPORT(op_ints_AO == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

        err = stdl_operator_int1e_dsp(ctx->bs, STDL_OP_DIPL, (rh->ops[iop] == STDL_OP_DIPL? -1. :1.), op_ints_AO);
        STDL_ERROR_CODE_HANDLE(err, free(op_ints_AO); goto _end);

        for (size_t cpt = 0; cpt < dim; ++cpt) {
            char buff[128];
            sprintf(buff, "Component %ld", cpt);
            stdl_matrix_dsp_print(ctx->original_wf->nao, op_ints_AO + cpt * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao), buff);
        }

        // AO to MO
        rh->ops_integrals[iop] = malloc(dim * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));

        stdl_log_msg(1, "+ ");
        stdl_log_msg(0, "Converting integrals from AO to MO basis (dim=%ld) >", dim);

        for (size_t cpt = 0; cpt < dim; ++cpt) {
            stdl_log_msg(1, "\n  | component %ld ", cpt);
            err = stdl_wavefunction_dsp_ao_to_dsp_mo(
                    ctx->original_wf->nao,
                    ctx->nmo,
                    STDL_OPERATOR_HERMITIAN[rh->ops[iop]],
                    ctx->C_ptr,
                    op_ints_AO + cpt * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao),
                    rh->ops_integrals[iop] + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo));
            STDL_ERROR_CODE_HANDLE(err, free(op_ints_AO); goto _end);
            stdl_log_msg(0, "-");
        }


        stdl_log_msg(0, "< done\n");

        for (size_t cpt = 0; cpt < dim; ++cpt) {
            char buff[128];
            sprintf(buff, "Component %ld", cpt);
            stdl_matrix_dsp_print(ctx->nmo, rh->ops_integrals[iop] + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo), buff);
        }

        STDL_FREE_IF_USED(op_ints_AO);
    }

    if(rh->nlrvreqs > 0) {
        stdl_log_msg(0, "~~ Compute all linear response vectors\n");
        for (size_t ilrv = 0; ilrv < rh->nlrvreqs; ++ilrv) {

            stdl_lrv_request* lrvreq = rh->lrvreqs[ilrv];

            // find the correct operator
            for (size_t iop = 0; iop < rh->nops; ++iop) {
                if (rh->ops[iop] == lrvreq->op) {
                    lrvreq->op_integrals = rh->ops_integrals[iop];
                    break;
                }
            }

            // compute perturbed gradient
            err = stdl_lrv_request_compute(lrvreq, ctx);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            // dump
            err = stdl_lrv_request_dump_h5(lrvreq, ctx, res_group_id);
            STDL_ERROR_CODE_HANDLE(err, goto _end);
        }
    }

    if(rh->nexci > 0) {
        stdl_log_msg(0, "~~ Compute amplitude vectors\n");

        if(ctx->B == NULL)
            stdl_response_TDA_casida(ctx, rh->nexci, rh->eexci, rh->Xamp);
        else
            stdl_response_TD_casida(ctx, rh->nexci, rh->eexci, rh->Xamp, rh->Yamp);

        _dump_amplitudes_h5(res_group_id, ctx, rh);
    }

    _end:
    H5Gclose(res_group_id);
    H5Fclose(file_id);

    return err;
}*/

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

    *amp_sz = rh->nexci * (1 + (rh->Yamp == NULL? 1 : 2) * ncsfs) * sizeof(float);

    *sz = sizeof(stdl_responses_handler)
            + *ops_sz
            + *amp_sz;

    return STDL_ERR_OK;
}
