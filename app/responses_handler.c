#include <stdlite/logging.h>
#include <stdlite/helpers.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/response.h>
#include <assert.h>
#include <string.h>

#include "responses_handler.h"


int stdl_responses_handler_new(size_t nops, size_t nlrvreqs, size_t nexci, stdl_context *ctx, stdl_responses_handler **rh_ptr) {
    assert(nops > 0 && ctx != NULL && rh_ptr != NULL);

    *rh_ptr = malloc(sizeof(stdl_responses_handler));
    STDL_ERROR_HANDLE_AND_REPORT(*rh_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create responses_hanlder %p", *rh_ptr);

    (*rh_ptr)->nops = nops;
    (*rh_ptr)->nlrvreqs = nlrvreqs;
    (*rh_ptr)->nexci = nexci;

    (*rh_ptr)->eexci = NULL;
    (*rh_ptr)->Xamp = NULL;
    (*rh_ptr)->Yamp = NULL;

    (*rh_ptr)->ops = malloc(nops * sizeof(stdl_operator));
    (*rh_ptr)->ev_matrices =calloc(nops, sizeof(double*));
    (*rh_ptr)->lrvreqs = calloc(nops, sizeof(stdl_lrv_request*));

    STDL_ERROR_HANDLE_AND_REPORT((*rh_ptr)->ops == NULL || (*rh_ptr)->ev_matrices == NULL || (*rh_ptr)->lrvreqs == NULL, stdl_responses_handler_delete(*rh_ptr); return STDL_ERR_MALLOC, "malloc");

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

    for (size_t i = 0; i < rh->nlrvreqs; ++i) {
        if(rh->lrvreqs[i] != NULL)
            stdl_lrv_request_delete(rh->lrvreqs[i]);
    }

    for (size_t i = 0; i < rh->nops; ++i) {
        STDL_FREE_IF_USED(rh->ev_matrices[i]);
    }

    STDL_FREE_ALL(rh->ops, rh->ev_matrices, rh->lrvreqs, rh->eexci, rh->Xamp, rh->Yamp, rh);

    return STDL_ERR_OK;
}

int stdl_responses_handler_compute(stdl_responses_handler* rh, stdl_context* ctx) {
    assert(rh != NULL && rh->nops > 0 && ctx != NULL);

    int err;

    stdl_log_msg(0, "~~ Compute EV matrices in MO basis\n");

    for (size_t iop = 0; iop < rh->nops; ++iop) {
        size_t dim;
        err = stdl_operator_dim(rh->ops[iop], &dim);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        double* op_AO = malloc(dim * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao) * sizeof(double ));
        STDL_ERROR_HANDLE_AND_REPORT(op_AO == NULL, err = STDL_ERR_MALLOC; goto _end, "malloc");

        if(rh->ops[iop] == STDL_OP_DIPL) {
            err = stdl_basis_dsp_dipole(ctx->bs, op_AO);
            STDL_ERROR_CODE_HANDLE(err, free(op_AO); goto _end);
        } else {
            err = STDL_ERR_INPUT;
            goto _end;
        }

        // AO to MO
        rh->ev_matrices[iop] = malloc(dim * STDL_MATRIX_SP_SIZE(ctx->nmo) * sizeof(double));

        for (size_t cpt = 0; cpt < dim; ++cpt) {
            err = stdl_wavefunction_dsp_ao_to_dsp_mo(
                    ctx->original_wf->nao,
                    ctx->nmo,
                    ctx->C_ptr,
                    op_AO + cpt * STDL_MATRIX_SP_SIZE(ctx->original_wf->nao),
                    rh->ev_matrices[iop] + cpt * STDL_MATRIX_SP_SIZE(ctx->nmo)
            );
            STDL_ERROR_CODE_HANDLE(err, free(op_AO); goto _end);
        }

        STDL_FREE_IF_USED(op_AO);
    }

    if(rh->nlrvreqs > 0) {
        stdl_log_msg(0, "~~ Compute all linear response vectors\n");
        for (size_t ilrv = 0; ilrv < rh->nlrvreqs; ++ilrv) {

            stdl_lrv_request* lrvreq = rh->lrvreqs[ilrv];
            float* egrad = NULL;
            size_t dim = 0;
            err = stdl_operator_dim(lrvreq->op, &dim);
            STDL_ERROR_CODE_HANDLE(err, goto _end);

            // compute perturbed gradient
            egrad = malloc(dim * ctx->ncsfs * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(egrad == NULL,  err = STDL_ERR_MALLOC; goto _end, "malloc");

            for (size_t iop = 0; iop < rh->nops; ++iop) { // find the correct operator
                if (rh->ops[iop] == lrvreq->op) {
                    lrvreq->eta_MO = rh->ev_matrices[iop];
                    stdl_response_perturbed_gradient(ctx, dim, lrvreq->eta_MO, egrad);
                    break;
                }
            }

            // compute
            if(ctx->B == NULL)
                err = stdl_response_TDA_linear(ctx, lrvreq->nw, lrvreq->w, dim, egrad, lrvreq->X, lrvreq->Y);
            else
                err = stdl_response_TD_linear(ctx, lrvreq->nw, lrvreq->w, dim, egrad, lrvreq->X, lrvreq->Y);

            STDL_FREE_IF_USED(egrad);
            STDL_ERROR_CODE_HANDLE(err, goto _end);
        }
    }

    if(rh->nexci > 0) {
        stdl_log_msg(0, "~~ Compute amplitude vectors\n");

        if(ctx->B == NULL)
            stdl_response_TDA_casida(ctx, rh->nexci, rh->eexci, rh->Xamp);
        else
            stdl_response_TD_casida(ctx, rh->nexci, rh->eexci, rh->Xamp, rh->Yamp);
    }

    _end:
    return err;
}
