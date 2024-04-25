#include "stdlite/property_tensor.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/utils/permutations.h"
#include "stdlite/helpers.h"

typedef struct _l_elm_ {
    stdl_lrv* lrv;
    size_t cpt;
} _l_elm;

int _property_tensor_linear_element(size_t components[2], stdl_context* ctx, stdl_lrv* lrvs[2], float* value) {

    size_t nvirt = ctx->nmo - ctx->nocc;

    float v = 0;
    size_t dim1 = STDL_OPERATOR_DIM[lrvs[1]->op];

    // time-reversal symmetry of the whole thing dictates which matrix is chosen
    float* XxY = STDL_OPERATOR_TRS[lrvs[0]->op] * STDL_OPERATOR_TRS[lrvs[1]->op] > 0 ? lrvs[1]->XpYw: lrvs[1]->XmYw;

    #pragma omp parallel for reduction(+:v)
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

        float s2, d1;

        d1 = (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)] * (STDL_OPERATOR_ISSYM[lrvs[0]->op]? 1.f: -1.f);
        s2 = XxY[lia * dim1 + components[1]];

        v -= 2 * d1 * s2;
    }

    *value = v;
    return STDL_ERR_OK;
}

int stdl_property_tensor_linear(stdl_context *ctx, stdl_lrv *lrvs[2], float *tensor) {
    assert(ctx != NULL && lrvs != NULL && tensor != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute linear response tensor >");
    stdl_log_msg(1, "\n  | Preparing ");

    stdl_permutations* set;
    int err;

    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0]->op], dim1 = STDL_OPERATOR_DIM[lrvs[1]->op];

    int* isset = calloc(dim0 * dim1, sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT(isset == NULL, return STDL_ERR_MALLOC, "malloc");

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        for (size_t sigma = 0; sigma < dim1; ++sigma) {
            if(isset[zeta * dim0 + sigma])
                continue;

            stdl_log_msg(0, "-");
            stdl_log_msg(1, "\n  | Computing (%d,%d) ", zeta, sigma);
            float value = .0f;

            err = _property_tensor_linear_element((size_t[]) {zeta, sigma}, ctx, lrvs, &value);
            STDL_ERROR_CODE_HANDLE(err, STDL_FREE_IF_USED(isset); return err);

            err = stdl_permutations_new((_l_elm[]) {
                    {lrvs[0], zeta},
                    {lrvs[1], sigma},
            }, 2, sizeof(_l_elm), &set);
            STDL_ERROR_CODE_HANDLE(err, STDL_FREE_IF_USED(isset); return err);

            stdl_permutations* element = set;
            while (element != NULL) {
                _l_elm* e0 = (_l_elm *) element->perm, *e1 = e0 + 1;

                if(e0->lrv == lrvs[0] && e1->lrv == lrvs[1]) {

                    STDL_DEBUG("Set element (%ld,%ld) via permutation", e0->cpt, e1->cpt);

                    isset[e0->cpt * dim0 + e1->cpt] = 1;
                    tensor[e0->cpt * dim0 + e1->cpt] = value;
                }

                element = element->next;
            }

            stdl_permutations_delete(set);
        }
    }

    stdl_log_msg(0, "< done\n");

    STDL_FREE_IF_USED(isset);

    return STDL_ERR_OK;
}

int stdl_property_tensor_g2e_moments(stdl_context *ctx, stdl_operator ops[2], double* ops_ints_MO[2], size_t nexci, float* XpYamp, float* XmYamp, float * tg2e) {
    assert(ctx != NULL && nexci > 0 && ops_ints_MO != NULL && XpYamp != NULL && tg2e != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute ground to excited transition moments for `%s` and `%s` >", STDL_OPERATOR_NAME[ops[0]], STDL_OPERATOR_NAME[ops[1]]);
    stdl_log_msg(1, "\n  | Looping through CSFs ");

    size_t nvirt = ctx->nmo - ctx->nocc, dim0 = STDL_OPERATOR_DIM[ops[0]], dim1 = STDL_OPERATOR_DIM[ops[1]];
    float s2 = sqrtf(2);

    #pragma omp parallel for
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        for (size_t cpt = 0; cpt < dim0 + dim1; ++cpt) {
            tg2e[cpt * nexci + iexci] = .0f;
        }

        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
            float amplitude0 = XpYamp[iexci * ctx->ncsfs + lia], amplitude1 = amplitude0;

            if(XmYamp != NULL) {
                amplitude0 = (STDL_OPERATOR_TRS[ops[0]] > 0 ? XpYamp : XmYamp)[iexci * ctx->ncsfs + lia];
                amplitude1 = (STDL_OPERATOR_TRS[ops[1]] > 0 ? XpYamp : XmYamp)[iexci * ctx->ncsfs + lia];
            }

            for (size_t cpt = 0; cpt < dim0; ++cpt)
                tg2e[cpt * nexci + iexci] += s2 * amplitude0 * ((float) ops_ints_MO[0][cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]) * (STDL_OPERATOR_ISSYM[ops[0]] ? 1.f : -1.f);

            for (size_t cpt = 0; cpt < dim1; ++cpt)
                tg2e[(dim0 + cpt) * nexci + iexci] += s2 * amplitude1 * ((float) ops_ints_MO[1][cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]) * (STDL_OPERATOR_ISSYM[ops[1]] ? 1.f : -1.f);
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}

int _property_tensor_quadratic_element(size_t components[3], stdl_context* ctx, stdl_lrv* lrvs[3], float* value) {
    size_t nvirt = ctx->nmo - ctx->nocc;
    *value = 0;

    float AKp = .0f, AKm = .0f, XpYPpQ = .0f, XmYPmQ = .0f;
    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0]->op], dim1 = STDL_OPERATOR_DIM[lrvs[1]->op], dim2 = STDL_OPERATOR_DIM[lrvs[2]->op];

    #pragma omp parallel for reduction(+:AKp) reduction(+:AKm) reduction(+:XpYPpQ) reduction(+:XmYPmQ)
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt;

        for (size_t ljb = 0; ljb < ctx->ncsfs; ++ljb) {
            size_t j = ctx->csfs[ljb] / nvirt, b = ctx->csfs[ljb] % nvirt;

            if(a == b) { // ia,ja
                AKp += -.25f * (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)] * (
                               lrvs[1]->XpYw[lia * dim1 + components[1]] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                               - lrvs[1]->XmYw[lia * dim1 + components[1]] * lrvs[2]->XmYw[ljb * dim2 + components[2]]
                               + lrvs[2]->XpYw[lia * dim2 + components[2]] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                               - lrvs[2]->XmYw[lia * dim2 + components[2]] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                       );

                AKm += .25f * (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)] * (
                               lrvs[1]->XpYw[lia * dim1 + components[1]] * lrvs[2]->XmYw[ljb * dim2 + components[2]]
                               - lrvs[1]->XmYw[lia * dim1 + components[1]] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                               + lrvs[2]->XpYw[lia * dim2 + components[2]] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                               - lrvs[2]->XmYw[lia * dim2 + components[2]] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                       );

                XpYPpQ -= lrvs[0]->XpYw[lia * dim0 + components[0]] * (
                        (float) lrvs[1]->op_ints_MO[components[1] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(j, i)] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                        + (float) lrvs[2]->op_ints_MO[components[2] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(j, i)] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                );

                XmYPmQ += lrvs[0]->XmYw[lia * dim0 + components[0]] * (
                        (float) lrvs[1]->op_ints_MO[components[1] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(j, i)] *  lrvs[2]->XmYw[ljb * dim2 + components[2]]
                        + (float) lrvs[2]->op_ints_MO[components[2] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(j, i)] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                );
            }

            if(i == j) { // ia,ib
                AKp += .25f * (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * (
                        lrvs[1]->XpYw[lia * dim1 + components[1]] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                        - lrvs[1]->XmYw[lia * dim1 + components[1]] * lrvs[2]->XmYw[ljb * dim2 + components[2]]
                        + lrvs[2]->XpYw[lia * dim2 + components[2]] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                        - lrvs[2]->XmYw[lia * dim2 + components[2]] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                );

                AKm += .25f * (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * (
                        lrvs[1]->XpYw[lia * dim1 + components[1]] * lrvs[2]->XmYw[ljb * dim2 + components[2]]
                        - lrvs[1]->XmYw[lia * dim1 + components[1]] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                        + lrvs[2]->XpYw[lia * dim2 + components[2]] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                        - lrvs[2]->XmYw[lia * dim2 + components[2]] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                );

                XpYPpQ += lrvs[0]->XpYw[lia * dim0 + components[0]] * (
                        (float) lrvs[1]->op_ints_MO[components[1] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * lrvs[2]->XpYw[ljb * dim2 + components[2]]
                        + (float) lrvs[2]->op_ints_MO[components[2] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * lrvs[1]->XpYw[ljb * dim1 + components[1]]
                );

                XmYPmQ += - lrvs[0]->XmYw[lia * dim0 + components[0]] * (
                        (float) lrvs[1]->op_ints_MO[components[1] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * lrvs[2]->XmYw[ljb * dim2 + components[2]]
                        + (float) lrvs[2]->op_ints_MO[components[2] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(a, b)] * lrvs[1]->XmYw[ljb * dim1 + components[1]]
                );
            }
        }
    }

    printf("AK+ = %.3f, AK- = %.3f\n", AKp, AKm);
    printf("(X+Y)(P+Q) = %.3f, (X-Y)(P-Q) = %.3f\n", XpYPpQ, XmYPmQ);

    *value = AKp + XpYPpQ + XmYPmQ;

    return STDL_ERR_OK;
}

int stdl_property_tensor_quadratic(stdl_context *ctx, stdl_lrv *lrvs[3], float *tensor) {
    assert(ctx != NULL && lrvs != NULL && tensor != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute quadratic response tensor >");
    stdl_log_msg(1, "\n  | Preparing ");
    stdl_permutations* set = NULL;
    int err;

    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0]->op], dim1 = STDL_OPERATOR_DIM[lrvs[1]->op], dim2 = STDL_OPERATOR_DIM[lrvs[2]->op];
    int* isset = calloc(dim0 * dim1 * dim2, sizeof(int));
    STDL_ERROR_HANDLE_AND_REPORT(isset == NULL, return STDL_ERR_MALLOC, "malloc");

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        for (size_t sigma = 0; sigma < dim1; ++sigma) {
            for (size_t tau = 0; tau < dim2; ++tau) {
                if(!isset[zeta * dim0 * dim1 + sigma * dim1 + tau]) {
                    stdl_log_msg(0, "-");
                    stdl_log_msg(1, "\n  | Computing (%d,%d,%d) ", zeta, sigma, tau);

                    float value = zeta * dim0 * dim1 + sigma * dim1 + tau;

                    err = _property_tensor_quadratic_element((size_t[]) {zeta, sigma, tau}, ctx, lrvs, &value);
                    STDL_ERROR_CODE_HANDLE(err, STDL_FREE_IF_USED(isset); return err);

                    // permute over the set of coordinates and use intrinsic permutations if any
                    err = stdl_permutations_new((_l_elm[]) {
                            {lrvs[0], zeta},
                            {lrvs[1], sigma},
                            {lrvs[2], tau},
                    }, 3, sizeof(_l_elm), &set);
                    STDL_ERROR_CODE_HANDLE(err, STDL_FREE_IF_USED(isset); return err);

                    err = stdl_permutations_remove_duplicates(set, 3, sizeof(_l_elm));
                    STDL_ERROR_CODE_HANDLE(err, return err);
                    stdl_permutations* current = set;
                    while(current != NULL) {
                        _l_elm * e0 = (_l_elm *) current->perm, *e1 = e0 + 1, *e2 = e0 + 2;
                        size_t rzeta = e0->cpt, rsigma = e1->cpt, rtau = e2->cpt;

                        if(e0->lrv == lrvs[0] && e1->lrv == lrvs[1] && e2->lrv == lrvs[2]) {
                            tensor[rzeta * dim0 * dim1 + rsigma * dim1 + rtau] = value;
                            isset[rzeta * dim0 * dim1 + rsigma * dim1 + rtau] = 1;

                            STDL_DEBUG("Set element (%ld,%ld,%ld) via permutation", rzeta, rsigma, rtau);
                        }

                        current = current->next;
                    }

                    err = stdl_permutations_delete(set);
                    STDL_ERROR_CODE_HANDLE(err, return err);
                }
            }
        }
    }

    stdl_log_msg(0, "< done\n");

    STDL_FREE_IF_USED(isset);

    return STDL_ERR_OK;
}
