#include "stdlite/property_tensor.h"
#include "stdlite/utils/matrix.h"
#include "stdlite/utils/permutations.h"
#include "stdlite/helpers.h"

#include <omp.h>

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
    stdl_log_msg(0, "Compute %ld ground to excited transition moments for `%s` and `%s` >", nexci, STDL_OPERATOR_NAME[ops[0]], STDL_OPERATOR_NAME[ops[1]]);
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

    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0]->op], zeta = components[0];

    // permute over B and C operators
    stdl_permutations* set = NULL;
    int err = stdl_permutations_new((_l_elm[]) {
            {lrvs[1], components[1]},
            {lrvs[2], components[2]},
    }, 2, sizeof(_l_elm), &set);
    STDL_ERROR_CODE_HANDLE(err, return err);

    STDL_ERROR_CODE_HANDLE(err, return err);
    stdl_permutations* current = set;
    while(current != NULL) {
        _l_elm * e1 = (_l_elm *) current->perm, *e2 = e1 + 1;
        size_t sigma = e1->cpt, tau = e2->cpt, dim1 = STDL_OPERATOR_DIM[e1->lrv->op], dim2 = STDL_OPERATOR_DIM[e2->lrv->op];

        float Mp = .0f, Np = .0f;

        #pragma omp parallel for reduction(+:Mp) reduction(+:Np)
        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt;

            float Ta = lrvs[0]->XpYw[lia * dim0 + zeta], Ua = lrvs[0]->XmYw[lia * dim0 + zeta], Tb = e1->lrv->XpYw[lia * dim1 + sigma], Ub = e1->lrv->XmYw[lia * dim1 + sigma];

            for (size_t ljb = 0; ljb < ctx->ncsfs; ++ljb) {
                size_t j = ctx->csfs[ljb] / nvirt, b = ctx->csfs[ljb] % nvirt;
                float Tc = e2->lrv->XpYw[ljb * dim2 + tau], Uc = e2->lrv->XmYw[ljb * dim2 + tau];

                if(a == b) { // ia,ja → A
                    Mp += .25f * (float) lrvs[0]->op_ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)] * (!STDL_OPERATOR_ISSYM[lrvs[0]->op] && i < j ? -1.f : 1.f) * (STDL_OPERATOR_ISSYM[lrvs[0]->op] ? Ub * Uc - Tb * Tc : Ub * Tc - Tb * Uc);
                    Mp -= .5f * (float) e1->lrv->op_ints_MO[sigma * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)] * (!STDL_OPERATOR_ISSYM[e1->lrv->op] && i < j ? -1.f : 1.f) * (Ta * (STDL_OPERATOR_ISSYM[e1->lrv->op] ? Tc : -Uc) + Ua * (STDL_OPERATOR_ISSYM[e1->lrv->op] ? Uc : -Tc));
                }

                if(i == j) { // ia,ib → B
                    Np += .25f * (float) lrvs[0]->op_ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(ctx->nocc + a, ctx->nocc + b)] * (!STDL_OPERATOR_ISSYM[lrvs[0]->op] && a < b ? -1.f : 1.f) * (STDL_OPERATOR_ISSYM[lrvs[0]->op] ? Tb * Tc - Ub * Uc : Ub * Tc - Tb * Uc);
                    Np += .5f * (float) e1->lrv->op_ints_MO[sigma * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(ctx->nocc + a, ctx->nocc + b)] * (!STDL_OPERATOR_ISSYM[e1->lrv->op] && a < b ? -1.f : 1.f) * (Ta * (STDL_OPERATOR_ISSYM[e1->lrv->op] ? Tc : -Uc) + Ua * (STDL_OPERATOR_ISSYM[e1->lrv->op] ? Uc : -Tc));
                }
            }
        }

        (*value) += Mp + Np;
        current = current->next;
    }

    err = stdl_permutations_delete(set);
    STDL_ERROR_CODE_HANDLE(err, return err);

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


int stdl_property_tensor_e2e_moments(stdl_context *ctx, stdl_operator ops[3], double* ops_ints_MO[3], size_t nexci, float* XpYamp, float* XmYamp, float * te2e) {
    assert(ops != NULL && ops_ints_MO != NULL && nexci > 0 &&  XpYamp != NULL && te2e != NULL);

    size_t nvirt = ctx->nmo - ctx->nocc, dims[3] = {STDL_OPERATOR_DIM[ops[0]], STDL_OPERATOR_DIM[ops[1]], STDL_OPERATOR_DIM[ops[2]]};
    size_t nthreads = 1;

    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
            nthreads = omp_get_num_threads();
    }

    size_t sz = (dims[0] + dims[1] + dims[2]) * nthreads;
    double wrk_asz;
    char* wrk_usz;
    stdl_convert_size(sz * sizeof(float), &wrk_asz, &wrk_usz);
    stdl_log_msg(0, "Extra memory required (components): %.1f%s\n", wrk_asz, wrk_usz);

    float * buff = calloc(sz, sizeof(float ));
    STDL_ERROR_HANDLE_AND_REPORT(buff == NULL, return STDL_ERR_MALLOC, "malloc");

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute %ld excited to excited transition moments for `%s`, `%s`, and `%s` >", STDL_MATRIX_SP_SIZE(nexci), STDL_OPERATOR_NAME[ops[0]], STDL_OPERATOR_NAME[ops[1]], STDL_OPERATOR_NAME[ops[2]]);
    stdl_log_msg(1, "\n  | Looping through CSFs ");

    #pragma omp parallel for schedule(guided)
    for (size_t m = 0; m < nexci; ++m) {
        float* tmp = buff + omp_get_thread_num() * (dims[0] + dims[1] + dims[2]);

        for (size_t n = 0; n <= m; ++n) {

            for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
                size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt;

                float tm_ia = XpYamp[m * ctx->ncsfs + lia];
                float tn_ia = XpYamp[n * ctx->ncsfs + lia];
                float um_ia = tm_ia, un_ia = tn_ia;

                if(XmYamp != NULL) {
                    um_ia = XmYamp[m * ctx->ncsfs + lia];
                    un_ia = XmYamp[n * ctx->ncsfs + lia];
                }

                for (size_t ljb = 0; ljb < ctx->ncsfs; ++ljb) {
                    size_t j = ctx->csfs[ljb] / nvirt, b = ctx->csfs[ljb] % nvirt;

                    float tm_jb = XpYamp[m * ctx->ncsfs + ljb];
                    float tn_jb = XpYamp[n * ctx->ncsfs + ljb];
                    float um_jb = tm_jb, un_jb = tn_jb;

                    if(XmYamp != NULL) {
                        um_jb = XmYamp[m * ctx->ncsfs + ljb];
                        un_jb = XmYamp[n * ctx->ncsfs + ljb];
                    }

                    float amplitude = tm_ia * tn_jb + um_ia * un_jb + tn_ia * tm_jb + un_ia * um_jb;

                    if(b == a) { // ia,ja
                        size_t shift = 0;
                        for (size_t iop = 0; iop < 3; ++iop) {
                            for (size_t cpt = 0; cpt < dims[iop]; ++cpt)
                                tmp[shift + cpt] -= (float) ops_ints_MO[iop][cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, j)] * amplitude;

                            shift += dims[iop];
                        }
                    }

                    if(j == i) {// ia,ib
                        size_t shift = 0;
                        for (size_t iop = 0; iop < 3; ++iop) {
                            for (size_t cpt = 0; cpt < dims[iop]; ++cpt) {
                                tmp[shift + cpt] += (float) ops_ints_MO[iop][cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(ctx->nocc + a,ctx->nocc + b)] * amplitude;
                            }

                            shift += dims[iop];
                        }
                    }
                }
            }

            for (size_t cpt = 0; cpt < (dims[0] + dims[1] + dims[2]); ++cpt) {
                te2e[cpt * STDL_MATRIX_SP_SIZE(nexci) + STDL_MATRIX_SP_IDX(m, n)] = .125f * tmp[cpt];
                tmp[cpt] = .0f;
            }
        }
    }

    STDL_FREE_IF_USED(buff);

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}
