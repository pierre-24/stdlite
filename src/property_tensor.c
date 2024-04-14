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
    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0]->op], dim1 = STDL_OPERATOR_DIM[lrvs[1]->op];

    // time-reversal symmetry of the whole thing
    float trs = STDL_OPERATOR_TRS[lrvs[0]->op] * STDL_OPERATOR_TRS[lrvs[1]->op];

    #pragma omp parallel for reduction(+:v)
    for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
        size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

        float s1, s2, d1, d2;

        d1 = (float) lrvs[0]->op_ints_MO[components[0] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)] * (STDL_OPERATOR_ISSYM[lrvs[0]->op]? 1.f: -1.f);
        d2 = (float) lrvs[1]->op_ints_MO[components[1] * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)] * (STDL_OPERATOR_ISSYM[lrvs[1]->op]? 1.f: -1.f);

        s1 = lrvs[0]->Xw[lia * dim0 + components[0]] + trs * lrvs[0]->Yw[lia * dim0 + components[0]];
        s2 = lrvs[1]->Xw[lia * dim1 + components[1]] + trs * lrvs[1]->Yw[lia * dim1 + components[1]];

        v -= d1 * s2 + d2 * s1;
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

                    STDL_DEBUG("Set element (%ld,%ld) of property_tensor via permutation", e0->cpt, e1->cpt);

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

int stdl_property_tensor_g2e_moments(stdl_context *ctx, stdl_operator op, double* op_ints_MO, size_t nexci, float* Xamp, float* Yamp, float * tg2e) {
    assert(ctx != NULL && nexci > 0 && op_ints_MO != NULL && Xamp != NULL && tg2e != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute ground to excited transition moments for `%s` >", STDL_OPERATOR_NAME[op]);
    stdl_log_msg(1, "\n  | Looping through CSFs ");

    size_t nvirt = ctx->nmo - ctx->nocc, dim0 = STDL_OPERATOR_DIM[op];
    float s2 = sqrtf(2);

    #pragma omp parallel for
    for (size_t iexci = 0; iexci < nexci; ++iexci) {
        tg2e[0 * nexci + iexci] = tg2e[1 * nexci + iexci] = tg2e[2 * nexci + iexci] = .0f;

        for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
            size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;
            float amplitude = Xamp[iexci * ctx->ncsfs + lia];
            if(Yamp != NULL)
                amplitude += STDL_OPERATOR_TRS[op] * Yamp[iexci * ctx->ncsfs + lia];

            for (size_t cpt = 0; cpt < dim0; ++cpt)
                tg2e[cpt * nexci + iexci] += s2 * amplitude * ((float) op_ints_MO[cpt * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)]) * (STDL_OPERATOR_ISSYM[op] ? 1.f : -1.f);
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}
