#include "stdlite/property_tensor.h"
#include "stdlite/utils/matrix.h"

int stdl_property_tensor_linear(stdl_context *ctx, stdl_lrv lrvs[2], float *tensor) {
    assert(ctx != NULL && lrvs != NULL && tensor != NULL);

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Compute linear response tensor >");
    stdl_log_msg(1, "\n  | Preparing ");

    size_t nvirt = ctx->nmo - ctx->nocc;

    float s1, s2, d1, d2;
    size_t dim0 = STDL_OPERATOR_DIM[lrvs[0].op], dim1 = STDL_OPERATOR_DIM[lrvs[1].op];

    for (size_t zeta = 0; zeta < dim0; ++zeta) {
        for (size_t sigma = 0; sigma < dim1; ++sigma) {

            stdl_log_msg(0, "-");
            stdl_log_msg(1, "\n  | Computing (%d1,%d1) ", zeta, sigma);
            float value = .0f;

            #pragma omp parallel for reduction(+:value) private(d1, d2, s1, s2)
            for (size_t lia = 0; lia < ctx->ncsfs; ++lia) {
                size_t i = ctx->csfs[lia] / nvirt, a = ctx->csfs[lia] % nvirt + ctx->nocc;

                d1 = (float) lrvs[0].op_ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];
                d2 = (float) lrvs[1].op_ints_MO[zeta * STDL_MATRIX_SP_SIZE(ctx->nmo) + STDL_MATRIX_SP_IDX(i, a)];

                s1 = lrvs[0].Xw[lia * dim1 + sigma] + (STDL_OPERATOR_HERMITIAN[lrvs[0].op]? 1.f: -1.f) * lrvs[0].Yw[lia * dim1 + sigma];
                s2 = lrvs[1].Xw[lia * dim1 + sigma] + (STDL_OPERATOR_HERMITIAN[lrvs[1].op]? 1.f: -1.f) * lrvs[1].Yw[lia * dim1 + sigma];

                value -= d1 * s1 + d2 * s2;
            }

            tensor[zeta * dim0 + sigma] = value;
        }
    }

    stdl_log_msg(0, "< done\n");

    return STDL_ERR_OK;
}