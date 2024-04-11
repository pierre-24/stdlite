#ifndef STDLITE_TENSOR_H
#define STDLITE_TENSOR_H

#include <stdlite/context.h>
#include <stdlite/integrals.h>

/**
 * A linear response vector.
 * @warning this is just a convenient object to compute property_tensor, one should still free the components.
 * @ingroup property_tensor
 */
struct stdl_lrv_ {
    /// The operator used to compute LRV
    stdl_operator op;

    /// The values of `<p|op|q>`
    double* op_ints_MO;

    /// energy a which the LRV was computed
    float w;

    /// `float[dim,ncsfs]` $\mathbf x_{\zeta,ia}(\omega)$
    float* Xw;

    /// `float[dim,ncsfs]` $\mathbf y_{\zeta,ia}(\omega)$
    float* Yw;
};

typedef struct stdl_lrv_ stdl_lrv;

/**
 * Compute the linear response function tensor $$-\braket{\braket{\hat A; \hat B}}_\omega$ elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param lrvs `stdl_lrv[2]` the LRV (one for $\hat A$, one for $\hat B$) from which the property_tensor will be computed
 * @param[out] tensor`float[STDL_OPERATOR_DIM[lrvs[0].op], STDL_OPERATOR_DIM[lrvs[1].op]]` the resulting property_tensor
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_linear(stdl_context *ctx, stdl_lrv *lrvs[2], float *tensor);


/**
 * Compute ground to excited transition moments for `op`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param op the operator
 * @param op_ints_MO `float[STDL_OPERATOR_DIM[op],STDL_MATRIX_SP_SIZE(ctx->nmo)]`, operator integrals
 * @param nexci number of excitations computed
 * @param X `float[nexci,ncsfs]` amplitude vector $\mathbf x$
 * @param Y `float[nexci,ncsfs]` amplitude vector $\mathbf y$, might be `NULL` if TDA.
 * @param[out] tdips `float[STDL_OPERATOR_DIM[op]nexci]` the resultings transition moments
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_g2e_moments(stdl_context *ctx, stdl_operator op, double* op_ints_MO, size_t nexci, float* Xamp, float* Yamp, float * tg2e);


#endif //STDLITE_TENSOR_H
