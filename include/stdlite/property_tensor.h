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

    /// `float[dim,ncsfs]` $\mathbf x_{\zeta,ia}(\omega)+\mathbf y_{\zeta,ia}(\omega)$
    float* XpYw;

    /// `float[dim,ncsfs]` $\mathbf x_{\zeta,ia}(\omega)-\mathbf y_{\zeta,ia}(\omega)$
    float* XmYw;
};

typedef struct stdl_lrv_ stdl_lrv;

/**
 * Compute the linear response function tensor $$-\braket{\braket{\hat A; \hat B}}_\omega$ elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param lrvs `stdl_lrv[2]` the LRV (one for $\hat A$, one for $\hat B$) from which the tensor will be computed
 * @param[out] tensor`float[STDL_OPERATOR_DIM[lrvs[0].op], STDL_OPERATOR_DIM[lrvs[1].op]]` the resulting tensor
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_linear(stdl_context *ctx, stdl_lrv *lrvs[2], float *tensor);


/**
 * Compute ground to excited transition moments for `ops`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param ops the operators
 * @param ops_ints_MO `float[dim,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, operator integrals
 * @param nexci number of excitations computed
 * @param XpYamp `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$
 * @param XmYamp `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$, might be `NULL` if TDA.
 * @param[out] tg2e `float[dim0 + dim1, nexci]` the resulting transition moments
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_g2e_moments(stdl_context *ctx, stdl_operator ops[2], double* ops_ints_MO[2], size_t nexci, float* XpYamp, float* XmYamp, float * tg2e);

/**
 * Compute the quadratic response function tensor $$-\braket{\braket{\hat A; \hat B, \hat C}}_{\omega_B,\omega_C}$ elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param lrvs `stdl_lrv[3]` the LRV (one for $\hat A$, one for $\hat B$, and one for $\hat C$) from which the property tensor will be computed
 * @param[out] tensor`float[STDL_OPERATOR_DIM[lrvs[0].op], STDL_OPERATOR_DIM[lrvs[1].op], STDL_OPERATOR_DIM[lrvs[2].op]]` the resulting tensor
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_quadratic(stdl_context *ctx, stdl_lrv *lrvs[3], float *tensor);

/**
 * Compute excited to excited transition moments for `ops`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param ops the operators
 * @param ops_ints_MO `float[dim,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, operator integrals
 * @param nexci number of excitations computed
 * @param XpYamp `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$
 * @param XmYamp `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$, might be `NULL` if TDA.
 * @param[out] te2e `float[dim0 + dim1 + dim2, STDL_MATRIX_SP_SIZE(nexci)]` the resulting transition moments
 * @return error code
 * @ingroup property_tensor
 */
int stdl_property_tensor_e2e_moments(stdl_context *ctx, stdl_operator ops[3], double* ops_ints_MO[3], size_t nexci, float* XpYamp, float* XmYamp, float * te2e);


#endif //STDLITE_TENSOR_H
