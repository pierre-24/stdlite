#ifndef STDLITE_TENSOR_H
#define STDLITE_TENSOR_H

#include <stdlite/context.h>
#include <stdlite/integrals.h>

/**
 * A linear response vector.
 * @warning this is just a convenient object to compute tensor, one should still free the components.
 * @ingroup tensor
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
 * Compute the linear response function $\pm\braket{\braket{\hat A; \hat B}}_\omega$ tensor elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param lrvs `stdl_lrv[2]` the LRV (one for $\hat A$, one for $\hat B$) from which the tensor will be computed
 * @param[out] tensor`float[STDL_OPERATOR_DIM[lrvs[0].op], STDL_OPERATOR_DIM[lrvs[1].op]]` the resulting tensor
 * @return error code
 * @ingroup tensor
 */
int stdl_property_tensor_linear(stdl_context *ctx, stdl_lrv lrvs[2], float *tensor);


#endif //STDLITE_TENSOR_H
