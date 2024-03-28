#ifndef STDLITE_RESPONSE_H
#define STDLITE_RESPONSE_H

#include <stdlite/context.h>

/**
 * `ABSTOL` parameter for the precision of eigenvalue in bisect algorithms used by LAPACK.
 * @ingroup response
 */
#ifndef STDL_RESPONSE_EIGV_ABSTOL
#define STDL_RESPONSE_EIGV_ABSTOL 1e-5f
#endif

/**
 * Solve the Casida equation to get excitation energies and amplitudes vectors ($\mathbf x^m$) within the Tamm-Dancoff approximation.
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations requested. Must be `0 < nexci <= ctx->ncsfs`.
 * @param[out] e `float[nexci]` excitation energies $\omega_m$ .
 * @param[out] X `float[nexci,ncsfs]` amplitudes vector, $\mathbf x^m$, for each excitation $\ket{m}$.
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_casida(stdl_context *ctx, size_t nexci, float *e, float *X);

/**
 * Solve the Casida equation to get excitation energies and their corresponding amplitude vectors ($\mathbf x^m$, $\mathbf y^m$).
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nexci number of excitations requested. Must be `0 < nexci <= ctx->ncsfs`.
 * @param[out] e `float[nexci]` excitation energies.
 * @param[out] X `float[nexci,ncsfs]` amplitude vector $\mathbf x^m$ for each excitation $\ket{m}$.
 * @param[out] Y `float[nexci,ncsfs]` amplitude vector $\mathbf y^m$ for each excitation $\ket{m}$.
 * @return error code
 * @ingroup response
 */
int stdl_response_TD_casida(stdl_context *ctx, size_t nexci, float *e, float *X, float *Y);


/**
 * Create $-2\eta$, the perturbed electronic gradient matrix to be used in linear response equation (`stdl_response_TD_linear()`).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dim dimension of the expectation value `eta_MO`
 * @param eta_MO `double[dim,ctx->nmo,ctx->nmo]`, the value of $\eta$ in MO basis
 * @param[out] egrad `float[ctx->ncsfs,dim]` $-2\eta$, the resulting perturbed electronic gradient
 * @return error code
 * @ingroup response
 */
int stdl_response_perturbed_gradient(stdl_context* ctx, size_t dim, double* eta_MO, float *egrad);


/**
 * Solve the linear response equation at `nw` energies $\{\omega_i\}$  to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nw]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] X `float[nw,ncsfs,ndim]` response vector $\mathbf x(\omega)$, in each dimension.
 * @param[out] Y `float[nw,ncsfs,ndim]` response vector $\mathbf y(\omega)$, in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TD_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, float *egrad, float *X, float *Y);

/**
 * Solve the linear response equation at `nw` energies $\{\omega_i\}$, within the Tamm-Dancoff approximation, to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nw]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] X `float[nw,ncsfs,ndim]` response vector $\mathbf x(\omega)$, in each dimension.
 * @param[out] Y `float[nw,ncsfs,ndim]` response vector $\mathbf y(\omega)$ in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, float *egrad, float *X, float *Y);


/**
 * Compute the linear response function $\pm\braket{\braket{\hat A; \hat B}}_\omega$ tensor elements.
 * Assume that $\hat A$ is Hermitian.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dims `size_t[2]` dimensionality of $\hat A$ and $\hat B$
 * @param A_elmt_MO `float[dims[0],STDL_MATRIX_SP_SIZE(ctx->nmo)]`, elements of $\braket{p|\hat A|q}$
 * @param X `float[ncsfs,dims[1]]` linear response vector $\mathbf x(\omega)$ for $\hat B$
 * @param Y `float[ncsfs,dims[1]]` linear response vector $\mathbf y(\omega)$ for $\hat B$
 * @param get_rf if `0`, the tensor element are $-\braket{\braket{\hat A; \hat B}}_\omega$, the sign is reversed otherwise
 * @param[out] tensor`float[dims[0], dims[1]]` the resulting tensor
 * @return error code
 * @ingroup response
 */
int stdl_response_lr_tensor(stdl_context *ctx, size_t dims[2], double *A_elmt_MO, float *X, float *Y, int get_rf, float *tensor);


#endif //STDLITE_RESPONSE_H
