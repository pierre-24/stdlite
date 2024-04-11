#ifndef STDLITE_RESPONSE_H
#define STDLITE_RESPONSE_H

#include <stdlite/context.h>

/**
 * `ABSTOL` parameter for the precision of eigenvalue in bisect algorithms used by LAPACK.
 * @ingroup response
 */
#ifndef STDL_RESPONSE_EIGV_ABSTOL
#define STDL_RESPONSE_EIGV_ABSTOL 1e-6f
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
 * Create the perturbed electronic gradient matrix, $-2\eta$ (or $-2\Im(\eta)$ if `is_hermitian=0`)
 * to be used in linear response equation (`stdl_response_TD_linear()`).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dim dimension of the corresponding operator
 * @param is_hermitian whether the operator is hermitian (`1`) or not (`0`)
 * @param op_ints_MO `double[dim,ctx->nmo,ctx->nmo]`, the value of $\eta$ in MO basis
 * @param[out] egrad `float[ctx->ncsfs,dim]` $-2\eta$, the resulting perturbed electronic gradient
 * @return error code
 * @ingroup response
 */
int stdl_response_perturbed_gradient(stdl_context *ctx, size_t dim, int is_hermitian, double *op_ints_MO, float *egrad);


/**
 * Solve the linear response equation at `nlrvs` energies $\{\omega_i\}$  to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nlrvs]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param is_hermitian whether the operator is hermitian (`1`) or not (`0`)
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] X `float[nlrvs,ncsfs,ndim]` response vector $\mathbf x(\omega)$, in each dimension.
 * @param[out] Y `float[nlrvs,ncsfs,ndim]` response vector $\mathbf y(\omega)$, in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TD_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int is_hermitian, float *egrad, float *X, float *Y);

/**
 * Solve the linear response equation at `nlrvs` energies $\{\omega_i\}$, within the Tamm-Dancoff approximation, to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nlrvs]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param is_hermitian whether the operator is hermitian (`1`) or not (`0`)
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] X `float[nlrvs,ncsfs,ndim]` response vector $\mathbf x(\omega)$, in each dimension.
 * @param[out] Y `float[nlrvs,ncsfs,ndim]` response vector $\mathbf y(\omega)$ in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int is_hermitian, float *egrad, float *X, float *Y);


#endif //STDLITE_RESPONSE_H
