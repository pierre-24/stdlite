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
 * @param[out] Xamp `float[nexci,ncsfs]` amplitudes vector, $\mathbf x^m$, for each excitation $\ket{m}$.
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_casida(stdl_context *ctx, size_t nexci, float *e, float *Xamp);

/**
 * Solve the Casida equation to get excitation energies and their corresponding amplitude vectors ($\mathbf x^m$, $\mathbf y^m$).
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nexci number of excitations requested. Must be `0 < nexci <= ctx->ncsfs`.
 * @param[out] e `float[nexci]` excitation energies.
 * @param[out] XpYamp `float[nexci,ncsfs]` $\mathbf x^m+\mathbf y^m$ for each excitation $\ket{m}$.
 * @param[out] XmYamp `float[nexci,ncsfs]` $\mathbf x^m - \mathbf y^m$ for each excitation $\ket{m}$.
 * @return error code
 * @ingroup response
 */
int stdl_response_TD_casida(stdl_context *ctx, size_t nexci, float *e, float *XpYamp, float *XmYamp);


/**
 * Create the perturbed electronic gradient matrix, $-2\eta$.
 * to be used in linear response equation (`stdl_response_TD_linear()`).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dim dimension of the corresponding operator
 * @param issym whether the operator is symmetric (`1`, $\eta_{pq} = \eta_{qp}$) or not (`0`, $\eta_{pq} = -\eta_{qp}$)
 * @param op_ints_MO `double[dim,ctx->nmo,ctx->nmo]`, the value of $\eta$ in MO basis
 * @param[out] egrad `float[ctx->ncsfs,dim]` $-2\eta$, the resulting perturbed electronic gradient
 * @return error code
 * @ingroup response
 */
int stdl_response_perturbed_gradient(stdl_context *ctx, size_t dim, int issym, double *op_ints_MO, float *egrad);


/**
 * Solve the linear response equation at `nlrvs` energies $\{\omega_i\}$  to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nlrvs]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param isherm whether the operator is hermitian (`1`) or anti-hermitian (`0`)
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] XpY `float[nlrvs,ncsfs,ndim]` $\mathbf x(\omega)+\mathbf y(\omega)$, in each dimension.
 * @param[out] XmY `float[nlrvs,ncsfs,ndim]` $\mathbf x(\omega)-\mathbf y(\omega)$ in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TD_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int isherm, float *egrad, float *XpY, float *XmY);

/**
 * Solve the linear response equation at `nlrvs` energies $\{\omega_i\}$, within the Tamm-Dancoff approximation, to get the corresponding response vectors ($\mathbf x(\omega_i)$, $\mathbf y(\omega_i)$).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nw number of energies at which linear response should be computed
 * @param w `float[nlrvs]` energies at which linear response should be computed
 * @param ndim dimension of the electronic gradient
 * @param isherm whether the operator is hermitian (`1`) or anti-hermitian (`0`)
 * @param egrad `float[ncsfs,ndim]` $-2\eta$, the perturbed electronic gradient in each dimension.
 * @param[out] XpY `float[nlrvs,ncsfs,ndim]` $\mathbf x(\omega)+\mathbf y(\omega)$, in each dimension.
 * @param[out] XmY `float[nlrvs,ncsfs,ndim]` $\mathbf x(\omega)-\mathbf y(\omega)$ in each dimension.
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_linear(stdl_context *ctx, size_t nw, float *w, size_t ndim, int isherm, float *egrad, float *XpY, float *XmY);


#endif //STDLITE_RESPONSE_H
