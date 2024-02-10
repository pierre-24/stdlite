#ifndef STDLITE_RESPONSE_H
#define STDLITE_RESPONSE_H

#include <stdlite/context.h>

/**
 * `ABSTOL` parameter for the precision of eigenvalue with bisect algorithm.
 * @ingroup response
 */
#ifndef STDL_RESPONSE_EIGV_ABSTOL
#define STDL_RESPONSE_EIGV_ABSTOL 1e-5f
#endif

/**
 * Solve the Casida equation to get excitation energies and amplitudes ($x^\omega$), within the Tamm-Dancoff approximation.
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @warning the `ctx->A` matrix is irreversibly modified in the process.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations requested. Must be `0 < nexci <= ctx->ncsfs`.
 * @param[out] e `float[nexci]` excitation energies (eigenvalues).
 * @param[out] X `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_TDA_casida(stdl_context *ctx, size_t nexci, float *e, float *X);

/**
 * Solve the Casida equation to get excitation energies and their corresponding response vectors ($x^\omega$, $y^\omega$).
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @warning the `ctx->A` and `ctx->B` matrices are irreversibly modified in the process.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0 && ctx->B != NULL`.
 * @param nexci number of excitations requested. Must be `0 < nexci <= ctx->ncsfs`.
 * @param[out] e `float[nexci]` excitation energies (eigenvalues).
 * @param[out] X `float[nexci,ncsfs]` response vector X for each excitation.
 * @param[out] Y `float[nexci,ncsfs]` response vector Y for each excitation.
 * @return error code
 * @ingroup response
 */
int stdl_response_RPA_casida(stdl_context *ctx, size_t nexci, float *e, float *X, float *Y);


/**
 * Create $-2\eta$, the perturbed electronic gradient matrix to be used in linear response equation.
 *
 * @param ctx a valid context
 * @param dim dimension of the expectation value `eta_MO`
 * @param eta_MO `double[dim,ctx->nmo,ctx->nmo]`, the value of $\eta$ in MO basis
 * @param[out] egrad `float[ctx->ncsfs,dim]` $-2\eta$, the resulting perturbed electronic gradient
 * @return error code
 * @ingroup response
 */
int stdl_response_perturbed_gradient(stdl_context* ctx, size_t dim, double* eta_MO, float *egrad);


/**
 * Solve the linear response equation at a frequency $\omega$ to get response vectors ($x^\omega$, $y^\omega$).
 * @param ctx a valid context
 * @param w frequency at which linear response should be computed
 * @param dim dimension of the electronic gradient
 * @param egrad `float[ncsfs,dim]` input $-2\eta$, the perturbed electronic gradient.
 * @param[out] X `float[ncsfs,dim]` response vector X for each excitation.
 * @param[out] Y `float[ncsfs,dim]` response vector Y for each excitation.
 * @return error code
 * @ingroup response
 */
int stdl_response_RPA_linear(stdl_context* ctx, float w, size_t dim, float* egrad, float* X, float* Y);

#endif //STDLITE_RESPONSE_H
