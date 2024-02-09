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
 * Solve the Casida equation to get excitation energies and responses vector ($x^\omega$, $y^\omega$).
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
 * Print (in stdout) the excitations that were computed.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param amplitudes `float[nexci,ncsfs]` amplitudes for each excitation: either `X` (TDA) or `X+Y` (TD).
 * @param dipoles_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @return error code
 * @ingroup response
 */
int stdl_response_print_excitations(stdl_context *ctx, size_t nexci, float *energies, float *amplitudes, double* dipoles_MO);

/**
 * Print (in stdout) the contribution of each CSF to the excitations that were computed.
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param amplitudes `float[nexci,ncsfs]` amplitudes for each excitation: either `X` (TDA) or `X+Y` (TD).
 * @param thresh threshold (on the square of the coefficient) above which a contribution is printed.
 * @return
 */
int stdl_response_print_excitations_contribs(stdl_context *ctx, size_t nexci, float* energies, float *amplitudes, float thresh);

#endif //STDLITE_RESPONSE_H
