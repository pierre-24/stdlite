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
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 * Gives all (i.e., `N=ncsfs`) possible excitations, so it might result in a big `amplitudes` array.
 *
 * @warning the `ctx->A` matrix is irreversibly modified in the process.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param[out] energies `float[ncsfs]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[ncsfs,ncsfs]` amplitudes for each excitation (eigenvectors). Note that
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA_full(stdl_context *ctx, float *energies, float *amplitudes);


/**
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @warning the `ctx->A` matrix is irreversibly modified in the process.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations requested.
 * @param[out] energies `float[nexci]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA(stdl_context *ctx, size_t nexci, float *energies, float *amplitudes);


/**
 * Print (in stdout) the excitations that were computed.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param amplitudes `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
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
 * @param amplitudes `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @param thresh threshold (on the square of the coefficient) above which a contribution is printed.
 * @return
 */
int stdl_response_print_excitations_contribs(stdl_context *ctx, size_t nexci, float* energies, float *amplitudes, float thresh);

#endif //STDLITE_RESPONSE_H
