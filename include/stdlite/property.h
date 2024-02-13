#ifndef STDLITE_PROPERTY_H
#define STDLITE_PROPERTY_H

#include <stdlite/context.h>

/**
 * Compute the transition dipoles.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param X `float[nexci,ncsfs]` linear response vector
 * @param Y `float[nexci,ncsfs]` linear response vector, might be `NULL` if TDA.
 * @param[out] tdips `float[3,nexci]` the transition tdips
 * @return error code
 * @ingroup property
 */
int stdl_property_transition_dipoles(stdl_context *ctx, size_t nexci, double* dips_MO, float* X, float* Y, float * tdips);

/**
 * Print (in stdout) the excitations that were computed.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param tdips `float[3,nexci]` transition dipoles.
 * @return error code
 * @ingroup property
 */
int stdl_property_print_excitations(stdl_context *ctx, size_t nexci, float *energies, float *tdips);

/**
 * Print (in stdout) the contribution ($c^\wp_{ia} = (x^\omega_{ia})^2-(y^\omega_{ia})^2$) of each CSF to the excitations that were computed.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param X `float[nexci,ncsfs]` linear response vector
 * @param Y `float[nexci,ncsfs]` linear response vector, might be `NULL` if TDA.
 * @param thresh threshold (on the square of the coefficient) above which a contribution is printed.
 * @return error code
 * @ingroup property
 */
int stdl_property_print_excitations_contribs(stdl_context *ctx, size_t nexci, float *energies, float *X, float *Y, float thresh);

/**
 * Compute the polarizability $\alpha(-\omega;\omega)$ tensor elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param X `float[ncsfs,3]` linear response vector
 * @param Y `float[ncsfs,3]` linear response vector, might be `NULL` if TDA.
 * @param[out] alpha `float[STDL_MATRIX_SP_SIZE(3)]` the polarizability tensor
 * @return error code
 * @ingroup property
 */
int stdl_property_polarizability(stdl_context* ctx, double* dips_MO, float* X, float* Y, float* alpha);

/**
 * Compute the mean polarizability
 *
 * @param alpha `float[STDL_MATRIX_SP_SIZE(3)]` the polarizability tensor
 * @param[out] mean the mean polarizability
 * @return error code
 * @ingroup property
 */
int stdl_property_mean_polarizability(float * alpha, float* mean);

#endif //STDLITE_PROPERTY_H
