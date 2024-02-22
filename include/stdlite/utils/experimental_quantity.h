#ifndef STDLITE_EXPERIMENTAL_QUANTITIES_H
#define STDLITE_EXPERIMENTAL_QUANTITIES_H

#include <stdlite/context.h>

/**
 * Compute the experimental quantities related to electric polarizability: iso- and anisotropy.
 *
 * @param alpha `float[STDL_MATRIX_SP_SIZE(3)]` the polarizability tensor
 * @param[out] iso the isotropic polarizability value
 * @param[out] aniso the anisotropic polarizability value
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_polarizability(float * alpha, float* iso, float* aniso);

/**
 * Compute the experimental quantities related to HRS electric first hyperpolarizability.
 *
 * @param beta `float[3,3,3]` the first hyperpolarizability
 * @param[out] beta2_ZZZ macroscopic HRS invariant
 * @param[out] beta2_ZXX macroscopic HRS invariant
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_first_hyperpolarizability_hrs(float beta[3][3][3], float* beta2_ZZZ, float* beta2_ZXX);


/**
 * Print (in stdout) the excitations that were computed.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param energies `float[nexci]` excitation energies (eigenvalues).
 * @param tdips `float[3,nexci]` transition dipoles.
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_excitations_print(stdl_context *ctx, size_t nexci, float *energies, float *tdips);

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
 * @ingroup qexp
 */
int stdl_qexp_excitations_contribs_print(stdl_context *ctx, size_t nexci, float *energies, float *X, float *Y, float thresh);

#endif //STDLITE_EXPERIMENTAL_QUANTITIES_H
