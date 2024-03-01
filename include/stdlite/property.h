#ifndef STDLITE_PROPERTY_H
#define STDLITE_PROPERTY_H

#include <stdlite/context.h>

/**
 * Compute the transition dipoles.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param X `float[nexci,ncsfs]` amplitude vector $\mathbf x$
 * @param Y `float[nexci,ncsfs]` amplitude vector $\mathbf y$, might be `NULL` if TDA.
 * @param[out] tdips `float[3,nexci]` the transition tdips
 * @return error code
 * @ingroup property
 */
int stdl_property_transition_dipoles(stdl_context *ctx, size_t nexci, double* dips_MO, float* X, float* Y, float * tdips);

/**
 * Compute the polarizability $\alpha(-\omega;\omega)$ tensor elements.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param X `float[ncsfs,3]` linear response vector
 * @param Y `float[ncsfs,3]` linear response vector.
 * @param[out] alpha `float[STDL_MATRIX_SP_SIZE(3)]` the polarizability tensor
 * @return error code
 * @ingroup property
 */
int stdl_property_polarizability(stdl_context* ctx, double* dips_MO, float* X, float* Y, float* alpha);

/**
 * Compute the hyperpolarizability $\beta(-\omega_\sigma;\omega_1,\omega_2)$ tensor elements.
 * Use intrinsic permutations to alleviate some costs if possible.
 *
 * @note Since permutations relies on the addresses of the response vectors, the same $\mathbf x(\omega)$ (and $\mathbf y(\omega)$) should be used if some frequencies are the same (e.g., for a static calculation, the input for `Xs` should be `(float*[]) {Xs, Xs, Xs}`).
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param Xs `float*[3]` the 3 linear response vectors corresponding to $\mathbf x(-\omega_\sigma)$, $\mathbf x(\omega_1)$, and $\mathbf x(\omega_2)$.
 * @param Ys `float*[3]` the 3 linear response vectors corresponding to $\mathbf y(-\omega_\sigma)$, $\mathbf y(\omega_1)$, and $\mathbf y(\omega_2)$.
 * @param[out] beta `float[3,3,3]` the hyperpolarizability tensor
 * @return error code
 * @ingroup property
 */
int stdl_property_first_hyperpolarizability(stdl_context* ctx, double* dips_MO, float *Xs[3], float *Ys[3], float* beta);

/**
 * Evaluate the transition dipole/fluctuation operator between two excited states, $\braket{u|\hat\mu_\zeta - \delta_{uv}\,\braket{0|\hat\mu_\zeta|0}|v}$.
 *
 * @param ctx a valid context, with `ctx->ncsfs > 0`.
 * @param nexci number of excitations computed
 * @param dips_MO `float[3,STDL_MATRIX_SP_SIZE(ctx->nmo)]`, the dipole moment matrix, **in MO basis**.
 * @param X `float[nexci,ncsfs]` amplitude vector $\mathbf x$
 * @param Y `float[nexci,ncsfs]` amplitude vector $\mathbf y$, might be `NULL` if TDA.
 * @param e2etdips `float[3,STDL_MATRIX_SP_SIZE(nexci)]` the values of the different components of the transition dipoles.
 * @return error code
 * @ingroup property
 */
int stdl_property_e2e_transition_dipoles(stdl_context* ctx, size_t nexci, double* dips_MO, float * X, float * Y, float* e2etdips);

#endif //STDLITE_PROPERTY_H
