#ifndef STDLITE_EXPERIMENTAL_QUANTITIES_H
#define STDLITE_EXPERIMENTAL_QUANTITIES_H

#include <stdlite/context.h>

/**
 * Compute the experimental quantities related to electric polarizability: iso- and anisotropy.
 *
 * @param alpha `float[STDL_MATRIX_SP_SIZE(3)]` the polarizability property_tensor
 * @param[out] iso the isotropic polarizability value
 * @param[out] aniso the anisotropic polarizability value
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_polarizability(float * alpha, float* iso, float* aniso);


/**
 * Compute macroscopic quantities associated with first hyperpolarizability.
 *
 * @param beta `float[3,3,3]` the first hyperpolarizability
 * @param[out] beta_vec the beta vector.
 * @return
 */
int stdl_qexp_first_hyperpolarizability(float beta[27], float beta_vec[3]);

/**
 * Compute the experimental macroscopic quantities related to HRS electric first hyperpolarizability:
 *
 * + $\braket{\beta^2_{ZZZ}}$,
 * + $\braket{\beta^2_{ZXX}}$,
 * + $|\beta_{J=1}|^2$ (assuming Kleinman conditions), and
 * + $|\beta_{J=3}|^2$ (assuming Kleinman conditions).
 *
 * @param beta `float[3,3,3]` the first hyperpolarizability
 * @param[out] beta2_ZZZ macroscopic HRS invariant
 * @param[out] beta2_ZXX macroscopic HRS invariant
 * @param[out] beta2_J1 macroscopic HRS rotational invariant
 * @param[out] beta2_J3 macroscopic HRS rotational invariant
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_first_hyperpolarizability_hrs(float beta[27], float *beta2_ZZZ, float *beta2_ZXX, float *beta2_J1, float *beta2_J3);

#endif //STDLITE_EXPERIMENTAL_QUANTITIES_H
