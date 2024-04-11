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
 * Compute the experimental quantities related to HRS electric first hyperpolarizability.
 *
 * @param beta `float[3,3,3]` the first hyperpolarizability
 * @param[out] beta2_ZZZ macroscopic HRS invariant
 * @param[out] beta2_ZXX macroscopic HRS invariant
 * @return error code
 * @ingroup qexp
 */
int stdl_qexp_first_hyperpolarizability_hrs(float beta[3][3][3], float* beta2_ZZZ, float* beta2_ZXX);

#endif //STDLITE_EXPERIMENTAL_QUANTITIES_H
