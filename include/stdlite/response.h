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
 * Gives all possible excitations, so it might result in a big `amplitudes` array.
 *
 * @param ctx a valid context
 * @param ncsfs number of CSFs that are considered (or, in other words, the length of `A`).
 * @param A `float[STLD_MATRIX_SP_SIZE(ncsfs)]`, an approximate electronic Hessian matrix.
 * @param[out] energies `float[ncsfs]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[ncsfs,ncsfs]` amplitudes for each excitation (eigenvectors). Note that
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA_full(stdl_context *ctx, size_t ncsfs, float *A, float **energies, float **amplitudes);


/**
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 * Only returns the first `nexci` first excitation energies. Works well if `nexci << ncsfs`.
 * The precision on the eigenvalues is given by `STDL_RESPONSE_EIGV_ABSTOL`.
 *
 * @param ctx a valid context
 * @param ncsfs number of CSFs that are considered (or, in other words, the length of `A`).
 * @param A `float[STLD_MATRIX_SP_SIZE(ncsfs)]`, an approximate electronic Hessian matrix.
 * @param nexci number of excitations requested.
 * @param[out] energies `float[nexci]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[ncsfs,nexci]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA(stdl_context* ctx, size_t ncsfs, float *A, size_t nexci, float** energies, float** amplitudes);

#endif //STDLITE_RESPONSE_H
