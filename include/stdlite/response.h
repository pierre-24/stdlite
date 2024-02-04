#ifndef STDLITE_RESPONSE_H
#define STDLITE_RESPONSE_H

#include <stdlite/context.h>

/**
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 * Gives all possible excitations, so it might result in a big `amplitudes` array.
 *
 * @param ctx a valid context
 * @param ncsfs number of csfs that are considered (or, in other word, length of `A`).
 * @param A `float[ncsfs,ncsfs]`, an approximate electronic Hessian matrix.
 * @param[out]  energies `float[ncsfs]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[ncsfs,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA_full(stdl_context *ctx, size_t ncsfs, float *A, float **energies, float **amplitudes);


/**
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 * Only returns the subset of exictations energies (i.e., eigenvalues) that are below a given energy threshold.
 *
 * @param ctx a valid context
 * @param ncsfs number of csfs that are considered (or, in other word, length of `A`).
 * @param ethr energy
 * @param A `float[ncsfs,ncsfs]`, an approximate electronic Hessian matrix.
 * @param[out] nexci number of resulting excitation
 * @param[out]  energies `float[nexci]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA(stdl_context* ctx, size_t ncsfs, float ethr, float *A, size_t* nexci, float** energies, float** amplitudes);

#endif //STDLITE_RESPONSE_H
