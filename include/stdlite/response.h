#ifndef STDLITE_RESPONSE_H
#define STDLITE_RESPONSE_H

#include <stdlite/context.h>

/**
 * Solve the Casida equation to get excitation energies and amplitudes, within the Tamm-Dancoff approximation.
 *
 * @param ctx a valid context
 * @param ncsfs number of csfs that are considered (or, in other word, length of `A`).
 * @param A `float[ncsfs,ncsfs]`, an approximate electronic Hessian matrix.
 * @param[out] nexci number of resulting excitation
 * @param[out]  energies `float[nexci]` excitation energies (eigenvalues).
 * @param[out] amplitudes `float[nexci,ncsfs]` amplitudes for each excitation (eigenvectors).
 * @return error code
 * @ingroup response
 */
int stdl_response_casida_TDA(stdl_context* ctx, size_t ncsfs, float *A, size_t* nexci, float** energies, float** amplitudes);


#endif //STDLITE_RESPONSE_H
