#ifndef STDLITE_RESPONSES_HANDLER_H
#define STDLITE_RESPONSES_HANDLER_H

#include "response_requests.h"

/**
 * Handler for responses
 * @ingroup responses_handler
 */
struct stdl_responses_handler_ {
    /// Number of operators
    size_t nops;

    /// `stdl_operator[nops]` list of operators
    stdl_operator* ops;

    /// `double[nops]`, <r|op|s> (in MO basis)
    double** integrals;

    /// Number of LRV requests
    size_t nlrvreqs;

    /// `stdl_lrv_request*[nlrvreqs]` Corresponding linear response vectors requests
    stdl_lrv_request** lrvreqs;

    /// Number of excited states requested
    size_t nexci;

    /// `float[nexci]` the excitation energies
    float* eexci;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf x^m$ for each excitation $\ket{m}$
    float* Xamp;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf y^m$ for each excitation $\ket{m}$, might be `NULL` if TDA
    float* Yamp;

    /// place to save responses
    char* data_output;
};

typedef struct stdl_responses_handler_ stdl_responses_handler;

/**
 * Create a new responses handler.
 *
 * @param nops Number of operators
 * @param nlrvreqs number of linear response vectors (LRV) requests
 * @param nexci number of excitations
 * @param[out] rh_ptr resulting handler
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_new(size_t nops, size_t nlrvreqs, size_t nexci, char *data_output, stdl_context *ctx, stdl_responses_handler **rh_ptr);


/**
 * Delete handler.
 *
 * @param rh a valid handler
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_delete(stdl_responses_handler* rh);

/**
 * Compute linear/amplitude vectors.
 *
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_compute(stdl_responses_handler* rh, stdl_context* ctx);

/**
 * Get the approximate space in memory
 * @param rh a valid responses handler
 * @param[out] sz the total size
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_approximate_size(stdl_responses_handler *rh, size_t nmo, size_t ncsfs, size_t *sz, size_t *ev_sz, size_t *lrv_sz, size_t *amp_sz);


#endif //STDLITE_RESPONSES_HANDLER_H
