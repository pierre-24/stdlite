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
int stdl_responses_handler_new(size_t nops, size_t nlrvreqs, size_t nexci, stdl_context *ctx,  stdl_responses_handler **rh_ptr);


/**
 * Delete handler.
 *
 * @param rh a valid handler
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_delete(stdl_responses_handler* rh);


#endif //STDLITE_RESPONSES_HANDLER_H
