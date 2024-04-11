#ifndef STDLITE_RESPONSE_REQUESTS_H
#define STDLITE_RESPONSE_REQUESTS_H


#include <stdlite/context.h>
#include <stdlite/integrals.h>
#include <stdlite/property_tensor.h>

/**
 * Any response request from the TOML file. It is a chained list.
 *
 * @ingroup requests
 */
struct stdl_response_request_ {
    /// order of the response (1 = linear, 2 = quadratic, etc)
    size_t resp_order;

    /// number of residues (0 = none, 1 = first residue, 2 = second residue)
    size_t res_order;

    /// Number of operators
    size_t nops; // = resp_order-res_order+1

    /// `stdl_operator[nops]` the operators associated with the request
    stdl_operator *ops;

    /// Number of frequencies and LRVs
    size_t nlrvs; // == (resp_order == res_order)? 0: resp_order-res_order+1

    /// `size_t[nlrvs]` the energies at which linear response vectors should be computed
    size_t* iw;

    /// number of amplitude vectors (*i.e.*, excitations) requested, a negative value means "all"
    int nroots;

    /// Tensor. Its shape depends on the property
    float* property_tensor;

    /// Next request
    struct stdl_response_request_* next;
};

typedef struct stdl_response_request_ stdl_response_request;

/**
 * Create a new response request.
 *
 * @param resp_order Order of the response, must be > 0.
 * @param res_order Order of the residue, must be `res_order < resp_order`
 * @param ops operators
 * @param iw index of energies at which linear response vectors should be computed
 * @param nroots number of excitations (and, thus, amplitude vectors) that should be computed. Negative number means "all"
 * @param[out] req_ptr the resulting request
 * @return error code
 * @ingroup requests
 */
int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, size_t *iw, int nroots, stdl_response_request** req_ptr);

/**
 * Delete a request.
 *
 * @param req a valid request
 * @return error code
 * @ingroup requests
 */
int stdl_response_request_delete(stdl_response_request* req);

/**
 * Get the approximate space in memory
 * @param req a valid request
 * @param[out] sz the total size
 * @return error code
 * @ingroup requests
 */
int stdl_response_request_approximate_size(stdl_response_request* req, size_t* sz);


#endif //STDLITE_RESPONSE_REQUESTS_H
