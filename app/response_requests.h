#ifndef STDLITE_RESPONSE_REQUESTS_H
#define STDLITE_RESPONSE_REQUESTS_H


#include <stdlite/context.h>
#include <stdlite/integrals.h>
#include <stdlite/property_tensor.h>

/**
 * Linear response vectors (LRV) requests, in an orderly manner. Also store linear vectors when computed.
 * @ingroup requests
 */
struct stdl_lrv_request_ {
    /// Operator
    stdl_operator op;

    /// pointer to the MO representation of `op`
    double* op_integrals;

    /// `float[]` perturbed gradient
    float* egrad;

    /// Number of energies at which the vectors should be computed
    size_t nw;

    /// `float[nw]`, the energies
    float* w;

    /// `float[nw,nscfs,dim]` The linear response vector $\mathbf x(\omega)$
    float* X;

    /// `float[nw,nscfs,dim]` The linear response vector $\mathbf y(\omega)$
    float* Y;
};

typedef struct stdl_lrv_request_ stdl_lrv_request;

/**
 * Create a new linear response vectors request for `op` at `nw` energies.
 *
 * @param op operator
 * @param nw number of energies
 * @param[out] req_ptr resulting request
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_new(stdl_operator op, size_t nw, size_t ncsfs, stdl_lrv_request **req_ptr);

/**
 * Actually compute the requested linear response vectors
 *
 * @param lrvreq a valid request
 * @param op_MO `double[dim,ctx->nmo,ctx->nmo]` elements of $\braket{p|op|q}$.
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_compute(stdl_lrv_request *lrvreq, stdl_context *ctx);

/**
 * Delete a LRV request.
 *
 * @param req a valid request
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_delete(stdl_lrv_request* req);

/**
 * Dump a LRV request.
 *
 * @param req a valid request
 * @param group_id a valid H5 group
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_dump_h5(stdl_lrv_request *req, stdl_context *ctx, hid_t group_id);

/**
 * Get the approximate space in memory
 * @param req a valid request
 * @param[out] sz the total size
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_approximate_size(stdl_lrv_request *req, size_t ncsfs, size_t *sz);

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

    /// Number of frequencies (and LRVs)
    size_t nw; // == (resp_order == res_order)? 0: resp_order-res_order+1

    /// `size_t[nw]` the energies at which linear response vectors should be computed
    size_t* iw;

    /// `stdl_lrv[nw]` for the corresponding linear response vector
    stdl_lrv* lrvs;

    /// number of amplitude vectors (*i.e.*, excitations) requested, a negative value means "all"
    int nroots;

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
