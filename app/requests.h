#ifndef STDLITE_REQUESTS_H
#define STDLITE_REQUESTS_H


#include <stdlite/context.h>

/**
 * Operators for the linear responses
 * @ingroup requests
 */
enum stdl_operator_ {
    /// Dipole length operator
    STDL_OP_DIPL,
};

typedef enum stdl_operator_ stdl_operator;


/**
 * Linear response vectors (LRV) requests, in an orderly manner. Also store linear vectors when computed.
 * @ingroup requests
 */
struct stdl_lrv_request_ {
    /// Operator
    stdl_operator op;

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
 * @param[out] req_prt resulting request
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_new(stdl_operator op, size_t nw, stdl_lrv_request** req_prt);

/**
 * Actually compute the requested linear response vectors
 *
 * @param req a valid request
 * @param op_MO `double[dim,ctx->nmo,ctx->nmo]` elements of $\braket{p|op|q}$.
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_compute(stdl_lrv_request* req, stdl_context* ctx, double* op_MO);

/**
 * Delete a LRV request.
 *
 * @param req a valid request
 * @return error code
 * @ingroup requests
 */
int stdl_lrv_request_delete(stdl_lrv_request* req);

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

    /// `stdl_operator[resp_order-res_order+1]` the operators associated with the request
    stdl_operator *ops;

    /// `float[resp_order-res_order]` the energies at which linear response vectors should be computed
    float* w;

    /// number of amplitude vectors (*i.e.*, excitations) requested, any negatives number means "all"
    int nroot;

    /// `stdl_linear_response_vectors_request[resp_order-res_order]` for each operator (except the first one), the corresponding linear response vector request
    struct stdl_lrv_request_* requests;

    /// `size_t[resp_order-res_order]` For each frequency, its corresponding position in the linear response vector request
    size_t* wpos;

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
 * @param w energies at which linear response vectors shoudl be computed
 * @param nroot number of excitations (and, thus, amplitude vectors) that should be computed. Negative number means "all"
 * @param[out] req_ptr the resulting request
 * @return error code
 * @ingroup requests
 */
int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroot, stdl_response_request** req_ptr);

/**
 * Delete a request.
 *
 * @param req a valid request
 * @return error code
 * @ingroup requests
 */
int stdl_response_request_delete(stdl_response_request* req);


#endif //STDLITE_REQUESTS_H