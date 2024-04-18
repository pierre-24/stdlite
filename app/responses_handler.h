#ifndef STDLITE_RESPONSES_HANDLER_H
#define STDLITE_RESPONSES_HANDLER_H

#include "user_input_handler.h"
#include "response_requests.h"

/**
 * Handler for the LRVs associated to a given operator.
 * @ingroup responses_handler
 */
struct stdl_op_data_ {
    /// Operator
    stdl_operator op;

    /// `double[STDL_MATRIX_SP_SIZE(nmo)]`
    double* op_ints_MO;

    /// Number LRVs that are requested for `op`
    size_t nlrvs;

    /// `float[dim,ncsfs]` perturbed gradient
    float* egrad;

    /// `float[nlrvs]`, the energies
    float* w;

    /// `size_t[nlrvs]`, the index in `inp->res_w`
    size_t* iw;

    /// `float[nlrvs,nscfs,dim]` The linear response vector $\mathbf x(\omega)+\mathbf y(\omega)$
    float* XpY;

    /// `float[nlrvs,nscfs,dim]` The linear response vector $\mathbf x(\omega)-\mathbf y(\omega)$
    float* XmY;

    /// `stdl_lrv[nlrv]`, ready-to-use objects for property property_tensor
    stdl_lrv* lrvs;
};

typedef struct stdl_op_data_ stdl_op_data;

/**
 * Create a new data storage for `op`, requesting for at `nlrvs` LRVs.
 *
 * @param op operator
 * @param nlrvs number of lrvs that are requested
 * @param[out] data_ptr resulting data
 * @return error code
 * @ingroup responses_handler
 */
int stdl_op_data_new(stdl_operator op, size_t nmo, size_t ncsfs, size_t nlrvs, float *w, size_t *iw, stdl_op_data **data_ptr);

/**
 * Actually compute the linear response vectors
 *
 * @param data a valid data
 * @return error code
 * @ingroup responses_handler
 */
int stdl_op_data_compute_lrvs(stdl_op_data *data, stdl_context *ctx);

/**
 * Delete a LRV request.
 *
 * @param data a valid data
 * @return error code
 * @ingroup responses_handler
 */
int stdl_op_data_delete(stdl_op_data* data);

/**
 * Dump a LRV request.
 *
 * @param data a valid data
 * @param group_id a valid H5 group
 * @return error code
 * @ingroup responses_handler
 */
int stdl_op_data_dump_h5(stdl_op_data *data, stdl_context *ctx, hid_t group_id);

/**
 * Get the approximate space in memory
 * @param data a valid data
 * @param[out] sz the total size
 * @return error code
 * @ingroup responses_handler
 */
int stdl_op_data_approximate_size(stdl_op_data *data, size_t nmo, size_t ncsfs, size_t *sz);

/**
 * Handler for responses
 * @ingroup responses_handler
 */
struct stdl_responses_handler_ {
    /// Data associated for each operator (can be NULL)
    stdl_op_data* lrvs_data[STDL_OP_COUNT];

    /// Number of excited states requested
    size_t nexci;

    /// `float[nexci]` the excitation energies
    float* eexci;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$ for each excitation $\ket{m}$
    float* XpYamp;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf x^m+\mathbf y^m$ for each excitation $\ket{m}$, might be `NULL` if TDA
    float* XmYamp;
};

typedef struct stdl_responses_handler_ stdl_responses_handler;

/**
 * Create a new responses handler.
 *
 * @param ctx a valid context
 * @param nexci number of excitations
 * @param[out] rh_ptr resulting handler
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_new(stdl_context *ctx, size_t nexci, stdl_responses_handler **rh_ptr);

/**
 * Create a new responses handler from user input.
 *
 * @param inp a valid user input
 * @param ctx a valid context
 * @param[out] rh_ptr resulting handler
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_new_from_input(stdl_user_input_handler* inp, stdl_context* ctx, stdl_responses_handler **rh_ptr);

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
 * @param inp a valid user input
 * @param ctx a valid context
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_compute(stdl_responses_handler *rh, stdl_user_input_handler *inp, stdl_context *ctx);

/**
 * Get the approximate space in memory
 * @param rh a valid responses handler
 * @param[out] sz the total size
 * @param[out] ops_sz size associated to operators and their LRVs
 * @param[out] amp_sz size associated to amplitudes
 * @return error code
 * @ingroup responses_handler
 */
int stdl_responses_handler_approximate_size(stdl_responses_handler *rh, size_t nmo, size_t ncsfs, size_t *sz, size_t *ops_sz, size_t *amp_sz);

/**
 * Compute the property tensors corresponding to each user's request.
 *
 * @param rh a valid responses handler
 * @param inp a valid user input
 * @param ctx a valid context
 * @return error code
 * @ingroup responses_handler
 */
int stdl_response_handler_compute_properties(stdl_responses_handler* rh, stdl_user_input_handler* inp, stdl_context* ctx);


#endif //STDLITE_RESPONSES_HANDLER_H
