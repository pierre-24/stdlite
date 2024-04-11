#ifndef STDLITE_LOG_PROPERTY_H
#define STDLITE_LOG_PROPERTY_H

#include "response_requests.h"
#include "responses_handler.h"

/**
 * Print the polarizability property_tensor and its corresponding experimental quantities.
 *
 * @param req a valid response request
 * @param alpha `float[9]` the polarizability property_tensor
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_polarizability(stdl_response_request* req, float* alpha, float w);

/**
 * Print a linear response property_tensor.
 *
 * @param req a valid response request
 * @param tensor `float[dim0 * dim1]` the property_tensor
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_linear_tensor(stdl_response_request* req, float* tensor, float w);


/**
 * Print the first hyperpolarizability property_tensor and its corresponding experimental quantities.
 *
 * @param req a valid response request
 * @param beta  `float[3,3,3]` the hyperpolarizability property_tensor
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_first_hyperpolarizability(stdl_response_request* req, float beta[3][3][3]);

/**
 * Print the ground to excited dipole moments, as well as the contribution of each CSFs (if above `thresh`)
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param tdips transition dipoles
 * @param thresh threshold for the contributions
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_g2e_dipoles(stdl_responses_handler *rh, stdl_context *ctx, float *tdips, float thresh);

/**
 * Print the ground to excited moments for `op`, as well as the contribution of each CSFs (if above `thresh`)
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param op the operator
 * @param tg2e transition dipoles
 * @param thresh threshold for the contributions
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_g2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator op, float *tg2e, float thresh);

#endif //STDLITE_LOG_PROPERTY_H
