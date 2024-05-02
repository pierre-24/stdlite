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
 * Print a linear response tensor $-\braket{\braket{\hat A; \hat B}}_{\omega_B}$.
 *
 * @param req a valid response request
 * @param tensor `float[dim0 * dim1]` the property_tensor
 * @param wB frequency for $\hat B$
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_linear_tensor(stdl_response_request* req, float* tensor, float wB);

/**
 * Print a quadratic response tensor  $-\braket{\braket{\hat A; \hat B, \hat C}}_{\omega_B,\omega_C}$.
 *
 * @param req a valid response request
 * @param tensor `float[dim0 * dim1 * dim2]` the property_tensor
 * @param wB frequency for $\hat B$
 * @param wC frequency for $\hat C$
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_quadratic_tensor(stdl_response_request* req, float* tensor, float wB, float wC);


/**
 * Print the first hyperpolarizability tensor, $\beta(\omega_\sigma,\omega_1,\omega_2)$ and its corresponding experimental quantities.
 *
 * @param req a valid response request
 * @param beta  `float[3,3,3]` the hyperpolarizability beta
 * @param wB first input frequency ($\omega_1$)
 * @param wC second input frequency ($\omega_2$)
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_first_hyperpolarizability(stdl_response_request* req, float* beta, float wB, float wC);


/**
 * Print the ground to excited moments for `ops`.
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param ops `stdl_operator[2]` the operators
 * @param tg2e `float[(dim0 + dim1), nexci]` transition moments for both operators
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_g2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator ops[2], size_t nexci, float *tg2e);

/**
 * Print the contribution of each CSFs (if above `thresh`).
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param thresh threshold for the contributions
 * @return error code
 * @ingroup log_property
 */
void stdl_log_property_amplitude_contributions(stdl_responses_handler *rh, stdl_context *ctx, float thresh);

/**
 * Print the excited to excited moments for `ops`.
 *
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param ops `stdl_operator[3]` the operators
 * @param tg2e `float[(dim0 + dim1 + dim2), STDL_MATRIX_SP_SIZE(nexci)]` transition moments for all operators
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_e2e_moments(stdl_responses_handler *rh, stdl_context *ctx, stdl_operator ops[3], size_t nexci, float *te2e);

#endif //STDLITE_LOG_PROPERTY_H
