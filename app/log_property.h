#ifndef STDLITE_LOG_PROPERTY_H
#define STDLITE_LOG_PROPERTY_H

#include "response_requests.h"
#include "responses_handler.h"

/**
 * Print the polarizability tensor and its corresponding experimental quantities.
 *
 * @param req a valid response request
 * @param alpha `float[6]` the polarizability tensor
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_polarizability(stdl_response_request* req, float* alpha);

/**
 * Print the first hyperpolarizability tensor and its corresponding experimental quantities.
 *
 * @param req a valid response request
 * @param beta  `float[3,3,3]` the hyperpolarizability tensor
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_first_hyperpolarizability(stdl_response_request* req, float beta[3][3][3]);

/**
 * Print the ground to excited dipole moments, as well as the contributin of each CSFs (if above `thresh`, to a given dipole moment)
 * @param rh a valid responses handler
 * @param ctx a valid context
 * @param tdips transition dipoles
 * @param thresh threshold for the contributions
 * @return error code
 * @ingroup log_property
 */
int stdl_log_property_g2e_dipoles(stdl_responses_handler *rh, stdl_context *ctx, float *tdips, float thresh);

#endif //STDLITE_LOG_PROPERTY_H
