#ifndef STDLITE_USER_INPUT_HANDLER_H
#define STDLITE_USER_INPUT_HANDLER_H

#include <toml.h>
#include <stdlite/logging.h>
#include <stdlite/context.h>

#include "response_requests.h"

#define STDL_APP_ARG_DOC "Run a sTDA/sTD-DFT calculation.\n\n"

/**
 * Type of input for wavefunction & basis set
 * @ingroup user_input_handler
 */
enum stdl_source_type_ {
    STDL_SRC_FCHK,
    STDL_SRC_MOLDEN,
    STDL_SRC_CTX,
    STDL_SRC_CTX_WB,
};

typedef enum stdl_source_type_ stdl_source_type;

/**
 * Method to select CSFs and build A/B matrices
 * @ingroup user_input_handler
 */
enum stdl_method_ {
    /// Monopole approximation
    STDL_METHOD_MONOPOLE,
};

typedef enum stdl_method_ stdl_method;

/**
 * Gauge origin for basis set
 * @ingroup user_input_handler
 */
enum  stdl_gauge_origin_ {
    /// Origin at center of mass (default)
    STDL_GAUGE_ORIGIN_CM,

    /// Custom
    STDL_GAUGE_ORIGIN_CUSTOM,
};

typedef enum stdl_gauge_origin_ stdl_gauge_origin;


/**
 * User input structure
 * @ingroup user_input_handler
 */
struct stdl_user_input_handler_ {
    /// Title
    char* title;

    // ----- context:
    /// Source path
    char* ctx_source;

    /// Source type
    stdl_source_type ctx_source_type;

    /// H5 file output
    char* data_output;

    /// Method
    stdl_method ctx_method;

    /// Use the Tamm-Dancoff approximation
    int ctx_tda;

    /// gammaJ for monopole approximation
    float ctx_gammaJ;

    /// gammaK for monopole approximation
    float ctx_gammaK;

    /// Energy threshold for MO and CSFs
    float ctx_ethr;

    /// Perturbation energy threshold for secondary CSFs
    float ctx_e2thr;

    /// Amount of HF exchange
    float ctx_ax;

    /// Common basis origin
    stdl_gauge_origin ctx_gauge_origin;

    /// custom basis origin, valid if `ctx_gauge_origin` is `STDL_GAUGE_ORIGIN_CUSTOM`
    double ctx_gauge_origin_custom[3];

    // ----- responses:
    /// Response requests
    stdl_response_request* res_resreqs;

    /// number of frequencies
    size_t res_nw;

    /// frequencies
    float* res_w;
};

typedef struct stdl_user_input_handler_ stdl_user_input_handler;

/**
 * Create an empty user input structure with default parameters.
 *
 * @param[out] inp_ptr resulting structure
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_new(stdl_user_input_handler** inp_ptr);

/**
 * Delete handler
 *
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_delete(stdl_user_input_handler* inp);

/**
 * Change user input from options found in a TOML-formatted file.
 *
 * @param inp a valid user input structure
 * @param f path to a TOML file
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_fill_from_toml(stdl_user_input_handler* inp, FILE *f);

/**
 * Helper function to create user input directly from program input.
 *
 * @param argc number of arguments
 * @param argv arguments
 * @param[out] inp
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_new_from_args(int argc, char* argv[], stdl_user_input_handler** inp);

/**
 * Parse frequency given a a string of the form `NUMBER UNIT`, where `NUMBER` is a valid `double` and `UNIT` is either nothing (atomic units are assumed), `au`, `eV` or `nm`.
 * The result is in atomic units.
 *
 * @param input a valid `\0`-terminated string
 * @param[out] result the resulting frequency in atomic unit, if any
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_parse_frequency(char* input, double* result);

/**
 * Change user input from program arguments.
 * If the `<file>` argument is set, this implies a call to `stdl_user_input_handler_fill_from_toml()`.
 *
 * @param inp a valid user input structure
 * @param argc number of arguments
 * @param argv actual arguments
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_fill_from_args(stdl_user_input_handler* inp, int argc, char* argv[]);

/**
 * Check that the user input is correct.
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_check(stdl_user_input_handler* inp);

/**
 * Print (through `stdl_log_msg()`) the user input, as a TOML-compatible.
 *
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_log(stdl_user_input_handler* inp);

/**
 * Create context from user input.
 * @param inp a valid user input
 * @param[out] ctx_ptr context to be created
 * @return error code
 * @ingroup user_input_handler
 */
int stdl_user_input_handler_make_context(stdl_user_input_handler* inp, stdl_context **ctx_ptr);

/**
 * Get the approximate space in memory
 * @param inp a valid user input
 * @param[out] sz the total size
 * @return error code
 * @ingroup user_input_handler
 */
int
stdl_user_input_handler_approximate_size(stdl_user_input_handler *inp, size_t nexci, size_t *sz, size_t *respreq_sz);

#endif //STDLITE_USER_INPUT_HANDLER_H
