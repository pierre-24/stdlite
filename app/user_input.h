#ifndef STDLITE_USER_INPUT_H
#define STDLITE_USER_INPUT_H

#include <toml.h>
#include <stdlite/logging.h>
#include <stdlite/context.h>

#include "requests.h"

#define STDL_APP_ARG_DOC "Run a sTDA/sTD-DFT calculation.\n\n"

/**
 * Type of input for wavefunction & basis set
 * @ingroup user_input
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
 * @ingroup user_input
 */
enum stdl_method_ {
    /// Monopole approximation
    STDL_METHOD_MONOPOLE,

    /// Monopole approximation, direct calculation
    STDL_METHOD_MONOPOLE_DIRECT,
};

typedef enum stdl_method_ stdl_method;


/**
 * User input structure
 * @ingroup user_input
 */
struct stdl_user_input_ {
    /// Title
    char* title;

    // ----- context:
    /// Source path
    char* ctx_source;

    /// Source type
    stdl_source_type ctx_source_type;

    /// H5 file output
    char* ctx_output;

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

    // ----- responses:
    /// Response requests
    stdl_response_request* res_resreqs;

    /// Number of operators
    size_t res_nops;

    /// `stdl_operator[res_nops]` list of operators
    stdl_operator* res_ops;

    /// Number of LRV requests
    size_t res_nlrvreq;

    /// `stdl_lrv_request*[nlrv_req]` Corresponding linear response vectors requests
    stdl_lrv_request** res_lrvreqs;

    /// Number of excited states requested
    size_t res_nexci;

    /// `float[res_nexci]` the excitation energies
    float* res_eexci;

    /// `float[res_nexci,ncsfs]` amplitude vector $\mathbf x^m$ for each excitation $\ket{m}$
    float* res_Xamp;

    /// `float[res_nexci,ncsfs]` amplitude vector $\mathbf y^m$ for each excitation $\ket{m}$, might be `NULL` if TDA
    float* res_Yamp;
};

typedef struct stdl_user_input_ stdl_user_input;

/**
 * Create an empty user input structure with default parameters.
 *
 * @param[out] inp_ptr resulting structure
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_new(stdl_user_input** inp_ptr);

/**
 * Delete user input structure
 *
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_delete(stdl_user_input* inp);

/**
 * Change user input from options found in a TOML-formatted file.
 *
 * @param inp a valid user input structure
 * @param f path to a TOML file
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_fill_from_toml(stdl_user_input* inp, FILE *f);

/**
 * Helper function to create user input directly from program input.
 *
 * @param argc number of arguments
 * @param argv arguments
 * @param[out] inp
 * @return error code
 * @ingroup app
 */
int stdl_user_input_new_from_args(int argc, char* argv[], stdl_user_input** inp);

/**
 * Parse frequency given a a string of the form `NUMBER UNIT`, where `NUMBER` is a valid `double` and `UNIT` is either nothing (atomic units are assumed), `au`, `eV` or `nm`.
 * The result is in atomic units.
 *
 * @param input a valid `\0`-terminated string
 * @param[out] result the resulting frequency in atomic unit, if any
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_parse_frequency(char* input, double* result);

/**
 * Change user input from program arguments.
 * If the `<file>` argument is set, this implies a call to `stdl_user_input_fill_from_toml()`.
 *
 * @param inp a valid user input structure
 * @param argc number of arguments
 * @param argv actual arguments
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_fill_from_args(stdl_user_input* inp, int argc, char* argv[]);

/**
 * Check that the user input is correct.
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_check(stdl_user_input* inp);

/**
 * Print (through `stdl_log_msg()`) the user input, as a TOML-compatible.
 *
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_log(stdl_user_input* inp);

/**
 * Create context from user input.
 * @param inp a valid user input
 * @param[out] ctx_ptr context to be created
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_make_context(stdl_user_input* inp, stdl_context **ctx_ptr);

/**
 * Compute responses
 * @param inp a valid user input
 * @param ctx a valid context
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_prepare_responses(stdl_user_input* inp, stdl_context * ctx);

#endif //STDLITE_USER_INPUT_H
