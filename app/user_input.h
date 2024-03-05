#ifndef STDLITE_USER_INPUT_H
#define STDLITE_USER_INPUT_H

#include <toml.h>
#include <stdlite/logging.h>

/**
 * Type of input for wavefunction & basis set
 * @ingroup user_input
 */
enum stdl_source_type_ {
    STDL_SRC_FCHK,
    STDL_SRC_MOLDEN,
    STDL_SRC_CTX,
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
    char* ctx_source_path;

    /// Source type
    stdl_source_type ctx_source_type;

    /// H5 file output
    char* ctx_output;

    /// Method
    stdl_method ctx_method;

    /// Use the Tamm-Dancoff approximation
    int ctx_use_tda;

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
 * @param inp a valid user input structure
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_delete(stdl_user_input* inp);

/**
 * Change user input from options found in a TOML-formatted file
 * @param inp a valid user input structure
 * @param path path to a TOML file
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_fill_from_toml(stdl_user_input* inp, char* path);

/**
 * Parse frequency given a a string of the form `NUMBER UNIT`, where `NUMBER` is a valid `double` and `UNIT` is either `au`, `eV` or `nm`.
 * The result is in atomic units.
 *
 * @param input a valid `\0`-terminated string
 * @param[out] result the resulting frequency in atomic unit, if any
 * @return error code
 * @ingroup user_input
 */
int stdl_user_input_parse_frequency(char* input, double* result);

#endif //STDLITE_USER_INPUT_H
