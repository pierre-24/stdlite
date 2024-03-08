#ifndef STDLITE_USER_INPUT_H
#define STDLITE_USER_INPUT_H

#include <toml.h>
#include <stdlite/logging.h>
#include <stdlite/context.h>

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
 * Operators for the linear responses
 * @ingroup user_input
 */
enum stdl_operator_ {
    /// Dipole length operator
    STDL_OP_DIPL,
};

typedef enum stdl_operator_ stdl_operator;


/**
 * Linear response vectors (LRV) requests, in an orderly manner. Also store linear vectors when computed.
 * @ingroup user_input
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
 * @param op operator
 * @param nw number of energies
 * @param[out] req_prt resulting request
 * @return error code
 * @ingroup user_input
 */
int stdl_lrv_request_new(stdl_operator op, size_t nw, stdl_lrv_request** req_prt);

/**
 * Actually compute the requested linear response vectors
 * @param req a valid request
 * @param op_MO `double[dim,ctx->nmo,ctx->nmo]` elements of $\braket{p|op|q}$.
 * @return error code
 * @ingroup user_input
 */
int stdl_lrv_request_compute(stdl_lrv_request* req, stdl_context* ctx, double* op_MO);

/**
 * Delete a LRV request.
 * @param req a valid request
 * @return error code
 * @ingroup user_input
 */
int stdl_lrv_request_delete(stdl_lrv_request* req);

/**
 * Any response request from the TOML file. It is a chained list.
 * @ingroup user_input
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

     struct stdl_response_request_* next;
 };

 typedef struct stdl_response_request_ stdl_response_request;

/**
 * Create a new response request.
 * @param resp_order Order of the response, must be > 0.
 * @param res_order Order of the residue, must be `res_order < resp_order`
 * @param ops operators
 * @param w energies at which linear response vectors shoudl be computed
 * @param nroot number of excitations (and, thus, amplitude vectors) that should be computed. Negative number means "all"
 * @param[out] req_ptr the resulting request
 * @return error code
 * @ingroup user_input
 */
int stdl_response_request_new(size_t resp_order, size_t res_order, stdl_operator* ops, float* w, int nroot, stdl_response_request** req_ptr);

/**
 * Delete a request
 * @param req a valid request
 * @return error code
 * @ingroup user_input
 */
int stdl_response_request_delete(stdl_response_request* req);

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

    /// Response requests
    stdl_response_request* res_reqs;

    /// Number of LRV requests
    size_t nlrv_req;

    /// `stdl_lrv_request[nlrv_req]` Corresponding linear response vectors requests
    stdl_lrv_request* lrv_reqs;

    /// Number of excited states requested
    size_t nexci;

    /// `float[nexci]` the excitation energies
    float* eexci;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf x^m$ for each excitation $\ket{m}$
    float* Xamp;

    /// `float[nexci,ncsfs]` amplitude vector $\mathbf y^m$ for each excitation $\ket{m}$, might be `NULL` if TDA
    float* Yamp;
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

#endif //STDLITE_USER_INPUT_H
