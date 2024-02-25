#ifndef STDL_FCHK_PARSER_H
#define STDL_FCHK_PARSER_H

#include <stdlite/utils/base_parser.h>
#include <stdlite/wavefunction.h>
#include <stdlite/basis.h>

/**
 * Skip the beginning of FCHK.
 * Skip *i.e.*, the two first lines (beginning of title section + type, method, basis),
 * since this is info one can found in other places of the file anyway.
 * @param lx a valid lexer
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_skip_intro(stdl_lexer* lx);

/**
 * Parse the info of a section in FCHK.
 * Stops at the beginning of the value (if scalar) or of the size (if vector).
 * @param lx a valid lexer
 * @param[out] name the name of the section
 * @param[out] type the type of section. Valid outputs are `I`, `R`, and `C`.
 * @param[out] is_scalar `1` if the section is a scalar
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar);

/**
 * Parse a scalar integer.
 * @param lx a valid lexer
 * @param[out] value the value, if any.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_scalar_integer(stdl_lexer* lx, long *value);

/**
 * Parse a scalar real number.
 * @param lx a valid lexer
 * @param[out] value the value, if any.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_scalar_number(stdl_lexer* lx, double* value);

/**
 * Parse a vector of integers in FCHK.
 * Expect the lexer to have stopped right before the size of the vector.
 * Read lines in format `6I12`.
 * @param lx a valid lexer
 * @param[out] sz size of the vector
 * @param[out] vector the vector, if any. Caller is responsible for free'ing it.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_vector_integers(stdl_lexer* lx, size_t* sz, long **vector);

/**
 * Parse a vector of integers in FCHK.
 * Expect the lexer to have stopped right before the size of the vector.
 * Read lines in format `6I12`.
 *
 * @note Immediate version: `*vector` must be have been allocated.
 * This version is to be used when one knows in advance the shape of the data (because, *e.g.*, it was given by another section) and tries to avoid copy.
 *
 * @param lx a valid lexer
 * @param sz expected size of the vector, will be checked.
 * @param[in,out] vector the vector, which will be filled.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_vector_integers_immediate(stdl_lexer* lx, size_t sz, long **vector);

/**
 * Parse a vector of real numbers in FCHK.
 * Expect the lexer to have stopped right before the size of the vector.
 * Read lines in format `5E16.8`.
 * @param lx a valid lexer
 * @param[out] sz size of the vector
 * @param[out] vector the vector, if any. Caller is responsible for free'ing it.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_vector_numbers(stdl_lexer* lx, size_t* sz, double** vector);

/**
 * Parse a vector of real numbers in FCHK.
 * Expect the lexer to have stopped right before the size of the vector.
 * Read lines in format `5E16.8`.
 *
 * @note Immediate version: `*vector` must be have been allocated.
 * This version is to be used when one knows in advance the shape of the data (because, *e.g.*, it was given by another section) and tries to avoid copy.
 *
 * @param lx a valid lexer
 * @param sz expected size of the vector, will be checked.
 * @param[in,out] vector the vector, which will be filled.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_vector_numbers_immediate(stdl_lexer* lx, size_t sz, double** vector);

/**
 * Parse a vector of strings in FCHK and merge everything in one (long) string.
 * Expect the lexer to have stopped right before the size of the vector.
 * Read lines in format `5A12`.
 * @param lx a valid lexer
 * @param[out] sz size of the vector
 * @param[out] out the string, if any. Caller is responsible for free'ing it.
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, char **out);

/**
 * Skip the current section.
 * In practice, skip as much `NL` as required.
 * @param lx a valid lexer
 * @param type type of the value(s)
 * @param is_scalar `1` if scalar, 0 if vector
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar);


/**
 * Extract a wavefunction (`stdl_wavefunction`) and a basis set (`stdl_basis`) from a FCHK file.
 * Expects that the introduction has been skipped (with `stdl_fchk_parser_skip_intro()`).
 * Expects that the section comes in the order in which Gaussian prints them.
 * The FCHK is read up to the "Alpha MO coefficients" section.
 * Then, it computes the overlap matrix (`S`) as this is not available in the FCHK.
 * @param lx a valid lexer, opened on a FCHK
 * @param[out] wf_ptr a wavefunction to be created
 * @param[out] bs_ptr a basis set to be created
 * @return if everything went well, `STDL_ERR_OK`.
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_extract(stdl_lexer *lx, stdl_wavefunction **wf_ptr, stdl_basis **bs_ptr);


/// Structure that holds basis set data in a format that resemble the one used by Gaussian in its FCHK.
struct stdl_basis_data_ {
    /// number of basis functions
    size_t nbas;

    /// number of primitives, `nprim >= nbas`
    size_t nprim;

    /// `long[nbas]`, the angular moment of each basis function.
    /// Use the Gaussian specification, so `0=s, 1=p, -1=sp, 2=6d, -2=5d, 3=10f, -3=7f` (and so all).
    long *bas_types;

    /// `long[nbas]` the number of primitive in each basis function
    long *prims_per_bas; // [nbas]

    /// `long[nbas]`, 1-based list of correspondence between basis function and atom
    long *bastoatm;

    /// `double[3 * nprim]`, list of exponents, coefs, ans S=P coefs, as `[e0, e1, ..., eN, c0, c1, ..., cN, cp0, cp1, ..., cpN]`.
    double *benv;
};

typedef struct stdl_basis_data_ stdl_basis_data;

/**
 * Create a new basis set data holder.
 *
 * @param nbas Number of basis function, must be >0
 * @param nprims Number of primitives, must be `nprims >= nbas`
 * @param[out] dt_ptr Resulting data
 * @return error code
 */
int stdl_basis_data_new(size_t nbas, size_t nprims, stdl_basis_data **dt_ptr);

/**
 * Delete a basis set data holder.
 * @param dt a valid pointer to data
 * @return error code
 */
int stdl_basis_data_delete(stdl_basis_data* dt);

/**
 * Convert a basis set data to an actual `stdl_basis`.
 *
 * @param dt a valid pointer to data
 * @param natm number of atoms, must be > 0.
 * @param atm `double[4*natm]` list of atoms with nuclear charge (0) and coordinates (1:3)
 * @param[out] bs_ptr Resulting basis set.
 * @return
 */
int stdl_basis_data_to_basis(stdl_basis_data *dt, size_t natm, double *atm, stdl_basis **bs_ptr);

#endif //STDL_FCHK_PARSER_H
