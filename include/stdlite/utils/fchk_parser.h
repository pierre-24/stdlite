#ifndef STDL_FCHK_PARSER_H
#define STDL_FCHK_PARSER_H

#include <stdlite/utils/base_parser.h>

/**
 * Skip the beginning of FCHK.
 * Skip *i.e.*, the two first lines (beginning of title section + type, method, basis),
 * since this is info one can found in other places of the file anyway.
 * @param lx a valid lexer
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_skip_begin(stdl_lexer* lx);

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
int stdl_fchk_parser_get_scalar_int(stdl_lexer* lx, long *value);

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
int stdl_fchk_parser_get_vector_ints(stdl_lexer* lx, size_t* sz, long **vector);

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

#endif //STDL_FCHK_PARSER_H
