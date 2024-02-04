#ifndef STDL_BASE_PARSER_H
#define STDL_BASE_PARSER_H

#include <stdlite/utils/lexer.h>

/**
 * Multiplier for string growth.
 * @ingroup base_parser
 * */
#define STDL_STR_MULT 128

/**
 * Allocate a string, and grow it from time to time if its size (given by `sz`) gets too large.
 * Actually allocates `fac * STDL_STR_MULT` bytes, and increases `fac` by one each time the size gets too large.
 * @param[in,out] str_ptr str_ptr pointer to a non-`NULL` string. Caller is responsible for free'ing it.
 * @param sz current size of said string
 * @param fac scaling factor, increase periodically. Set to 0 to initialize the string
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup base_parser
 */
int stdl_grow_string(char** str_ptr, int sz, int* fac);

/**
 * Output a specific error message indicating the current token if `DEBUG_LVL` is above or equal to 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param lx a valid lexer
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup base_parser
 */
void stdl_error_msg_parser(char *file, int line, stdl_lexer* lx, char *format, ...);

/**
 * Copycat of `STDL_ERROR_HANDLE_AND_REPORT`, but using the `stdl_error_msg_parser` function instead.
 * @ingroup base_parser
 */
#define STDL_LEXER_ERROR_HAR(lx, assertion, error_action, ...)       \
    if(assertion) {                                                  \
        stdl_error_msg_parser(__FILE__, __LINE__, lx, __VA_ARGS__);  \
        {error_action;}                                              \
    }

/**
 * Store the current token value in a string, increase its size by 1, then grow it.
 * Also advance the lexer to the next token
 * @param lx a valid lexer
 * @param str str pointer to a non-`NULL` string. Caller is responsible for free'ing it.
 * @param sz size of said string.
 * @param fac scaling factor, increase periodically.
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup base_parser
 */
int stdl_parser_store_value_and_grow_string(stdl_lexer* lx, char** str, int* sz, int* fac);

/**
* Parse an integer, if any.
* @param lx a valid lexer
* @param[out] result the resulting integer if there was something to read
* @return `STDL_ERR_OK` if integer was read, `STDL_ERR_UTIL_PARSER` otherwise.
* @ingroup base_parser
*/
int stdl_parser_get_integer(stdl_lexer* lx, long *result);

/**
 * Parse a real number matching `(PLUS|DASH)? DIGIT* (DOT DIGIT*)? (('E'|'e') (PLUS|MINUS)* DIGIT*)?` (as a `double`), if any.
 * @param lx a valid lexer
 * @param[out] result the resulting real number if there was something to read
 * @return `STDL_ERR_OK` if real number was read, `STDL_ERR_UTIL_PARSER` otherwise.
 * @ingroup base_parser
 */
int stdl_parser_get_number(stdl_lexer* lx, double* result);

/**
 * Put the stream of tokens in `result` as long as `predicate` is true.
 * @param lx a valid lexer
 * @param predicate a predicate, called for each token, of the form `int predicate(int c)`, where `c` is the current token value.
 * @param[out] result the resulting string. Caller is responsible for free'ing it.
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup base_parser
 */
int stdl_parser_get_literal(stdl_lexer* lx, int (*predicate)(int), char** result);

#endif //STDL_BASE_PARSER_H
