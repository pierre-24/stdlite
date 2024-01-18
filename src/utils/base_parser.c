#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>

#include "stdlite/utils/base_parser.h"
#include "stdlite/errors.h"

void stdl_error_msg_parser(char *file, int line, stdl_lexer* lx, char *format, ...) {
    assert(file != NULL && lx != NULL && format != NULL);

    va_list arglist;

    char buff[64];

    fprintf(stderr, "ERROR (%s:%d) :: ", file, line);

    sprintf(buff, isgraph(lx->current_tk_value) ? "`%c`": "0x%x", lx->current_tk_value);

    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, " (from token@%d:%d = {type=%d, value=%s})", lx->current_line, lx->current_pos_in_line, lx->current_tk_type, buff);
    fprintf(stderr, "\n");
}

/**
 * Allocate a string, and grow it from time to time if its size (given by `sz`) gets too large.
 * Actually allocates `fac * STR_MULT` bytes, and increases `fac` by one each time the size gets too large.
 * @param str str pointer to a non-`NULL` string. Caller is responsible for free'ing it.
 * @param sz current size of said string
 * @param fac scaling factor, increase periodically. Set to 0 to initialize the string
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup base_parser
 */
int stdl_grow_string(char** str, int sz, int* fac) {
    assert(str != NULL && fac != NULL && sz >= 0 && *fac >= 0);

    char* buff;
    if(*fac == 0) {
        (*fac)++;
        buff = malloc(*fac * STR_MULT * sizeof(char));
        if (buff == NULL) {
            *str = NULL;
            return STDL_ERR_MALLOC;
        } else
            *str = buff;

    } else if(sz == *fac * STR_MULT) {
        (*fac)++;
        buff = realloc(*str, (*fac) * STR_MULT * sizeof(char));
        if (buff == NULL) {
            *str = NULL;
            return STDL_ERR_MALLOC;
        } else
            *str = buff;
    }

    return STDL_ERR_OK;
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
int stdl_parser_store_value_and_grow_string(stdl_lexer* lx, char** str, int* sz, int* fac) {
    assert(lx != NULL && str != NULL && *str != NULL && sz != NULL && *sz >= 0 && fac != NULL && *fac > 0);

    (*str)[*sz] = lx->current_tk_value;
    (*sz)++;

    int err;

    // advance lexer
    err = stdl_lexer_advance(lx, 1);
    RETURN_ON_ERROR(err);

    // grow string
    err = stdl_grow_string(str, *sz, fac);
    RETURN_ON_ERROR(err);

    return STDL_ERR_OK;
}

/**
 * Parse an integer, if any.
 * @param lx a valid lexer
 * @param[out] result the resulting integer if there was something to read
 * @return `STDL_ERR_OK` if integer was read, `STDL_ERR_UTIL_PARSER` otherwise.
 * @ingroup base_parser
 */
int stdl_parser_get_integer(stdl_lexer* lx, long *result) {
    assert(lx != NULL && result != NULL);

    if(lx->current_tk_type != STDL_TK_DIGIT && lx->current_tk_type != STDL_TK_PLUS && lx->current_tk_type != STDL_TK_DASH) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected integer to start by DIGIT|PLUS|DASH");
        return STDL_ERR_UTIL_PARSER;
    }

    int err;
    char* str;
    int sz = 0, fac = 0;

    err = stdl_grow_string(&str, sz, &fac);
    RETURN_ON_ERROR(err);

    // store first token
    err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
    if(err != STDL_ERR_OK) {
        free(str);
        return err;
    }

    // store next digits
    while(lx->current_tk_type == STDL_TK_DIGIT) {
        err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
        if(err != STDL_ERR_OK) {
            free(str);
            return err;
        }
    }

    // interpret
    str[sz] = '\0';
    char* end;
    *result = strtol(str, &end, 10);
    free(str);

    if ((int) (end-str) != sz) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "strtol() error");
        return STDL_ERR_UTIL_PARSER;
    }

    return STDL_ERR_OK;
}

/**
 * Parse a real number matching `(PLUS|DASH)? DIGIT* (DOT DIGIT*)? (('E'|'e') (PLUS|MINUS)* DIGIT*)?` (as a `double`), if any.
 * @param lx a valid lexer
 * @param[out] result the resulting real number if there was something to read
 * @return `STDL_ERR_OK` if real number was read, `STDL_ERR_UTIL_PARSER` otherwise.
 * @ingroup base_parser
 */
int stdl_parser_get_number(stdl_lexer* lx, double* result) {
    assert(lx != NULL && result != NULL);

    if(lx->current_tk_type != STDL_TK_DIGIT && lx->current_tk_type != STDL_TK_PLUS && lx->current_tk_type != STDL_TK_DASH && lx->current_tk_type != STDL_TK_DOT) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected real to start by DIGIT|PLUS|DASH|DOT");
        return STDL_ERR_UTIL_PARSER;
    }

    int err;
    char* str;
    int sz = 0, fac = 0, dot_found = lx->current_tk_type == STDL_TK_DOT, exp_found = 0;

    err = stdl_grow_string(&str, sz, &fac);
    RETURN_ON_ERROR(err);

    // store first token
    err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
    if(err != STDL_ERR_OK) {
        free(str);
        return err;
    }

    // store next tokens
    while (lx->current_tk_type == STDL_TK_DIGIT || lx->current_tk_type == STDL_TK_DOT  || lx->current_tk_type == STDL_TK_ALPHA) {
        if(lx->current_tk_type == STDL_TK_DIGIT)
            err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
        else if (lx->current_tk_type == STDL_TK_DOT) {
            if(dot_found) // that would be two dots, so break
                break;
            if(exp_found) // that would be a dot in exp, so break
                break;

            dot_found = 1;
            err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
        } else if (lx->current_tk_type == STDL_TK_ALPHA) {
            if (exp_found) // that would be two exps, break!
                break;

            if (lx->current_tk_value != 'e' && lx->current_tk_value != 'E')
                break;

            exp_found = 1;
            err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);

            if(err == STDL_ERR_OK) {
                if(lx->current_tk_type == STDL_TK_PLUS || lx->current_tk_type == STDL_TK_DASH) // check for PLUS|DASH
                    err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);

                if(lx->current_tk_type != STDL_TK_DIGIT) {
                    stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected DIGIT after exp mark");
                    err = STDL_ERR_UTIL_PARSER;
                }
            }
        }
        if(err != STDL_ERR_OK) {
            free(str);
            return err;
        }
    }

    // interpret
    str[sz] = '\0';
    char* end;
    *result = strtod(str, &end);
    free(str);

    if ((int) (end-str) != sz) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "strtod() error");
        return STDL_ERR_UTIL_PARSER;
    }

    return STDL_ERR_OK;
}

/**
 * Put the stream of tokens in `result` as long as `predicate` is true.
 * @param lx a valid lexer
 * @param predicate a predicate, called for each token, of the form `bool predicate(int c)`, where `c` is the current token value.
 * @param[out] result the resulting string. Caller is responsible for free'ing it.
 * @return `STDL_ERR_OK` if everything went well.
 * @ingroup base_parser
 */
int stdl_parser_get_literal(stdl_lexer* lx, int (*predicate)(int), char** result) {
    assert(lx != NULL && predicate != NULL && result != NULL);

    int err;
    int sz = 0, fac = 0;

    err = stdl_grow_string(result, sz, &fac);
    RETURN_ON_ERROR(err);

    while(lx->current_tk_type != STDL_TK_EOF && predicate((int) lx->current_tk_value)) {
        err = stdl_parser_store_value_and_grow_string(lx, result, &sz, &fac);
        if(err != STDL_ERR_OK) {
            free(*result);
            return err;
        }
    }

    (*result)[sz] = '\0';

    return STDL_ERR_OK;
}
