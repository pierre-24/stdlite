#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "stdlite/errors.h"
#include "stdlite/utils/fchk_parser.h"

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
int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar) {
    assert(lx != NULL && name != NULL && type != NULL && is_scalar != NULL);

    if(lx->current_tk_type != STDL_TK_ALPHA) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected FCHK section to start with ALPHA");
        return STDL_ERR_UTIL_FCHK;
    }

    // read up the 40 next characters (name)
    int err;
    char buff[41]; // 40 + '\0'
    int i = 0, last_alnum = 0;
    while (i < 40) {
        if(lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF)  {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "unexpected end in section name");
            return STDL_ERR_UTIL_FCHK;
        }

        buff[i] = lx->current_tk_value;

        if(isalnum(lx->current_tk_value))
            last_alnum = i;

        err = stdl_lexer_advance(lx, 1);
        RETURN_ON_ERROR(err);
        i++;
    }

    // read out the next 7 character (only the fourth is relevant)
    i = 0;
    while (i < 7) {
        if (i == 3) {
            if(lx->current_tk_type != STDL_TK_ALPHA || (lx->current_tk_value != 'I' && lx->current_tk_value != 'R' && lx->current_tk_value != 'C')) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expecting data type to be I/R/C");
                return STDL_ERR_UTIL_FCHK;
            }

            *type = lx->current_tk_value;
        }

        err = stdl_lexer_advance(lx, 1);
        RETURN_ON_ERROR(err);

        i++;
    }

    // time to know if it is a scalar or a vector
    if (lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "line is too short to check for scalar/vector");
        return STDL_ERR_UTIL_FCHK;
    }

    *is_scalar = lx->current_tk_type != STDL_TK_ALPHA || lx->current_tk_value != 'N';

    // just advance to the value/size
    if(!(*is_scalar)) { // skip '='
        err = stdl_lexer_eat(lx, STDL_TK_ALPHA);
        RETURN_ON_ERROR(err);

        err = stdl_lexer_eat(lx, STDL_TK_EQ);
        if(err != STDL_ERR_OK) {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected EQ after 'N'");
            return STDL_ERR_UTIL_FCHK;
        }
    }

    err = stdl_lexer_skip(lx, isspace);
    RETURN_ON_ERROR(err);

    // ok, normally we are good. Time to copy the name
    *name = malloc((last_alnum + 2) * sizeof(char));
    if(*name == NULL)
        return STDL_ERR_MALLOC;

    memcpy(*name, buff, last_alnum + 1);
    (*name)[last_alnum + 1] = '\0';

    return STDL_ERR_OK;
}


int _get_vec_sz(stdl_lexer* lx, size_t* sz) {
    assert(lx != NULL && sz != NULL);

    if (lx->current_tk_type != STDL_TK_DIGIT) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected FCHK vector to start with DIGIT");
        return STDL_ERR_UTIL_FCHK;
    }

    long sz_read = -1;
    int err;

    err = stdl_parser_get_integer(lx, &sz_read);
    RETURN_ON_ERROR(err);

    if (sz_read < 0) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "size of the vector is <0");
        return STDL_ERR_UTIL_FCHK;
    }

    err = stdl_lexer_eat(lx, STDL_TK_NL); // TODO: thus, position the cursor at the beginning of next line!
    if(err != STDL_ERR_OK)  {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "size must be followed by NL");
        return STDL_ERR_UTIL_FCHK;
    }

    *sz = (size_t) sz_read;
    return STDL_ERR_OK;
}

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
int stdl_fchk_parser_get_vector_ints(stdl_lexer* lx, size_t* sz, long **vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    RETURN_ON_ERROR(err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(long));
    if(*vector == NULL)
        return STDL_ERR_MALLOC;

    size_t i = 0;
    long x;
    while (i < *sz) {
        err = stdl_parser_get_integer(lx, &x);

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 6 == 5 && lx->current_tk_type != STDL_TK_NL) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected 6 integer per line!");
                err = STDL_ERR_UTIL_FCHK;
            } else
                err = stdl_lexer_skip_whitespace_and_nl(lx);
        }

        if(err != STDL_ERR_OK) {
            free(*vector);
            return err;
        }

        i++;
    }

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

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
int stdl_fchk_parser_get_vector_numbers(stdl_lexer* lx, size_t* sz, double** vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    RETURN_ON_ERROR(err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(double));
    if(*vector == NULL)
        return STDL_ERR_MALLOC;

    size_t i = 0;
    double x;
    while (i < *sz) {
        err = stdl_parser_get_number(lx, &x);

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 5 == 4 && lx->current_tk_type != STDL_TK_NL) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected 5 real numbers per line!");
                err = STDL_ERR_UTIL_FCHK;
            } else
                err = stdl_lexer_skip_whitespace_and_nl(lx);
        }

        if(err != STDL_ERR_OK) {
            free(*vector);
            return err;
        }

        i++;
    }

    // now, we should be at the next section
    return STDL_ERR_OK;
}

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
int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, char **out) {
    assert(lx != NULL && sz != NULL && out != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    RETURN_ON_ERROR(err);

    // allocate string
    *out = malloc(((*sz) * 12 + 1) * sizeof(char));
    if(*out == NULL)
        return STDL_ERR_MALLOC;

    size_t i = 0, j;
    while (i < *sz) {
        j = 0;
        while (j < 12) {
            (*out)[i * 12 + j] = lx->current_tk_value;
            err = stdl_lexer_advance(lx, 1);

            if(err != STDL_ERR_OK) {
                free(*out);
                return err;
            }

            j++;
        }

        if(i % 5 == 4) {
            if (lx->current_tk_type != STDL_TK_NL) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected 5 packs of 12 characters per line!");
                err = STDL_ERR_UTIL_FCHK;
            } else
                err = stdl_lexer_advance(lx, 1);
        }

        if(err != STDL_ERR_OK) {
            free(*out);
            return err;
        }

        i++;
    }

    // almost there
    if(lx->current_tk_type != STDL_TK_EOF) {
        if(lx->current_tk_type != STDL_TK_NL) {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "string should ends with NL");
            err = STDL_ERR_UTIL_FCHK;
        } else
            err =stdl_lexer_advance(lx, 1);

        if(err != STDL_ERR_OK) {
            free(*out);
            return err;
        }
    }

    (*out)[(*sz) * 12] = '\0';

    return STDL_ERR_OK;
}

/**
 * Skip the current section.
 * In practice, skip as much `NL` as required.
 * @param lx a valid lexer
 * @param type type of the value(s)
 * @param is_scalar `1` if scalar, 0 if vector
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar) {
    assert(lx != NULL);

    int err, nl_to_skip = 1;
    if(!is_scalar) {
        // fetch vector size
        size_t sz = -1;
        err = _get_vec_sz(lx, &sz);
        RETURN_ON_ERROR(err);

        switch (type) {
            case 'R':
            case 'C':
                nl_to_skip = (int) sz / 5 + ((sz % 5 == 0)? 0 : 1);
                break;
            case 'I':
                nl_to_skip = ((int) sz / 6) + ((sz % 6 == 0)? 0 : 1);
                break;
        }
    }

    int i = 0;
    while (i < nl_to_skip) {
        while (lx->current_tk_type != STDL_TK_NL && lx->current_tk_type != STDL_TK_EOF) {
            err = stdl_lexer_advance(lx, 1);
            RETURN_ON_ERROR(err);
        }

        err = stdl_lexer_eat(lx, STDL_TK_NL);
        RETURN_ON_ERROR(err);

        i++;
    }

    return STDL_ERR_OK;
}

/**
 * Skip the beginning of FCHK.
 * @param lx a valid lexer
 * @return `STDL_ERR_OK` if everything was ok, error code otherwise
 * @ingroup fchk_parser
 */
int stdl_fchk_parser_skip_begin(stdl_lexer* lx) {
    assert(lx != NULL);

    int i = 0, err;
    while (i < 2) {
        while (lx->current_tk_type != STDL_TK_NL && lx->current_tk_type != STDL_TK_EOF) {
            err = stdl_lexer_advance(lx, 1);
            RETURN_ON_ERROR(err);
        }

        if(lx->current_tk_type == STDL_TK_NL) {
            err = stdl_lexer_advance(lx, 1);
            RETURN_ON_ERROR(err);
        } else { // met EOF, so stop there!
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "FCHK is too short!");
            return STDL_ERR_UTIL_FCHK;
        }

        i++;
    }

    return STDL_ERR_OK;
}
