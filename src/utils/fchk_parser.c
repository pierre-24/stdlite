#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "stdlite/errors.h"
#include "stdlite/utils/fchk_parser.h"


int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar) {
    assert(lx != NULL && name != NULL && type != NULL && is_scalar != NULL);

    if(lx->current_tk_type != STDL_TK_ALPHA) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected FCHK section to start with ALPHA");
        return STDL_ERR_UTIL_FCHK;
    }

    // read up the 40 next characters (name)
    int err = STDL_ERR_OK;
    char buff[41]; // 40 + '\0'
    int i = 0, last_alnum = 0;

    while (i < 40) {
        if(lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF)  {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "unexpected end in section name");
            err = STDL_ERR_UTIL_FCHK;
        }

        if(err == STDL_ERR_OK) {
            buff[i] = lx->current_tk_value;

            if(isalnum(lx->current_tk_value))
                last_alnum = i;

            err = stdl_lexer_advance(lx, 1);
        }

        RETURN_ON_ERROR(err);
        i++;
    }

    // read out the next 7 character (only the fourth is relevant)
    i = 0;
    while (i < 7) {
        if (i == 3) {
            if(lx->current_tk_type != STDL_TK_ALPHA || (lx->current_tk_value != 'I' && lx->current_tk_value != 'R' && lx->current_tk_value != 'C')) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expecting data type to be I/R/C");
                err = STDL_ERR_UTIL_FCHK;
            } else
                *type = lx->current_tk_value;
        }

        if(err == STDL_ERR_OK)
            err = stdl_lexer_advance(lx, 1);

        RETURN_ON_ERROR(err);
        i++;
    }

    // time to know if it is a scalar or a vector
    if (lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "line is too short to check for kind (scalar/vector)");
        err = STDL_ERR_UTIL_FCHK;
    } else {
        *is_scalar = lx->current_tk_type != STDL_TK_ALPHA || lx->current_tk_value != 'N';

        // just advance to the value/size
        if(!(*is_scalar)) { // skip '='
            err = stdl_lexer_eat(lx, STDL_TK_ALPHA);

            if(err == STDL_ERR_OK)
                err = stdl_lexer_eat(lx, STDL_TK_EQ);

            if(err != STDL_ERR_OK) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected EQ after 'N'");
                err = STDL_ERR_UTIL_FCHK;
            }
        }
    }

    if(err == STDL_ERR_OK)
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

int _skip_NL(stdl_lexer* lx) {
    /* Skips a NL. Stops right after.
     */
    int err = STDL_ERR_OK;

    if(lx->current_tk_type != STDL_TK_EOF) {
        err = stdl_lexer_eat(lx, STDL_TK_NL);
        if(err != STDL_ERR_OK) {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected NL");
            err = STDL_ERR_UTIL_FCHK;
        }
    }

    return err;
}


int stdl_fchk_parser_get_scalar_int(stdl_lexer* lx, long *value) {
    int err = stdl_parser_get_integer(lx, value);
    if(err == STDL_ERR_OK)
        err = _skip_NL(lx);

    return err;
}


int stdl_fchk_parser_get_scalar_number(stdl_lexer* lx, double* value) {
    int err = stdl_parser_get_number(lx, value);
    if(err == STDL_ERR_OK)
        err = _skip_NL(lx);

    return err;
}

int _get_vec_sz(stdl_lexer* lx, size_t* sz) {
    /* Reads the size. Stops after NL.
     */
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

    err = _skip_NL(lx);

    if(err == STDL_ERR_OK)
        *sz = (size_t) sz_read;

    return err;
}


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
        err = stdl_lexer_skip(lx, isspace);

        if(err == STDL_ERR_OK)
            err = stdl_parser_get_integer(lx, &x);

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 6 == 5 && i != *sz - 1)
                err = _skip_NL(lx);
        }

        if(err != STDL_ERR_OK) {
            free(*vector);
            return err;
        }

        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);

    if(err != STDL_ERR_OK) {
        free(*vector);
        return err;
    }

    // ok, we should be at the next section
    return STDL_ERR_OK;
}


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
        err = stdl_lexer_skip(lx, isspace);
        if(err == STDL_ERR_OK) {
            err = stdl_parser_get_number(lx, &x);
        }

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 5 == 4 && i != *sz - 1)
                err = _skip_NL(lx);
        }

        if(err != STDL_ERR_OK) {
            free(*vector);
            return err;
        }

        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);

    if(err != STDL_ERR_OK) {
        free(*vector);
        return err;
    }

    // now, we should be at the next section
    return STDL_ERR_OK;
}

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
        // read the 12 characters of a pack
        j = 0;
        while (j < 12 && err == STDL_ERR_OK) {
            (*out)[i * 12 + j] = lx->current_tk_value;
            err = stdl_lexer_advance(lx, 1);
            j++;
        }

        if(err == STDL_ERR_OK && i % 5 == 4 && i != *sz - 1)
            err = _skip_NL(lx);

        if(err != STDL_ERR_OK) {
            free(*out);
            return err;
        }

        i++;
    }

    // skip last NL if any
    err = _skip_NL(lx);

    if(err != STDL_ERR_OK) {
        free(*out);
        return err;
    }

    (*out)[(*sz) * 12] = '\0';

    return STDL_ERR_OK;
}

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
