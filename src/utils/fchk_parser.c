#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "stdlite/errors.h"
#include "stdlite/utils/fchk_parser.h"

/**
 * Parse the info of a section in FCHK.
 * @param lx a valid lexer
 * @param[out] name the name of the section
 * @param[out] type the type of section. Valid outputs are `I`, `R`, and `C`.
 * @param[out] is_scalar `1` if the section is a scalar
 * @return `STDL_ERR_OK`
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
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expecting data type");
                return STDL_ERR_UTIL_FCHK;
            }

            *type = lx->current_tk_value;
        }

        err = stdl_lexer_advance(lx, 1);
        RETURN_ON_ERROR(err);

        i++;
    }

    // time to know if it is a scalar or a vector
    *is_scalar = lx->current_tk_type != STDL_TK_ALPHA || lx->current_tk_value != 'N';

    // just advance to the value/size
    if(!(*is_scalar)) { // skip '='
        err = stdl_lexer_advance(lx, 1);
        RETURN_ON_ERROR(err);

        if(lx->current_tk_type != STDL_TK_EQ) {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected EQ after 'N'");
            return STDL_ERR_UTIL_FCHK;
        }

        err = stdl_lexer_advance(lx, 1);
        RETURN_ON_ERROR(err);
    }

    err = stdl_lexer_skip_whitespace_and_nl(lx);
    RETURN_ON_ERROR(err);

    // ok, normally we are good. Time to copy the name
    *name = malloc((last_alnum + 2) * sizeof(char));
    if(*name == NULL)
        return STDL_ERR_MALLOC;

    memcpy(*name, buff, last_alnum + 1);
    (*name)[last_alnum + 1] = '\0';

    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_int(stdl_lexer* lx, size_t* sz, int** vector) {
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_real(stdl_lexer* lx, size_t* sz, double** vector) {
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, double** vector) {
    return STDL_ERR_OK;
}

int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar) {
    return STDL_ERR_OK;
}