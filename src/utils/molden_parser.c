#include <assert.h>

#include "stdlite/utils/molden_parser.h"
#include "stdlite/logging.h"

int stdl_molden_parser_read_section_title(stdl_lexer* lx, char** title) {

    assert(lx != NULL && title != NULL);

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != '[',
            return STDL_ERR_UTIL_MOLDEN,
            "expected MOLDEN section to start with `[`"
    );

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_pos_in_line != 1,
            return STDL_ERR_UTIL_MOLDEN,
            "expected MOLDEN section to start with `[` as the first character of a line"
    );

    stdl_lexer_advance(lx, 1);

    // read until `]`
    int err;
    char* str;
    int sz = 0, fac = 0;

    err = stdl_grow_string(&str, sz, &fac);
    STDL_ERROR_CODE_HANDLE(err, return err);

    while (lx->current_tk_type != STDL_TK_EOF  && lx->current_tk_type != STDL_TK_NL && lx->current_tk_value != ']') {
        err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
        STDL_ERROR_CODE_HANDLE(err, free(str); return err);
    }

    STDL_LEXER_ERROR_HAR(lx, lx->current_tk_value != ']', free(str); return STDL_ERR_UTIL_MOLDEN, "expected MOLDEN section to end with `]`");

    str[sz] = '\0';
    *title = str;

    return STDL_ERR_OK;
}

int _pred_not_begining_of_section(int c) {
    return c != '[';
}

int stdl_molden_parser_skip_section(stdl_lexer* lx) {
    assert(lx != NULL);

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to start section"
    );

    return stdl_lexer_skip(lx, _pred_not_begining_of_section);
}
