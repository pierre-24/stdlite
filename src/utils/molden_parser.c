#include <assert.h>
#include <ctype.h>
#include <string.h>

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
            "expected `]` to continue section"
    );

    return stdl_lexer_skip(lx, _pred_not_begining_of_section);
}

// just a linked-chain to temporary store atoms
struct _atom_info_ {
    double atm[4];
    struct _atom_info_* next;
};

int _pred_notnl(int c) {
    return c != '\n';
}

int stdl_molden_parser_read_atoms_section(stdl_lexer* lx, size_t* natm, double** atm) {
    assert(lx != NULL && natm != NULL && atm != NULL);

    STDL_DEBUG("Read [Atoms] section");

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to continue section"
    );

    int err;

    err = stdl_lexer_eat(lx, STDL_TK_CHAR);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // Get unit
    int is_bohr = 0;
    err = stdl_lexer_skip(lx, isspace);
    STDL_ERROR_CODE_HANDLE(err, return err);

    char* unit = NULL;
    err = stdl_parser_get_literal(lx, isalpha, &unit);
    STDL_ERROR_CODE_HANDLE(err, return err);

    if(strcmp(unit, "AU") == 0)
        is_bohr = 1;

    free(unit);

    // go to next line
    err = stdl_lexer_skip_whitespace_and_nl(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // read atoms
    *natm = 0;
    struct _atom_info_* first = NULL;
    struct _atom_info_* current = NULL;
    while (lx->current_tk_type != STDL_TK_EOF && lx->current_tk_value != '[') {

        char* tmp;
        struct _atom_info_* inf = malloc(sizeof(struct _atom_info_));
        inf->next = NULL;

        err = stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, return err);
        free(tmp);

        err = stdl_lexer_skip(lx, isspace);
        STDL_ERROR_CODE_HANDLE(err, return err);

        double n;
        err = stdl_parser_get_number(lx, &n);
        STDL_ERROR_CODE_HANDLE(err, return err);
        STDL_LEXER_ERROR_HAR(lx, n != *natm + 1, return STDL_ERR_UTIL_MOLDEN, "expected atoms that are in consecutive order");

        err = stdl_lexer_skip(lx, isspace);
        STDL_ERROR_CODE_HANDLE(err, return err);

        err = stdl_parser_get_number(lx, &n);
        STDL_ERROR_CODE_HANDLE(err, return err);
        inf->atm[0] = (double ) n;

        double d;
        for (int i = 0; i < 3; ++i) {
            err = stdl_lexer_skip(lx, isspace);
            STDL_ERROR_CODE_HANDLE(err, return err);

            err = stdl_parser_get_number(lx, &d);
            STDL_ERROR_CODE_HANDLE(err, return err);
            inf->atm[i + 1] = d;
        }

        err = stdl_lexer_skip_whitespace_and_nl(lx);
        STDL_ERROR_CODE_HANDLE(err, return err);

        if(first == NULL) {
            first = inf;
            current = first;
        } else {
            current->next = inf;
            current = current->next;
        }

        *natm += 1;
    }

    STDL_DEBUG("Got %ld atoms", *natm);

    // actually create the `atm` array, now that we know its size
    *atm = malloc(*natm * 4 * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(*atm == NULL, return STDL_ERR_MALLOC, "malloc");

    size_t iatm = 0;
    current = first;
    struct _atom_info_* prev;

    while (current != NULL) {
        memcpy(*atm + iatm * 4, current->atm, 4 * sizeof (double ));
        prev = current;
        current = current->next;
        free(prev); // free linked list while we're at it

        iatm += 1;
    }

    return STDL_ERR_OK;
}
