#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>

#include "stdlite/helpers.h"
#include "stdlite/logging.h"
#include "stdlite/utils/lexer.h"

int lexer_translator[] = {
    0x20, STDL_TK_WHITESPACE,
    0x0a, STDL_TK_NL,
    0x0d, STDL_TK_CR,
    0x09, STDL_TK_WHITESPACE,
    '0', STDL_TK_DIGIT,
    '1', STDL_TK_DIGIT,
    '2', STDL_TK_DIGIT,
    '3', STDL_TK_DIGIT,
    '4', STDL_TK_DIGIT,
    '5', STDL_TK_DIGIT,
    '6', STDL_TK_DIGIT,
    '7', STDL_TK_DIGIT,
    '8', STDL_TK_DIGIT,
    '9', STDL_TK_DIGIT,
    '=', STDL_TK_EQ,
    ',', STDL_TK_COMMA,
    '.', STDL_TK_DOT,
    '\\', STDL_TK_ESCAPE,
    '"', STDL_TK_QUOTE,
    '-', STDL_TK_DASH,
    '+', STDL_TK_PLUS,
    '\0', STDL_TK_EOF,

    -1
};


void stdl_error_msg_lexer(char *file, int line, stdl_lexer* lx, char *format, ...) {
    assert(file != NULL && lx != NULL && format != NULL);

    if(stdl_get_debug_level() < 0)
        return;

    va_list arglist;

    char buff[64];

    fprintf(stderr, "ERROR (%s:%d):", file, line);

    sprintf(buff, isgraph(lx->current_tk_value) ? "%c": "0x%x", lx->current_tk_value);
    fprintf(stderr, "file@%d+%d(%s=%d): ", lx->current_line, lx->current_pos_in_line, buff, lx->current_tk_type);

    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, "\n");
}


int stdl_lexer_new(FILE *input, stdl_lexer **lx_ptr) {
    assert(lx_ptr != NULL && input != NULL);

    // go ahead
    *lx_ptr = malloc(sizeof(stdl_lexer));
    STDL_ERROR_HANDLE_AND_REPORT(*lx_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create lexer %p", *lx_ptr);

    (*lx_ptr)->file = input;

    (*lx_ptr)->stream = malloc(STDL_LEXER_STREAM_BUFF_SIZE * sizeof(char));
    STDL_ERROR_HANDLE_AND_REPORT((*lx_ptr)->stream == NULL, stdl_lexer_delete(*lx_ptr); return STDL_ERR_MALLOC, "malloc");

    // read first bytes
    size_t n = fread((*lx_ptr)->stream, 1, STDL_LEXER_STREAM_BUFF_SIZE, input);
    if (n != STDL_LEXER_STREAM_BUFF_SIZE) // shorter than expected!
        (*lx_ptr)->stream[n] = '\0';

    (*lx_ptr)->pos_in_stream = 0;
    (*lx_ptr)->current_line = 1;
    (*lx_ptr)->current_pos_in_line = 1;
    (*lx_ptr)->current_tk_value = (*lx_ptr)->stream[0];

    stdl_lexer_advance(*lx_ptr, 0);

    return STDL_ERR_OK;
}


int stdl_lexer_delete(stdl_lexer *lx) {
    assert(lx != NULL);

    STDL_DEBUG("delete lexer %p", lx);

    STDL_FREE_ALL(lx->stream, lx);

    return STDL_ERR_OK;
}


int stdl_lexer_advance(stdl_lexer *lx, int shift) {
    assert(lx != NULL && lx->file != NULL && lx->stream != NULL);
    assert(shift == 0 || shift == 1);

    STDL_LEXER_ERROR_HAR(lx, lx->current_tk_value == '\0', lx->current_tk_type = STDL_TK_EOF; return STDL_ERR_UTIL_LEXER, "reading past EOF");

    if (lx->pos_in_stream == STDL_LEXER_STREAM_BUFF_SIZE - 1) {
        // need to read next bytes
        lx->stream[0] = lx->stream[STDL_LEXER_STREAM_BUFF_SIZE - 1]; // wrap around
        lx->pos_in_stream = 0;

        size_t n = fread(lx->stream + 1, 1, STDL_LEXER_STREAM_BUFF_SIZE - 1, lx->file);
        if (n != STDL_LEXER_STREAM_BUFF_SIZE - 1) // shorter than expected!
            lx->stream[n + 1] = '\0';
    }

    char c = lx->stream[lx->pos_in_stream + shift];
    stdl_token_type t = STDL_TK_CHAR;

    if (isalpha(c))
        t = STDL_TK_ALPHA;
    else {
        int *tr = lexer_translator;
        while (*tr != -1) {
            if (c == *tr) {
                t = *(tr + 1);
                break;
            }

            tr += 2;
        }
    }

    lx->pos_in_stream += shift;
    lx->current_tk_type = t;
    lx->current_tk_value = c;

    if (lx->current_tk_type == STDL_TK_NL) {
        lx->current_line += 1;
        lx->current_pos_in_line = 0;
    } else {
        lx->current_pos_in_line += shift;
    }

    return STDL_ERR_OK;
}


int stdl_lexer_eat(stdl_lexer *lx, stdl_token_type t) {
    assert(lx != NULL);
    STDL_LEXER_ERROR_HAR(lx, lx->current_tk_type != t, return STDL_ERR_UTIL_LEXER, "expected token type %d, got %d", t, lx->current_tk_type);

    return stdl_lexer_advance(lx, 1);
}


int stdl_lexer_skip(stdl_lexer *lx, int (*predicate)(int)) {
    assert(lx != NULL);

    int err = STDL_ERR_OK;
    while (lx->current_tk_type != STDL_TK_EOF && predicate(lx->current_tk_value) && err == STDL_ERR_OK)
        err = stdl_lexer_advance(lx, 1);

    return STDL_ERR_OK;
}


int stdl_lexer_skip_whitespace_and_nl(stdl_lexer *lx) {
    assert(lx != NULL);

    int err = STDL_ERR_OK;
    while ((lx->current_tk_type == STDL_TK_WHITESPACE || lx->current_tk_type == STDL_TK_NL) && err == STDL_ERR_OK) {
        err = stdl_lexer_advance(lx, 1);
    }

    return err;
}
