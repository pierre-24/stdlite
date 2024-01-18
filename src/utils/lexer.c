#include <assert.h>
#include <ctype.h>
#include <stdlib.h>

#include "stdlite/errors.h"
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

/**
 * Create a new lexer object (`stdl_lexer*`), provided that `input` is a valid `FILE*`.
 * This object has to be free'd after use.
 * @param input A valid `FILE*`.
 * @return a `stdl_lexer` object.
 * @ingroup lexer
 */
stdl_lexer *stdl_lexer_new(FILE *input) {
    assert(input != NULL);

    // go ahead
    stdl_lexer *lx = malloc(sizeof(stdl_lexer));
    if (lx != NULL) {
        lx->stream = NULL;
        lx->file = input;

        lx->stream = malloc(STDL_LEXER_STREAM_BUFF_SIZE * sizeof(char));
        if (lx->stream == NULL) {
            stdl_lexer_delete(lx);
            lx = NULL;
        } else {
            // read first bytes
            size_t n = fread(lx->stream, 1, STDL_LEXER_STREAM_BUFF_SIZE, input);
            if (n != STDL_LEXER_STREAM_BUFF_SIZE) // shorter than expected!
                lx->stream[n] = '\0';

            lx->pos_in_stream = 0;
            lx->current_line = 1;
            lx->current_pos_in_line = 0;
            lx->current_tk_value = lx->stream[0];

            stdl_lexer_advance(lx, 0);
        }
    }

    return lx;
}

/**
 * Free the lexer.
 * @param lx a valid lexer
 * @return `STDL_ERR_OK`
 * @ingroup lexer
 */
int stdl_lexer_delete(stdl_lexer *lx) {
    assert(lx != NULL);

    if (lx->stream != NULL)
        free(lx->stream);

    free(lx);

    return STDL_ERR_OK;
}

/**
 * Advance to the next token, except if there are no more token to extract (current token is EOF).
 * @param lx a valid lexer
 * @param shift How much to advance. Valid choices are 0 or 1.
 * @return `STDL_ERR_OK` if current token is not `STDL_TK_EOF`, `STDL_ERR_UTIL_LEXER` otherwise.
 * @ingroup lexer
 */
int stdl_lexer_advance(stdl_lexer *lx, int shift) {
    assert(lx != NULL && lx->file != NULL && lx->stream != NULL);
    assert(shift == 0 || shift == 1);

    if(lx->current_tk_value == '\0') {
        lx->current_tk_type = STDL_TK_EOF;
        return STDL_ERR_UTIL_LEXER;
    }

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

/**
 * If the current token type is of `type`, advance to the next token. Otherwise, return an error.
 * @param lx A valid lexer.
 * @param t Type of the token
 * @return `STDL_ERR_OK` if the current token was of the correct type, `STDL_ERR_UTIL_LEXER` otherwise.
 * @ingroup lexer
 */
int stdl_lexer_eat(stdl_lexer *lx, stdl_token_type t) {
    assert(lx != NULL);

    if (lx->current_tk_type != t)
        return STDL_ERR_UTIL_LEXER;

    return stdl_lexer_advance(lx, 1);
}

/**
 * If the current token type is of `type`, advance to the next token. Otherwise, do nothing.
 * @param lx A valid lexer.
 * @param t Type of the token
 * @return `STDL_ERR_OK`.
 * @ingroup lexer
 */
int stdl_lexer_skip(stdl_lexer *lx, stdl_token_type t) {
    assert(lx != NULL);

    while (lx->current_tk_type == t)
        stdl_lexer_advance(lx, 1);

    return STDL_ERR_OK;
}

/**
 * Skip all tokens that are `STDL_TK_WHITESPACE` or `STDL_TK_NL`.
 * @param lx A valid lexer.
 * @return `STDL_ERR_OK`.
 * @ingroup lexer
 */
int stdl_lexer_skip_whitespace_and_nl(stdl_lexer *lx) {
    assert(lx != NULL);

    while (lx->current_tk_type == STDL_TK_WHITESPACE || lx->current_tk_type == STDL_TK_NL) {
        stdl_lexer_advance(lx, 1);
    }

    return STDL_ERR_OK;
}
