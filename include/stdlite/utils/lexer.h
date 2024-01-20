#ifndef STDL_LEXER_H
#define STDL_LEXER_H

#include <stdio.h>

/**
 * Buffer size for the lexer.
 * @ingroup lexer
 * */
#define STDL_LEXER_STREAM_BUFF_SIZE 1024

/**
 * Enum for the different token types.
 * To be used in the context of parsing files.
 * @ingroup lexer
 */
enum stdl_token_type_ {
    /// Space (`0x20`), horizontal tab (`0x9`)
    STDL_TK_WHITESPACE,
    /// Linefeed (`\n`, `U+000A`)
    STDL_TK_NL,
    /// Carriage return (`\r`, `U+000D`)
    STDL_TK_CR,
    /// Digit, `[0-9]`
    STDL_TK_DIGIT,
    /// Alpha character, `[a-zA-Z]`
    STDL_TK_ALPHA,
    /// Comma, `,`
    STDL_TK_COMMA,
    /// Dot, `.`
    STDL_TK_DOT,
    /// Escape character, `\`
    STDL_TK_ESCAPE,
    /// Quote, `"`
    STDL_TK_QUOTE,
    /// Dash, `-`
    STDL_TK_DASH,
    /// Plus, `+`
    STDL_TK_PLUS,
    /// Equal, `=`
    STDL_TK_EQ,
    /// End of file, `0x0`
    STDL_TK_EOF,
    /// Anything else.
    STDL_TK_CHAR,

    STDL_TK_LAST
};

typedef enum stdl_token_type_ stdl_token_type;


/**
 * A structure that represent a lexer.
 * This object contains a buffer, `stream`, that is automatically filled by `stdl_lexer_advance` when it gets close to the end.
 * @ingroup lexer
 */
struct stdl_lexer_ {
    /** Current position of the token in the temporary buffer. */
    int pos_in_stream;

    /** Current position of the token in the line. */
    int current_pos_in_line;

    /** Current line in the file. */
    int current_line;

    /** Current token value, equals to `lx->stream[lx->pos_in_stream]`. */
    char current_tk_value;

    /** Current token type. */
    stdl_token_type current_tk_type;

    /** The file from which character are read. */
    FILE* file;

    /** Temporary buffer. */
    char* stream;
};

typedef struct stdl_lexer_ stdl_lexer;

/**
 * Create a new lexer object (`stdl_lexer*`), provided that `input` is a valid `FILE*`.
 * This object has to be free'd after use.
 * @param input A valid `FILE*`.
 * @return a `stdl_lexer` object.
 * @ingroup lexer
 */
stdl_lexer* stdl_lexer_new(FILE* input);

/**
 * Free the lexer.
 * @param lx a valid lexer
 * @return `STDL_ERR_OK`
 * @ingroup lexer
 */
int stdl_lexer_delete(stdl_lexer *lx);

/**
 * Advance to the next token, except if there are no more token to extract (current token is EOF).
 * @param lx a valid lexer
 * @param shift How much to advance. Valid choices are 0 or 1.
 * @return `STDL_ERR_OK` if current token is not `STDL_TK_EOF`, `STDL_ERR_UTIL_LEXER` otherwise.
 * @ingroup lexer
 */
int stdl_lexer_advance(stdl_lexer *lx, int shift);

/**
 * If the current token type is of `type`, advance to the next token. Otherwise, return an error.
 * @param lx A valid lexer.
 * @param t Type of the token
 * @return `STDL_ERR_OK` if the current token was of the correct type, `STDL_ERR_UTIL_LEXER` otherwise.
 * @ingroup lexer
 */
int stdl_lexer_eat(stdl_lexer *lx, stdl_token_type t);

/**
 * Skip tokens while `predicate` is true (and one does not hit `EOF`).
 * @param lx A valid lexer.
 * @param predicate predicate of the form `int predicate(int c)`, where `c` is the current token value
 * @return `STDL_ERR_OK`.
 * @ingroup lexer
 */
int stdl_lexer_skip(stdl_lexer *lx, int (*predicate)(int));

/**
 * Skip all tokens that are `STDL_TK_WHITESPACE` or `STDL_TK_NL`.
 * @param lx A valid lexer.
 * @return `STDL_ERR_OK`.
 * @ingroup lexer
 */
int stdl_lexer_skip_whitespace_and_nl(stdl_lexer *lx);

#endif //STDL_LEXER_H
