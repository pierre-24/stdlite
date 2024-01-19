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

/**
 * `typedef` for `stdl_token_type`.
 * @ingroup lexer
 */
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

/**
 * `typedef` for `stdl_lexer`.
 * @ingroup lexer
 */
typedef struct stdl_lexer_ stdl_lexer;

stdl_lexer* stdl_lexer_new(FILE* input);
int stdl_lexer_delete(stdl_lexer *lx);

int stdl_lexer_advance(stdl_lexer *lx, int shift);
int stdl_lexer_eat(stdl_lexer *lx, stdl_token_type t);
int stdl_lexer_skip(stdl_lexer *lx, int (*predicate)(int));
int stdl_lexer_skip_whitespace_and_nl(stdl_lexer *lx);

#endif //STDL_LEXER_H
