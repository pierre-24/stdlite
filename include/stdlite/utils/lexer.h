#ifndef STDL_LEXER_H
#define STDL_LEXER_H

#include <stdio.h>

#define STDL_LEXER_STREAM_BUFF_SIZE 1024

/**
 * Token types
 * @ingroup lexer
 */
typedef enum stdl_token_type_ {
    STDL_TK_WHITESPACE, // space (U+0020), horizontal tab (U+0009)
    STDL_TK_NL, // linefeed (\n, U+000A)
    STDL_TK_CR, // carriage return (\r, U+000D)
    STDL_TK_DIGIT, // [0-9]
    STDL_TK_ALPHA, // [a-zA-Z]
    STDL_TK_COMMA, // ","
    STDL_TK_DOT, // "."
    STDL_TK_ESCAPE, // "\"
    STDL_TK_QUOTE, // "
    STDL_TK_DASH, // "-"
    STDL_TK_PLUS, // "+"
    STDL_TK_EQ, // =

    STDL_TK_EOF, // → end of string

    STDL_TK_CHAR, // anything but what is before

    STDL_TK_LAST
} stdl_token_type;

/**
 * Lexer object
 * @ingroup lexer
 */
typedef struct stdl_lexer_ {
    // streaming
    FILE* file;
    char* stream;
    int pos_in_stream;

    // position in file
    int current_line;
    int current_pos_in_line;

    // current token
    char current_tk_value;
    stdl_token_type current_tk_type;
} stdl_lexer;


stdl_lexer* stdl_lexer_new(FILE* input);
int stdl_lexer_delete(stdl_lexer *lx);

int stdl_lexer_advance(stdl_lexer *lx, int shift);
int stdl_lexer_eat(stdl_lexer *lx, stdl_token_type t);
int stdl_lexer_skip(stdl_lexer *lx, stdl_token_type t);
int stdl_lexer_skip_whitespace_and_nl(stdl_lexer *lx);

#endif //STDL_LEXER_H
