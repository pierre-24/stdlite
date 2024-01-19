#ifndef STDL_BASE_PARSER_H
#define STDL_BASE_PARSER_H

#include <stdlite/utils/lexer.h>

/**
 * Multiplier for string growth.
 * @ingroup base_parser
 * */
#define STDL_STR_MULT 128

int stdl_grow_string(char** str, int sz, int* fac);

void stdl_error_msg_parser(char *file, int line, stdl_lexer* lx, char *format, ...);

int stdl_parser_get_integer(stdl_lexer* lx, long *result);
int stdl_parser_get_number(stdl_lexer* lx, double* result);
int stdl_parser_get_literal(stdl_lexer* lx, int (*predicate)(int), char** result);

#endif //STDL_BASE_PARSER_H
