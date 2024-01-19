#ifndef STDL_FCHK_PARSER_H
#define STDL_FCHK_PARSER_H

#include <stdlite/utils/base_parser.h>

int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar);

int stdl_fchk_parser_get_vector_ints(stdl_lexer* lx, size_t* sz, long **vector);
int stdl_fchk_parser_get_vector_numbers(stdl_lexer* lx, size_t* sz, double** vector);
int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, char **out);

int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar);
int stdl_fchk_parser_skip_begin(stdl_lexer* lx);

#endif //STDL_FCHK_PARSER_H
