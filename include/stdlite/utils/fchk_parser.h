#ifndef STDL_FCHK_PARSER_H
#define STDL_FCHK_PARSER_H

#include <stdlite/utils/base_parser.h>

int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar);

int stdl_fchk_parser_get_vector_int(stdl_lexer* lx, size_t* sz, int** vector);
int stdl_fchk_parser_get_vector_real(stdl_lexer* lx, size_t* sz, double** vector);
int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, double** vector);

int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar);

#endif //STDL_FCHK_PARSER_H
