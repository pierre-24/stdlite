#ifndef STDL_ERRORS_H
#define STDL_ERRORS_H

enum {
    STDL_ERR_OK,

    // unexpected
    STDL_ERR_MALLOC,
    STDL_ERR_READ,
    STDL_ERR_NOT_FOUND,
    STDL_ERR_ALREADY,

    // specific errors
    STDL_ERR_UTIL_LEXER,
    STDL_ERR_UTIL_PARSER,
    STDL_ERR_UTIL_FCHK,

    STDL_ERR_LAST
};

void stdl_set_debug_level(const int level);
void stdl_debug_msg(char *file, int line, char *format, ...);
void stdl_warning_msg(char *file, int line, char *format, ...);
void stdl_error_msg(char *file, int line, char *format, ...) ;


#define RETURN_ON_ERROR(a) if((a) != STDL_ERR_OK) \
return (a);


#endif //STDL_ERRORS_H
