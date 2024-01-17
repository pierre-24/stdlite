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

    STDL_ERR_LAST
};

#endif //STDL_ERRORS_H
