#ifndef STDL_ERRORS_H
#define STDL_ERRORS_H

/**
 * Enum for the errors codes.
 * @ingroup errors
 */
enum stdl_error_code_ {
    /// Everything went well
    STDL_ERR_OK = 0,

    /// `malloc()` failure
    STDL_ERR_MALLOC = 1,

    /// Read error
    STDL_ERR_READ = 2,

    /// Element not found
    STDL_ERR_NOT_FOUND = 3,

    /// Element already found
    STDL_ERR_ALREADY = 4,

    STDL_UNUSED5 = 5,
    STDL_UNUSED6 = 6,
    STDL_UNUSED7 = 7,
    STDL_UNUSED8 = 8,
    STDL_UNUSED9 = 9,
    STDL_UNUSED10 = 10,

    /// Error in the `lexer` module
    STDL_ERR_UTIL_LEXER = 11,

    /// Error in the `base_parser` module
    STDL_ERR_UTIL_PARSER = 12,

    /// Error in the `fchk_parser` module
    STDL_ERR_UTIL_FCHK = 13,

    STDL_ERR_LAST = 14
};

/**
 * Set `DEBUG_LVL`.
 * @param level any number, the larger the more messages one get. Setting this number to any negative totally shut down any message.
 * @ingroup errors
 */
void stdl_set_debug_level(const int level);

/**
 * Get the value of `DEBUG_LVL`.
 * @return the debug level
 * @ingroup errors
 */
int stdl_get_debug_level();

/**
 * Print a debug message in `stdout`, if `DEBUG_LVL` is above 1.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_debug_msg(char *file, int line, char *format, ...);

/**
 * Print a warning message in `stdout`, if `DEBUG_LVL` is above 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_warning_msg(char *file, int line, char *format, ...);

/**
 * Print an error message in `stderr`, if `DEBUG_LVL` is equal or above 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_error_msg(char *file, int line, char *format, ...) ;


#define RETURN_ON_ERROR(a) if((a) != STDL_ERR_OK) \
return (a);


#endif //STDL_ERRORS_H
