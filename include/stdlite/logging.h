#ifndef STDL_LOGGING_H
#define STDL_LOGGING_H

/// Handle error (if `assertion` is true) by printing a message and using an action.
/// Inspired by B. Klemens in *21st century C* (O'Reilly).
/// @ingroup logging
#define STDL_ERROR_HANDLE_AND_REPORT(assertion, error_action, ...)   \
    if(assertion) {                                                  \
        stdl_error_msg(__FILE__, __LINE__, __VA_ARGS__);             \
        {error_action;}                                              \
    }

/// Handle error (if `assertion` is true) by using an action.
/// @ingroup logging
#define STDL_ERROR_HANDLE(assertion, error_action) if(assertion){error_action;}

/// Handle error code (if different from `STDL_ERR_OK`) by using an action.
/// @ingroup logging
#define STDL_ERROR_CODE_HANDLE(code, error_action) STDL_ERROR_HANDLE(code != STDL_ERR_OK, error_action)


/**
 * Enum for the errors codes.
 * @ingroup logging
 */
enum stdl_error_code_ {
    /// Everything went well
    STDL_ERR_OK = 0,

    /// `malloc()` failure
    STDL_ERR_MALLOC = 1,

    /// Open error
    STDL_ERR_OPEN = 2,

    /// Read error
    STDL_ERR_READ = 3,

    /// Write error
    STDL_ERR_WRITE = 4,

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

    /// Error in the `molden_parser` module
    STDL_ERR_UTIL_MOLDEN = 14,

    /// Error in the `basis` module
    STDL_ERR_BASIS = 15,

    /// Error in the `response` module
    STDL_ERR_RESPONSE = 16,

    STDL_ERR_LAST = 17
};

/**
 * Set `LOG_LVL`.
 * @param level any number, the larger the more messages one get. Setting this number to any negative totally shut down any message.
 * @ingroup logging
 */
void stdl_set_log_level(const int level);

/**
 * Get the value of `LOG_LVL` (default is 1).
 * @return the log level
 * @ingroup logging
 */
int stdl_get_log_level();

/**
 * Print a log message in `stdout` if `loglevel >= LOG_LVL`.
 * @param loglevel logging level of the messafe
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup logging
 */
void stdl_log_msg(int loglevel, char *format, ...);

/**
 * Set `DEBUG_LVL`.
 * @param level any number, the larger the more messages one get. Setting this number to any negative totally shut down any message.
 * @ingroup logging
 */
void stdl_set_debug_level(const int level);

/**
 * Get the value of `DEBUG_LVL` (default is 1).
 * @return the debug level
 * @ingroup logging
 */
int stdl_get_debug_level();

/**
 * Print a debug message in `stdout`, if `DEBUG_LVL` is above 1.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup logging
 */
void stdl_debug_msg(char *file, int line, char *format, ...);

/**
 * Print a debug message
 * @ingroup logging
 */
#define STDL_DEBUG(...) stdl_debug_msg(__FILE__, __LINE__, __VA_ARGS__)

/**
 * Print a warning message in `stdout`, if `DEBUG_LVL` is above 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup logging
 */
void stdl_warning_msg(char *file, int line, char *format, ...);

/**
 * Print a warning message
 * @ingroup logging
 */
#define STDL_WARN(...) stdl_warning_msg(__FILE__, __LINE__, __VA_ARGS__)

/**
 * Print an error message in `stderr`, if `DEBUG_LVL` is equal or above 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup logging
 */
void stdl_error_msg(char *file, int line, char *format, ...) ;


#endif //STDL_LOGGING_H
