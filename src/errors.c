#include <stdio.h>
#include <assert.h>
#include <stdarg.h>

#include "stdlite/errors.h"

int DEBUG_LVL = 1; // by default, warnings

/**
 * Set `DEBUG_LVL`.
 * @param level any number, the larger the more messages one get.
 * @ingroup errors
 */
void stdl_set_debug_level(const int level) {
    DEBUG_LVL = level;
}

/**
 * Print a debug message in `stdout`, if `DEBUG_LVL` is above 1.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_debug_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    if(DEBUG_LVL <= 1)
        return;

    va_list arglist;

    printf("DEBUG (%s:%d): ", file, line);
    va_start(arglist, format);
    vprintf(format, arglist);
    va_end(arglist);
    printf("\n");
}

/**
 * Print a warning message in `stdout`, if `DEBUG_LVL` is above 0.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_warning_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    if(DEBUG_LVL < 1)
        return;

    va_list arglist;

    printf("WARN (%s:%d): ", file, line);
    va_start(arglist, format);
    vprintf(format, arglist);
    va_end(arglist);
    printf("\n");
}

/**
 * Print an error message in `stderr`.
 * @param file source file (use `__FILE__`)
 * @param line line (use `__LINE__`)
 * @param format format of the string
 * @param ... extra parameters
 * @ingroup errors
 */
void stdl_error_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    va_list arglist;

    fprintf(stderr, "ERROR (%s:%d): ", file, line);
    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, "\n");
}
