#include <stdio.h>
#include <assert.h>
#include <stdarg.h>

#include "stdlite/errors.h"

int DEBUG_LVL = 1; // by default, warnings & errors

void stdl_set_debug_level(const int level) {
    DEBUG_LVL = level;
}

int stdl_get_debug_level() {
    return DEBUG_LVL;
}

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

void stdl_error_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    if(DEBUG_LVL < 0)
        return;

    va_list arglist;

    fprintf(stderr, "ERROR (%s:%d): ", file, line);
    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
    va_end(arglist);
    fprintf(stderr, "\n");
}
