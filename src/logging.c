#include <stdio.h>
#include <assert.h>
#include <stdarg.h>

#include "stdlite/logging.h"

char* stdl_library_name() {
    return PROJECT_NAME;
}

char* stdl_library_version()  {
    return PROJECT_VERSION;
}

int LOG_LVL = 0; // by default, moderate logging

void stdl_set_log_level(const int level) {
    LOG_LVL = level;
}

int stdl_get_log_level() {
    return LOG_LVL;
}

int DEBUG_LVL = 1; // by default, warnings & errors but not debug

void stdl_set_debug_level(const int level) {
    DEBUG_LVL = level;
}

int stdl_get_debug_level() {
    return DEBUG_LVL;
}

void stdl_log_msg(int loglevel, char *format, ...) {
    assert(format != NULL);

    if(LOG_LVL < loglevel)
        return;

    va_list arglist;

    va_start(arglist, format);
    vfprintf(stdout, format, arglist);
    va_end(arglist);
}

void stdl_debug_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    if(DEBUG_LVL <= 1)
        return;

    va_list arglist;

    fprintf(stdout,"DEBUG (%s:%d): ", file, line);
    va_start(arglist, format);
    vfprintf(stdout, format, arglist);
    va_end(arglist);
    printf("\n");
}

void stdl_warning_msg(char *file, int line, char *format, ...) {
    assert(file != NULL && format != NULL);

    if(DEBUG_LVL < 1)
        return;

    va_list arglist;

    fprintf(stdout, "WARN (%s:%d): ", file, line);
    va_start(arglist, format);
    vfprintf(stdout, format, arglist);
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
