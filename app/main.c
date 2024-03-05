#include <stdlib.h>
#include <string.h>

#include "app.h"
#include "timer.h"

#define NDASHES 80

void title(char* title) {
    int sz = strlen(title);
    stdl_log_msg(0, "\n--- %s -", title);
    for (int i = 0; i < NDASHES - sz - 6; ++i) {
        stdl_log_msg(0, "-");
    }
    stdl_log_msg(0, "\n");
}

int main(int argc, char* argv[]) {
    struct timespec elapsed_prog;
    timer_start(&elapsed_prog);

    int err;
    stdl_user_input* input = NULL;

    // read input
    err = stdl_user_input_new(&input);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    err = stdl_user_input_fill_from_args(input, argc, argv);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    // scream aloud
    stdl_log_msg(0,
                 "==========================================\n"
                 "Running %s\n"
                 "using %s (version %s)\n"
                 "==========================================\n",
                 APP_NAME, stdl_library_name(), stdl_library_version()
                 );

    _end:
    if(input != NULL)
        stdl_user_input_delete(input);

    if(err < STDL_ERR_LAST) {
        title("End");
        stdl_log_msg(0, "Elapsed time in %s: %.2f secs", APP_NAME, timer_stop(&elapsed_prog));
    }

    if(err == STDL_ERR_OK || err > STDL_ERR_LAST) {
        return EXIT_SUCCESS;
    } else {
        stdl_log_msg(0, "\n\n** Something went bad :(");
        return err;
    }
}
