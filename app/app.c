#include "app.h"

int stdl_app_user_input(int argc, char* argv[], stdl_user_input** inp) {
    int err;

    err = stdl_user_input_new(inp);
    STDL_ERROR_CODE_HANDLE(err, return err);

    err = stdl_user_input_fill_from_args(*inp, argc, argv);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_delete(*inp); *inp = NULL; return err);

    err = stdl_user_input_check(*inp);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_delete(*inp); *inp = NULL; return err);

    return STDL_ERR_OK;
}
