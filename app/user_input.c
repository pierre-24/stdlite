#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <toml.h>

#include <stdlite/helpers.h>

#include "user_input.h"

int stdl_user_input_new(stdl_user_input** inp_ptr) {
    assert(inp_ptr != NULL);

    *inp_ptr = malloc(sizeof(stdl_user_input));
    STDL_ERROR_HANDLE_AND_REPORT(*inp_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*inp_ptr)->title = NULL;

    // context
    (*inp_ptr)->ctx_source_path = NULL;
    (*inp_ptr)->ctx_source_type = STDL_SRC_CTX;

    (*inp_ptr)->ctx_output = malloc(128 * sizeof(char ));
    STDL_ERROR_HANDLE_AND_REPORT((*inp_ptr)->ctx_output == NULL, stdl_user_input_delete(*inp_ptr); return STDL_ERR_MALLOC, "malloc");
    strcpy((*inp_ptr)->ctx_output, "context.h5\0");

    (*inp_ptr)->ctx_method = STDL_METHOD_MONOPOLE;
    (*inp_ptr)->ctx_use_tda = 0;
    (*inp_ptr)->ctx_gammaJ = 4.f;
    (*inp_ptr)->ctx_gammaK = 2.f;
    (*inp_ptr)->ctx_ethr = 10.f / STDL_CONST_AU_TO_EV;
    (*inp_ptr)->ctx_e2thr = 1e-4f;
    (*inp_ptr)->ctx_ax = 0.5f;

    return STDL_ERR_OK;
}

int stdl_user_input_delete(stdl_user_input* inp) {
    assert(inp != NULL);

    STDL_FREE_ALL(inp->title, inp->ctx_source_path, inp->ctx_output, inp);

    return STDL_ERR_OK;
}

int stdl_user_input_fill_from_toml(stdl_user_input* inp, char* path) {
    assert(inp != NULL && path != NULL);

    stdl_log_msg(0, "Reading `%s` >", path);

    FILE* f = fopen(path, "r");
    STDL_ERROR_HANDLE_AND_REPORT(f == NULL, return STDL_ERR_OPEN, "cannot open `%s`", path);

    char errbuff[200];
    toml_table_t* conf = toml_parse_file(f, errbuff, sizeof(errbuff));
    fclose(f);

    STDL_ERROR_HANDLE_AND_REPORT(conf == NULL, return STDL_ERR_READ, "error while parsing conf: %s", errbuff);

    int err = STDL_ERR_OK;

    // level 0
    toml_datum_t title = toml_string_in(conf, "title");
    if(title.ok) {
        STDL_DEBUG("- title");
        inp->title = title.u.s;
    }

    // context
    toml_table_t* ctx = toml_table_in(conf, "context");
    if(ctx != NULL) {
        toml_datum_t ctx_source = toml_string_in(ctx, "source");
        if(ctx_source.ok) {
            inp->ctx_source_path = ctx_source.u.s;
            STDL_DEBUG("- context.source");
        }

        toml_datum_t ctx_source_type = toml_string_in(ctx, "source_type");
        if(ctx_source_type.ok) {
            STDL_DEBUG("- context.source_type");

            if(strcmp(ctx_source_type.u.s, "FCHK") == 0) {
                inp->ctx_source_type = STDL_SRC_FCHK;
            } else if(strcmp(ctx_source_type.u.s, "MOLDEN") == 0) {
                inp->ctx_source_type = STDL_SRC_MOLDEN;
            } else if(strcmp(ctx_source_type.u.s, "STDL_CTX") == 0) {
                inp->ctx_source_type = STDL_SRC_CTX;
            } else {
                STDL_ERROR_HANDLE_AND_REPORT(1, free(ctx_source_type.u.s); goto _end, "unknown value for `context.source_type`");
            }

            free(ctx_source_type.u.s);
        }

        toml_datum_t ctx_output = toml_string_in(ctx, "output");
        if(ctx_output.ok) {
            STDL_DEBUG("- context.output");
            STDL_FREE_IF_USED(inp->ctx_output);

            inp->ctx_output = ctx_output.u.s;
        }

        toml_datum_t ctx_method = toml_string_in(ctx, "method");
        if(ctx_method.ok) {
            STDL_DEBUG("- context.method");
            if(strcmp(ctx_method.u.s, "monopole") == 0) {
                inp->ctx_method = STDL_METHOD_MONOPOLE;
            } else if(strcmp(ctx_method.u.s, "monopole_direct") == 0) {
                inp->ctx_method = STDL_METHOD_MONOPOLE_DIRECT;
            } else {
                STDL_ERROR_HANDLE_AND_REPORT(1, free(ctx_method.u.s); goto _end, "unknown value for `context.method`");
            }

            free(ctx_method.u.s);
        }

        toml_datum_t  ctx_use_tda = toml_bool_in(ctx, "use_tda");
        if(ctx_use_tda.ok) {
            STDL_DEBUG("- context.use_tda");
            inp->ctx_use_tda = ctx_use_tda.u.b;
        }

        toml_datum_t ctx_gammaJ = toml_double_in(ctx, "gammaJ");
        if(ctx_gammaJ.ok) {
            STDL_DEBUG("- context.gammaJ");
            inp->ctx_gammaJ = (float) ctx_gammaJ.u.d;
        }

        toml_datum_t ctx_gammaK = toml_double_in(ctx, "gammaK");
        if(ctx_gammaK.ok) {
            STDL_DEBUG("- context.gammaK");
            inp->ctx_gammaK = (float) ctx_gammaK.u.d;
        }

        toml_datum_t ctx_ax = toml_double_in(ctx, "ax");
        if(ctx_ax.ok) {
            STDL_DEBUG("- context.ax");
            inp->ctx_ax = (float) ctx_ax.u.d;
        }
    }

    stdl_log_msg(0, "-");
    toml_table_t* response = toml_table_in(conf, "response");
    if(response != NULL) {
        // TODO: that.
    }

    stdl_log_msg(0, "< done\n");

    _end:
    toml_free(conf);
    return err;
}
