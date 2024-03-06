#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <toml.h>
#include <argtable3.h>

#include <stdlite/helpers.h>

#include "user_input.h"

int stdl_user_input_new(stdl_user_input** inp_ptr) {
    assert(inp_ptr != NULL);

    *inp_ptr = malloc(sizeof(stdl_user_input));
    STDL_ERROR_HANDLE_AND_REPORT(*inp_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*inp_ptr)->title = NULL;

    // context
    (*inp_ptr)->ctx_source = NULL;
    (*inp_ptr)->ctx_source_type = STDL_SRC_CTX;

    (*inp_ptr)->ctx_output = malloc(128 * sizeof(char ));
    STDL_ERROR_HANDLE_AND_REPORT((*inp_ptr)->ctx_output == NULL, stdl_user_input_delete(*inp_ptr); return STDL_ERR_MALLOC, "malloc");
    strcpy((*inp_ptr)->ctx_output, "context.h5\0");

    // defaults of stda:
    (*inp_ptr)->ctx_method = STDL_METHOD_MONOPOLE;
    (*inp_ptr)->ctx_tda = 1;
    (*inp_ptr)->ctx_gammaJ = 4.f;
    (*inp_ptr)->ctx_gammaK = 2.f;
    (*inp_ptr)->ctx_ethr = 7.f / STDL_CONST_AU_TO_EV;
    (*inp_ptr)->ctx_e2thr = 1e-4f;
    (*inp_ptr)->ctx_ax = 0.5f;

    return STDL_ERR_OK;
}

int stdl_user_input_delete(stdl_user_input* inp) {
    assert(inp != NULL);

    STDL_FREE_ALL(inp->title, inp->ctx_source, inp->ctx_output, inp);

    return STDL_ERR_OK;
}

int _frequency_float_in(toml_table_t* table, char* field, float* result) {
    toml_datum_t freq = toml_double_in(table, field);
    if(freq.ok) {
        STDL_DEBUG("- (frequency) %s", field);
        *result = (float) freq.u.d;
    } else {
        freq = toml_string_in(table, field);
        if(freq.ok) {
            STDL_DEBUG("- (frequency) %s", field);
            double val;
            int err = stdl_user_input_parse_frequency(freq.u.s, &val);
            free(freq.u.s);
            STDL_ERROR_CODE_HANDLE(err, return err);

            *result = (float) val;
        }
    }

    return STDL_ERR_OK;
}

int stdl_user_input_fill_from_toml(stdl_user_input* inp, char* path) {
    assert(inp != NULL && path != NULL);

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
        STDL_DEBUG("Read [context]");

        toml_datum_t ctx_source = toml_string_in(ctx, "source");
        if(ctx_source.ok) {
            inp->ctx_source = ctx_source.u.s;
            STDL_DEBUG("- source");
        }

        toml_datum_t ctx_source_type = toml_string_in(ctx, "source_type");
        if(ctx_source_type.ok) {
            STDL_DEBUG("- source_type");

            if(strcmp(ctx_source_type.u.s, "FCHK") == 0) {
                inp->ctx_source_type = STDL_SRC_FCHK;
            } else if(strcmp(ctx_source_type.u.s, "MOLDEN") == 0) {
                inp->ctx_source_type = STDL_SRC_MOLDEN;
            } else if(strcmp(ctx_source_type.u.s, "STDL_CTX") == 0) {
                inp->ctx_source_type = STDL_SRC_CTX;
            } else if(strcmp(ctx_source_type.u.s, "STDL_CTX_WB") == 0) {
                inp->ctx_source_type = STDL_SRC_CTX_WB;
            } else {
                STDL_ERROR_HANDLE_AND_REPORT(1, err = STDL_ERR_INPUT; free(ctx_source_type.u.s); goto _end, "unknown value for `context.source_type`");
            }

            free(ctx_source_type.u.s);
        }

        toml_datum_t ctx_output = toml_string_in(ctx, "output");
        if(ctx_output.ok) {
            STDL_DEBUG("- output");
            STDL_FREE_IF_USED(inp->ctx_output);

            inp->ctx_output = ctx_output.u.s;
        }

        toml_datum_t ctx_method = toml_string_in(ctx, "method");
        if(ctx_method.ok) {
            STDL_DEBUG("- method");
            if(strcmp(ctx_method.u.s, "monopole") == 0) {
                inp->ctx_method = STDL_METHOD_MONOPOLE;
            } else if(strcmp(ctx_method.u.s, "monopole_direct") == 0) {
                inp->ctx_method = STDL_METHOD_MONOPOLE_DIRECT;
            } else {
                STDL_ERROR_HANDLE_AND_REPORT(1, free(ctx_method.u.s); goto _end, "unknown value for `context.method`");
            }

            free(ctx_method.u.s);
        }

        toml_datum_t  ctx_tda = toml_bool_in(ctx, "tda");
        if(ctx_tda.ok) {
            STDL_DEBUG("- tda");
            inp->ctx_tda = ctx_tda.u.b;
        }

        toml_datum_t ctx_gammaJ = toml_double_in(ctx, "gammaJ");
        if(ctx_gammaJ.ok) {
            STDL_DEBUG("- gammaJ");
            inp->ctx_gammaJ = (float) ctx_gammaJ.u.d;
        }

        toml_datum_t ctx_gammaK = toml_double_in(ctx, "gammaK");
        if(ctx_gammaK.ok) {
            STDL_DEBUG("- gammaK");
            inp->ctx_gammaK = (float) ctx_gammaK.u.d;
        }

        err = _frequency_float_in(ctx, "ethr", &(inp->ctx_ethr));
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        err = _frequency_float_in(ctx, "e2thr", &(inp->ctx_e2thr));
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        toml_datum_t ctx_ax = toml_double_in(ctx, "ax");
        if(ctx_ax.ok) {
            STDL_DEBUG("- ax");
            inp->ctx_ax = (float) ctx_ax.u.d;
        }
    }

    toml_table_t* response = toml_table_in(conf, "responses");
    if(response != NULL) {
        STDL_DEBUG("Read [responses]");
        // TODO: that.
    }

    _end:
    toml_free(conf);
    return err;
}

int stdl_user_input_parse_frequency(char* input, double* result) {
    assert(input != NULL && result != NULL);

    char *endptr;
    double value = strtod(input, &endptr);
    STDL_ERROR_HANDLE_AND_REPORT(endptr == input, return STDL_ERR_READ, "incorrect number in `%s`", input);

    if(strcmp(endptr, "eV") == 0 || strcmp(endptr, "ev") == 0) {
        *result = value / STDL_CONST_AU_TO_EV;
    } else if(strcmp(endptr, "nm") == 0) {
        *result = STDL_CONST_HC / value;
    } else if(strcmp(endptr, "au") == 0) {
        *result = value;
    } else {
        STDL_ERROR_HANDLE_AND_REPORT(1, return STDL_ERR_READ, "incorrect unit `%s`", endptr);
    }

    return STDL_ERR_OK;
}

int stdl_user_input_fill_from_args(stdl_user_input* inp, int argc, char* argv[]) {
    assert(inp != NULL && argc > 0 && argv != NULL);

    char* self = argv[0];
    struct arg_lit* arg_help;
    struct arg_end* arg_end_;
    struct arg_file* arg_input, *arg_ctx_source, *arg_ctx_output;
    struct arg_dbl* arg_ctx_gammaJ, *arg_ctx_gammaK, *arg_ctx_ax;
    struct arg_str* arg_ctx_source_type, *arg_ctx_ethr, *arg_ctx_e2thr;
    struct arg_int* arg_ctx_tda;

    int err = STDL_ERR_OK;

    void* args_table[] = {
            arg_help = arg_litn("h", "help", 0, 1, "display this help and exit"),
            arg_input = arg_file0(NULL, NULL, "<input>", "input file (TOML format)"),
            // context
            arg_ctx_source = arg_file0(NULL, "ctx_source", NULL, "source of wavefunction/basis"),
            arg_ctx_source_type = arg_str0(NULL, "ctx_source_type", "{FCHK,MOLDEN,STDL_CTX}", "type of source"),
            arg_ctx_output = arg_file0(NULL, "ctx_output", NULL, "output of context phase"),
            arg_ctx_gammaJ = arg_dbl0(NULL, "ctx_gammaJ", NULL, "gamma_J"),
            arg_ctx_gammaK = arg_dbl0(NULL, "ctx_gammaK", NULL, "gamma_K"),
            arg_ctx_ethr = arg_str0(NULL, "ctx_ethr", "<freq>", "ethr"),
            arg_ctx_e2thr = arg_str0(NULL, "ctx_e2thr", "<freq>", "e2thr"),
            arg_ctx_ax = arg_dbl0(NULL, "ctx_ax", NULL, "HF exchange"),
            arg_ctx_tda = arg_int0(NULL, "ctx_tda", "{0,1}", "use the Tamm-Dancoff approximation"),
            // end0
            arg_end_ = arg_end(20)
    };

    int nerrors = arg_parse(argc, argv, args_table);

    /* `--help` takes precedence over the rest */
    if (arg_help->count > 0) {
        printf("Usage: %s", self);
        arg_print_syntax(stdout, args_table, "\n");
        printf(STDL_APP_ARG_DOC);
        arg_print_glossary(stdout, args_table, "  %-40s %s\n");
        err = STDL_ERR_LAST + 1;
        goto _end;
    }

    /* report errors if any */
    if (nerrors > 0) {
        arg_print_errors(stderr, arg_end_, self);
        err = STDL_ERR_INPUT;
        goto _end;
    }

    // read input file
    if(arg_input->count > 0) {
        err = stdl_user_input_fill_from_toml(inp, arg_input->filename[0]);
        STDL_ERROR_CODE_HANDLE(err, goto _end);
    }

    // modify context
    size_t sz;
    double val;

    if(arg_ctx_source->count > 0) {
        sz = strlen(arg_ctx_source->filename[0]);
        inp->ctx_source = realloc(inp->ctx_source, (sz + 1) * sizeof(char ));
        STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_source == NULL, goto _end, "malloc");
        strcpy(inp->ctx_source, arg_ctx_source->filename[0]);
    }

    if(arg_ctx_source_type->count > 0) {
        if(strcmp(arg_ctx_source_type->sval[0], "FCHK") == 0) {
            inp->ctx_source_type = STDL_SRC_FCHK;
        } else if(strcmp(arg_ctx_source_type->sval[0], "MOLDEN") == 0) {
            inp->ctx_source_type = STDL_SRC_MOLDEN;
        } else if(strcmp(arg_ctx_source_type->sval[0], "STDL_CTX") == 0) {
            inp->ctx_source_type = STDL_SRC_CTX;
        } else if(strcmp(arg_ctx_source_type->sval[0], "STDL_CTX_WB") == 0) {
            inp->ctx_source_type = STDL_SRC_CTX_WB;
        } else {
            STDL_ERROR_HANDLE_AND_REPORT(1, err = STDL_ERR_INPUT; goto _end, "unknown source type `%s`", arg_ctx_source_type->sval[0]);
        }
    }

    if(arg_ctx_output->count > 0) {
        sz = strlen(arg_ctx_output->filename[0]);
        inp->ctx_output = realloc(inp->ctx_output, (sz + 1) * sizeof(char ));
        STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_output == NULL, goto _end, "malloc");
        strcpy(inp->ctx_output, arg_ctx_output->filename[0]);
    }

    if(arg_ctx_gammaJ->count > 0)
        inp->ctx_gammaJ = (float) arg_ctx_gammaJ->dval[0];

    if(arg_ctx_gammaK->count > 0)
        inp->ctx_gammaK = (float) arg_ctx_gammaK->dval[0];

    if(arg_ctx_ethr->count > 0) {
        err = stdl_user_input_parse_frequency(arg_ctx_ethr->sval[0], &val);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        inp->ctx_ethr = (float) val;
    }

    if(arg_ctx_e2thr->count > 0) {
        err = stdl_user_input_parse_frequency(arg_ctx_e2thr->sval[0], &val);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        inp->ctx_e2thr = (float) val;
    }

    if(arg_ctx_ax->count > 0)
        inp->ctx_ax = (float) arg_ctx_ax->dval[0];

    if(arg_ctx_tda->count > 0)
        inp->ctx_tda = arg_ctx_tda->ival[0];

    _end:
    arg_freetable(args_table, sizeof(args_table) / sizeof(args_table[0]));
    return err;
}


int stdl_user_input_check(stdl_user_input* inp) {
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_source == NULL, return STDL_ERR_INPUT, "missing context.source");

    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_gammaJ < .0, return STDL_ERR_INPUT, "context.gammaJ < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_gammaK < .0, return STDL_ERR_INPUT, "context.gammaK < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ethr < .0, return STDL_ERR_INPUT, "context.ethr < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_e2thr < .0, return STDL_ERR_INPUT, "context.e2thr < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ax < .0, return STDL_ERR_INPUT, "context.ax < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ax > 1, return STDL_ERR_INPUT, "context.ax > 1");

    return STDL_ERR_OK;
}

int stdl_user_input_log(stdl_user_input* inp) {

    if(inp->title != NULL)
        stdl_log_msg(0, "title = \"%s\"\n", inp->title);

    stdl_log_msg(0, "[context]\n");

    stdl_log_msg(0, "source = \"%s\"\n", inp->ctx_source);

    switch (inp->ctx_source_type) {
        case STDL_SRC_CTX:
            stdl_log_msg(0, "source_type = \"STDL_CTX\"\n");
            break;
        case STDL_SRC_CTX_WB:
            stdl_log_msg(0, "source_type = \"STDL_CTX_WB\"\n");
            break;
        case STDL_SRC_FCHK:
            stdl_log_msg(0, "source_type = \"FCHK\"\n");
            break;
        case STDL_SRC_MOLDEN:
            stdl_log_msg(0, "source_type = \"MOLDEN\"\n");
            break;
    }

    if(inp->ctx_source_type != STDL_SRC_CTX) {
        stdl_log_msg(0, "# parameters to build A' and B':\n");

        switch (inp->ctx_method) {
            case STDL_METHOD_MONOPOLE:
                stdl_log_msg(0, "method = \"monopole\"\n");
                break;
            case STDL_METHOD_MONOPOLE_DIRECT:
                stdl_log_msg(0, "method = \"monopole_direct\"\n");
                break;
        }

        stdl_log_msg(0, "use_tda = %s\n", inp->ctx_tda ? "true" : "false");
        stdl_log_msg(0, "gammaJ = %f\ngammaK = %f\nax = %f\n", inp->ctx_gammaJ, inp->ctx_gammaK, inp->ctx_ax);
        stdl_log_msg(0, "ethr = %f # au\ne2thr = %e # au\n", inp->ctx_ethr, inp->ctx_e2thr);

        stdl_log_msg(0, "output = \"%s\"\n", inp->ctx_output);
    } else {
        stdl_log_msg(0, "# remaining parameters will be obtained from %s\n", inp->ctx_source);
    }


    return STDL_ERR_OK;
}
