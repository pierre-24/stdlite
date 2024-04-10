#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <toml.h>
#include <argtable3.h>

#include <stdlite/helpers.h>
#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/molden_parser.h>
#include <stdlite/property.h>
#include <stdlite/response.h>
#include <stdlite/utils/matrix.h>
#include <stdlite/utils/experimental_quantity.h>

#include "user_input_handler.h"
#include "log_property.h"

int stdl_user_input_handler_new(stdl_user_input_handler** inp_ptr) {
    assert(inp_ptr != NULL);

    *inp_ptr = malloc(sizeof(stdl_user_input_handler));
    STDL_ERROR_HANDLE_AND_REPORT(*inp_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create user_input %p", *inp_ptr);

    (*inp_ptr)->title = NULL;

    // --- context:
    (*inp_ptr)->ctx_source = NULL;
    (*inp_ptr)->ctx_source_type = STDL_SRC_CTX;

    (*inp_ptr)->data_output = malloc(128 * sizeof(char ));
    STDL_ERROR_HANDLE_AND_REPORT((*inp_ptr)->data_output == NULL, stdl_user_input_handler_delete(*inp_ptr); return STDL_ERR_MALLOC, "malloc");
    strcpy((*inp_ptr)->data_output, "stdlite_calculation.h5\0");

    // defaults of stda:
    (*inp_ptr)->ctx_method = STDL_METHOD_MONOPOLE;
    (*inp_ptr)->ctx_tda = 1;
    (*inp_ptr)->ctx_gammaJ = 4.f;
    (*inp_ptr)->ctx_gammaK = 2.f;
    (*inp_ptr)->ctx_ethr = 7.f / STDL_CONST_AU_TO_EV;
    (*inp_ptr)->ctx_e2thr = 1e-4f;
    (*inp_ptr)->ctx_ax = 0.5f;

    // -- response:
    (*inp_ptr)->res_resreqs = NULL;

    return STDL_ERR_OK;
}

int stdl_user_input_handler_delete(stdl_user_input_handler* inp) {
    assert(inp != NULL);

    STDL_DEBUG("delete user_input %p", inp);

    if(inp->res_resreqs != NULL)
        stdl_response_request_delete(inp->res_resreqs);

    STDL_FREE_ALL(inp->title, inp->ctx_source, inp->data_output, inp);

    return STDL_ERR_OK;
}

int _energy_in(toml_table_t* table, char* field, float* result, int* isset) {
    toml_datum_t energy = toml_double_in(table, field);
    *isset = 0;

    if(energy.ok) {
        STDL_DEBUG("- (energy) %s", field);
        *result = (float) energy.u.d;
        *isset = 1;
    } else {
        energy = toml_string_in(table, field);
        if(energy.ok) {
            STDL_DEBUG("- (energy) %s", field);
            double val;
            int err = stdl_user_input_handler_parse_frequency(energy.u.s, &val);
            free(energy.u.s);
            STDL_ERROR_CODE_HANDLE(err, return err);

            *result = (float) val;
            *isset = 1;
        }
    }

    return STDL_ERR_OK;
}

int _operator_in(toml_table_t* table, char* field, stdl_operator * result, int* isset) {
    toml_datum_t op = toml_string_in(table, field);
    int err = STDL_ERR_OK;
    *isset = 0;

    if(op.ok) {
        STDL_DEBUG("- (operator) %s", field);
        if(strcmp(op.u.s, "dipl") == 0) {
            *result = STDL_OP_DIPL;
            *isset = 1;
        } else {
            err = STDL_ERR_INPUT;
        }

        free(op.u.s);
    }

    return err;
}

int stdl_user_input_handler_fill_from_toml(stdl_user_input_handler* inp, FILE *f) {
    assert(inp != NULL && f != NULL);

    char errbuff[200];
    int isset;
    toml_table_t* conf = toml_parse_file(f, errbuff, sizeof(errbuff));

    STDL_ERROR_HANDLE_AND_REPORT(conf == NULL, return STDL_ERR_READ, "error while parsing conf: %s", errbuff);

    int err = STDL_ERR_OK;

    // level 0
    toml_datum_t title = toml_string_in(conf, "title");
    if(title.ok) {
        STDL_DEBUG("- title");
        inp->title = title.u.s;
    }

    toml_datum_t ctx_output = toml_string_in(conf, "data_output");
    if(ctx_output.ok) {
        STDL_DEBUG("- data output");
        STDL_FREE_IF_USED(inp->data_output);

        inp->data_output = ctx_output.u.s;
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

        toml_datum_t ctx_method = toml_string_in(ctx, "method");
        if(ctx_method.ok) {
            STDL_DEBUG("- method");
            if(strcmp(ctx_method.u.s, "monopole") == 0) {
                inp->ctx_method = STDL_METHOD_MONOPOLE;
            } else {
                STDL_ERROR_HANDLE_AND_REPORT(1, free(ctx_method.u.s); goto _end, "unknown value for `context.method`");
            }

            free(ctx_method.u.s);
        }

        toml_datum_t  ctx_tda = toml_int_in(ctx, "tda");
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

        err = _energy_in(ctx, "ethr", &(inp->ctx_ethr), &isset);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        err = _energy_in(ctx, "e2thr", &(inp->ctx_e2thr), &isset);
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

        stdl_response_request* prev = NULL;

        // linear response
        toml_array_t* lresp = toml_array_in(response, "linear");
        if(lresp != NULL) {
            STDL_DEBUG("- linear");
            int i = 0;
            toml_table_t* t = toml_table_at(lresp, i);
            while (t != NULL) {
                STDL_DEBUG("linear[%d]", i);

                stdl_operator opA, opB;
                err = _operator_in(t, "opA", &opA, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `opA` in `linear[%d]`", i);

                err = _operator_in(t, "opB", &opB, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `opB` in `linear[%d]`", i);

                float w;
                err = _energy_in(t, "wB", &w, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `w` in `linear[%d]`", i);

                stdl_response_request* req = NULL;
                err = stdl_response_request_new(1, 0, (stdl_operator[]) {opA, opB}, (float[]) {w, w}, 0, &req);
                STDL_ERROR_CODE_HANDLE(err, goto _end);

                if(prev == NULL) {
                    inp->res_resreqs = req;
                } else {
                    prev->next = req;
                }

                prev = req;

                i += 1;
                t = toml_table_at(lresp, i);
            }
        }

        // quadratic response
        toml_array_t* qresp = toml_array_in(response, "quadratic");
        if(qresp != NULL) {
            STDL_DEBUG("- quadratic");
            int i = 0;
            toml_table_t* t = toml_table_at(qresp, i);
            while (t != NULL) {
                STDL_DEBUG("quadratic[%d]", i);

                stdl_operator opA, opB, opC;
                err = _operator_in(t, "opA", &opA, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `opA` in `quadratic[%d]`", i);

                err = _operator_in(t, "opB", &opB, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `opB` in `quadratic[%d]`", i);

                err = _operator_in(t, "opC", &opC, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `opC` in `quadratic[%d]`", i);

                float wB, wC;
                err = _energy_in(t, "wB", &wB, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `wB` in `quadratic[%d]`", i);

                err = _energy_in(t, "wC", &wC, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `wC` in `quadratic[%d]`", i);

                stdl_response_request* req = NULL;
                err = stdl_response_request_new(2, 0, (stdl_operator[]) {opA, opB, opC}, (float[]) {wB + wC, wB, wC}, 0, &req);
                STDL_ERROR_CODE_HANDLE(err, goto _end);

                if(prev == NULL) {
                    inp->res_resreqs = req;
                } else {
                    prev->next = req;
                }

                prev = req;

                i += 1;
                t = toml_table_at(qresp, i);
            }
        }

        // single residue of the linear response
        toml_array_t* lresp_sr = toml_array_in(response, "linear_sr");
        if(lresp_sr != NULL) {
            STDL_DEBUG("- linear_sr");
            int i = 0;
            toml_table_t* t = toml_table_at(lresp_sr, i);
            while (t != NULL) {
                STDL_DEBUG("lresp_sr[%d]", i);

                stdl_operator opA;
                err = _operator_in(t, "opA", &opA, &isset);
                STDL_ERROR_CODE_HANDLE(err, goto _end);
                STDL_ERROR_HANDLE_AND_REPORT(!isset, err = STDL_ERR_INPUT; goto _end, "missing `op` in `linear_sr[%d]`", i);

                int nroots;
                toml_datum_t root = toml_int_in(t, "nroots");
                STDL_ERROR_HANDLE_AND_REPORT(!root.ok, err = STDL_ERR_INPUT; goto _end, "missing `nroots` in `linear_sr[%d]`", i);
                nroots = (int) root.u.i;

                stdl_response_request* req = NULL;
                err = stdl_response_request_new(1, 1, (stdl_operator[]) {opA}, NULL, nroots, &req);
                STDL_ERROR_CODE_HANDLE(err, goto _end);

                if(prev == NULL) {
                    inp->res_resreqs = req;
                } else {
                    prev->next = req;
                }

                prev = req;

                i += 1;
                t = toml_table_at(lresp_sr, i);
            }
        }
    }

    _end:
    toml_free(conf);
    return err;
}

int stdl_user_input_handler_parse_frequency(char* input, double* result) {
    assert(input != NULL && result != NULL);

    char *endptr;
    double value = strtod(input, &endptr);
    STDL_ERROR_HANDLE_AND_REPORT(endptr == input, return STDL_ERR_READ, "incorrect number in `%s`", input);

    if(*endptr != '\0') {
        if(strcmp(endptr, "eV") == 0 || strcmp(endptr, "ev") == 0) {
            *result = value / STDL_CONST_AU_TO_EV;
        } else if(strcmp(endptr, "nm") == 0) {
            *result = STDL_CONST_HC / value;
        } else if(strcmp(endptr, "au") == 0) {
            *result = value;
        } else {
            STDL_ERROR_HANDLE_AND_REPORT(1, return STDL_ERR_READ, "incorrect unit `%s`", endptr);
        }
    } else
        *result = value;

    return STDL_ERR_OK;
}

int stdl_user_input_handler_fill_from_args(stdl_user_input_handler* inp, int argc, char* argv[]) {
    assert(inp != NULL && argc > 0 && argv != NULL);

    char* self = argv[0];
    struct arg_lit* arg_help;
    struct arg_end* arg_end_;
    struct arg_file* arg_input, *arg_ctx_source, *arg_data_output;
    struct arg_dbl* arg_ctx_gammaJ, *arg_ctx_gammaK, *arg_ctx_ax;
    struct arg_str* arg_ctx_source_type, *arg_ctx_ethr, *arg_ctx_e2thr, *arg_ctx_method;
    struct arg_int* arg_ctx_tda;

    int err = STDL_ERR_OK;

    void* args_table[] = {
            arg_help = arg_litn("h", "help", 0, 1, "display this help and exit"),
            arg_input = arg_file0(NULL, NULL, "<input>", "input file (TOML format)"),
            // context
            arg_ctx_source = arg_file0(NULL, "ctx_source", NULL, "source of wavefunction/basis"),
            arg_ctx_source_type = arg_str0(NULL, "ctx_source_type", "{FCHK,MOLDEN,STDL_CTX}", "type of source"),
            arg_data_output = arg_file0(NULL, "data_output", NULL, "checkpoint file"),
            arg_ctx_gammaJ = arg_dbl0(NULL, "ctx_gammaJ", NULL, "gamma_J"),
            arg_ctx_gammaK = arg_dbl0(NULL, "ctx_gammaK", NULL, "gamma_K"),
            arg_ctx_ethr = arg_str0(NULL, "ctx_ethr", "<freq>", "ethr"),
            arg_ctx_e2thr = arg_str0(NULL, "ctx_e2thr", "<freq>", "e2thr"),
            arg_ctx_method = arg_str0(NULL, "ctx_method", "str", "method"),
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
        FILE* f = fopen(arg_input->filename[0], "r");
        STDL_ERROR_HANDLE_AND_REPORT(f == NULL, err = STDL_ERR_OPEN; goto _end, "cannot open `%s`", arg_input->filename[0]);

        err = stdl_user_input_handler_fill_from_toml(inp, f);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        fclose(f);
    }

    size_t sz;
    if(arg_data_output->count > 0) {
        sz = strlen(arg_data_output->filename[0]);
        inp->data_output = realloc(inp->data_output, (sz + 1) * sizeof(char ));
        STDL_ERROR_HANDLE_AND_REPORT(inp->data_output == NULL, goto _end, "malloc");
        strcpy(inp->data_output, arg_data_output->filename[0]);
    }

    // modify context
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

    if(arg_ctx_gammaJ->count > 0)
        inp->ctx_gammaJ = (float) arg_ctx_gammaJ->dval[0];

    if(arg_ctx_gammaK->count > 0)
        inp->ctx_gammaK = (float) arg_ctx_gammaK->dval[0];

    if(arg_ctx_ethr->count > 0) {
        err = stdl_user_input_handler_parse_frequency(arg_ctx_ethr->sval[0], &val);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        inp->ctx_ethr = (float) val;
    }

    if(arg_ctx_e2thr->count > 0) {
        err = stdl_user_input_handler_parse_frequency(arg_ctx_e2thr->sval[0], &val);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        inp->ctx_e2thr = (float) val;
    }

    if(arg_ctx_method->count > 0) {
        if(strcmp(arg_ctx_method->sval[0], "monopole") == 0) {
            inp->ctx_method = STDL_METHOD_MONOPOLE;
        } else {
            STDL_ERROR_HANDLE_AND_REPORT(1, err = STDL_ERR_INPUT; goto _end, "unknown method `%s`", arg_ctx_method->sval[0]);
        }
    }

    if(arg_ctx_ax->count > 0)
        inp->ctx_ax = (float) arg_ctx_ax->dval[0];

    if(arg_ctx_tda->count > 0)
        inp->ctx_tda = arg_ctx_tda->ival[0];

    _end:
    arg_freetable(args_table, sizeof(args_table) / sizeof(args_table[0]));
    return err;
}


int stdl_user_input_handler_check(stdl_user_input_handler* inp) {
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_source == NULL, return STDL_ERR_INPUT, "missing context.source");

    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_gammaJ < .0, return STDL_ERR_INPUT, "context.gammaJ < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_gammaK < .0, return STDL_ERR_INPUT, "context.gammaK < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ethr < .0, return STDL_ERR_INPUT, "context.ethr < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_e2thr < .0, return STDL_ERR_INPUT, "context.e2thr < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ax < .0, return STDL_ERR_INPUT, "context.ax < 0");
    STDL_ERROR_HANDLE_AND_REPORT(inp->ctx_ax > 1, return STDL_ERR_INPUT, "context.ax > 1");

    return STDL_ERR_OK;
}

int stdl_user_input_handler_new_from_args(int argc, char* argv[], stdl_user_input_handler** inp) {
    int err;

    err = stdl_user_input_handler_new(inp);
    STDL_ERROR_CODE_HANDLE(err, return err);

    err = stdl_user_input_handler_fill_from_args(*inp, argc, argv);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_handler_delete(*inp); *inp = NULL; return err);

    err = stdl_user_input_handler_check(*inp);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_handler_delete(*inp); *inp = NULL; return err);

    return STDL_ERR_OK;
}

void _op_log(char* name, stdl_operator op) {
    stdl_log_msg(0, "%s=\"", name);
    switch (op) {
        case STDL_OP_DIPL:
            stdl_log_msg(0, "dipl");
            break;
        default:
            stdl_log_msg(0, "unk");
    }

    stdl_log_msg(0, "\"");
}

int stdl_user_input_handler_log(stdl_user_input_handler* inp) {

    if(inp->title != NULL)
        stdl_log_msg(0, "title = \"%s\"\n", inp->title);

    stdl_log_msg(0, "data_output = \"%s\"\n", inp->data_output);

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
        }

        stdl_log_msg(0, "tda = %d\n", inp->ctx_tda);
        stdl_log_msg(0, "gammaJ = %f\ngammaK = %f\nax = %f\n", inp->ctx_gammaJ, inp->ctx_gammaK, inp->ctx_ax);
        stdl_log_msg(0, "ethr = %f # au\ne2thr = %e # au\n", inp->ctx_ethr, inp->ctx_e2thr);
    } else {
        stdl_log_msg(0, "# A' and B' are obtained from %s\n", inp->ctx_source);
    }

    if(inp->res_resreqs != NULL) {
        stdl_log_msg(0, "[responses]\n");

        stdl_response_request* req = NULL;

        // linear
        stdl_log_msg(0, "linear = [\n");
        req = inp->res_resreqs;
        while(req != NULL) {
            if(req->resp_order == 1 && req->res_order == 0) {
                stdl_log_msg(0, "  {");
                _op_log("opA", req->ops[0]);
                stdl_log_msg(0, ", ");
                _op_log("opB", req->ops[1]);
                stdl_log_msg(0, ", wB=%f},\n", req->w[1]);
            }
            req = req->next;
        }
        stdl_log_msg(0, "]\n");

        // quadratic
        stdl_log_msg(0, "quadratic = [\n");
        req = inp->res_resreqs;
        while(req != NULL) {
            if(req->resp_order == 2 && req->res_order == 0) {
                stdl_log_msg(0, "  {");
                _op_log("opA", req->ops[0]);
                stdl_log_msg(0, ", ");
                _op_log("opB", req->ops[1]);
                stdl_log_msg(0, ", ");
                _op_log("opC", req->ops[2]);
                stdl_log_msg(0, ", wB=%f, wC=%f},\n", req->w[1], req->w[2]);
            }
            req = req->next;
        }
        stdl_log_msg(0, "]\n");

        // linear_sr
        stdl_log_msg(0, "linear_sr = [\n");
        req = inp->res_resreqs;
        while(req != NULL) {
            if(req->resp_order == 1 && req->res_order == 1) {
                stdl_log_msg(0, "  {");
                _op_log("opA", req->ops[0]);
                stdl_log_msg(0, ", nroots=%d},\n", req->nroots);
            }
            req = req->next;
        }
        stdl_log_msg(0, "]\n");
    }

    return STDL_ERR_OK;
}


int _from_h5(char* path, stdl_wavefunction** wf_ptr, stdl_basis** bs_ptr) {
    hid_t file_id;

    // open
    file_id = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", path);

    // read
    int err = stdl_wavefunction_load_h5(file_id, wf_ptr);
    STDL_ERROR_CODE_HANDLE(err, H5Fclose(file_id); return err);

    err = stdl_basis_load_h5(file_id, bs_ptr);

    // close
    H5Fclose(file_id);
    STDL_ERROR_CODE_HANDLE(err, stdl_wavefunction_delete(*wf_ptr); return err);

    size_t nprim = 0;
    for (int ibas = 0; ibas < (*bs_ptr)->nbas; ++ibas) {
        nprim += (*bs_ptr)->bas[ibas * 8 + 2];
    }

    stdl_log_msg(0, "Got %d atoms, %d AOs (%d primitives in %d basis functions), and %d MOs\n", (*wf_ptr)->natm, (*wf_ptr)->nao, nprim, (*bs_ptr)->nbas, (*wf_ptr)->nmo);

    return STDL_ERR_OK;
}

int _from_h5_full(char* path, stdl_context** ctx_ptr) {
    hid_t file_id;

    // open
    file_id = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", path);

    // read
    int err = stdl_context_load_h5(file_id, ctx_ptr);

    // close
    H5Fclose(file_id);

    return err;
}

int stdl_user_input_handler_make_context(stdl_user_input_handler* inp, stdl_context **ctx_ptr) {
    stdl_wavefunction* wf = NULL;
    stdl_basis* bs = NULL;
    int err;

    if(inp->ctx_source_type == STDL_SRC_CTX) { // directly read context from a previous calculation
        err = _from_h5_full(inp->ctx_source, ctx_ptr);
        STDL_ERROR_CODE_HANDLE(err, return err);
    } else { // read file for wavefunction & basis, then select normally
        if(inp->ctx_source_type == STDL_SRC_CTX_WB) {
            err = _from_h5(inp->ctx_source, &wf, &bs);
            STDL_ERROR_CODE_HANDLE(err, return err);
        } else { // MOLDEN or FCHK, requires lexer
            stdl_lexer* lx = NULL;
            FILE* f = fopen(inp->ctx_source, "r");
            STDL_ERROR_HANDLE_AND_REPORT(f == NULL, return STDL_ERR_OPEN, "cannot open `%s`", inp->ctx_source);

            err = stdl_lexer_new(f, &lx);
            STDL_ERROR_CODE_HANDLE(err, fclose(f); return err);

            if(inp->ctx_source_type == STDL_SRC_FCHK) {
                err = stdl_fchk_parser_skip_intro(lx) || stdl_fchk_parser_extract(lx, &wf, &bs);
            } else if(inp->ctx_source_type == STDL_SRC_MOLDEN)
                err = stdl_molden_parser_extract(lx, &wf, &bs);

            fclose(f);
            stdl_lexer_delete(lx);

            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        // create context
        err = stdl_context_new(wf, bs, inp->ctx_gammaJ, inp->ctx_gammaK, inp->ctx_ethr, inp->ctx_e2thr, inp->ctx_ax, ctx_ptr);
        STDL_ERROR_CODE_HANDLE(err, *ctx_ptr = NULL; return err);

        // select and build A' and B'
        if(inp->ctx_method == STDL_METHOD_MONOPOLE)
            err = stdl_context_select_csfs_monopole(*ctx_ptr, !inp->ctx_tda);

        STDL_ERROR_CODE_HANDLE(err, return err);

        // save context
        if(strcmp(inp->ctx_source, inp->data_output) == 0)
            STDL_WARN("`context.source` and `context.data_output` are the same, so the content of `%s` will be replaced", inp->data_output);

        hid_t file_id = H5Fcreate(inp->data_output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", inp->data_output)

        err = stdl_context_dump_h5(*ctx_ptr, file_id);
        H5Fclose(file_id);
    }

    return err;
}

struct _w_list {
    float w;
    struct _w_list* next;
};

int _w_list_new(float w, struct _w_list** elm) {
    *elm = malloc(sizeof(struct _w_list));
    STDL_ERROR_HANDLE_AND_REPORT(*elm == NULL, return STDL_ERR_MALLOC, "malloc");

    (*elm)->w = w;
    (*elm)->next = NULL;
    return STDL_ERR_OK;
}

int _w_list_delete(struct _w_list* lst) {
    assert(lst != NULL);

    if(lst->next != NULL)
        _w_list_delete(lst->next);

    STDL_FREE_IF_USED(lst);

    return STDL_ERR_OK;
}

int stdl_user_input_handler_prepare_responses(stdl_user_input_handler *inp, stdl_context *ctx, stdl_responses_handler **rh_ptr) {
    assert(inp != NULL && ctx != NULL && inp->res_resreqs != NULL && rh_ptr != NULL);

    int err;

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Preparing responses >");
    stdl_log_msg(1, "\n  | Count requests ");

    // count the number of operators, LRV requests, amplitudes, and freqs.
    size_t res_nops = 0, res_nlrvreq = 0, res_nexci = 0, totnw = 0;

    short operators[STDL_OP_COUNT] = {0};
    short islrvs[STDL_OP_COUNT] = {0};
    struct _w_list* lrvs_w[STDL_OP_COUNT] = {NULL};

    stdl_response_request* req = inp->res_resreqs;
    while (req != NULL) {
        // check out if it contains a new operator
        size_t nops = req->resp_order - req->res_order + 1;
        for (size_t iop = 0; iop < nops; ++iop) {
            stdl_operator op = req->ops[iop];
            if(!operators[op])
                res_nops += 1;

            if(!islrvs[op] && req->res_order < req->resp_order) {
                islrvs[op] = 1;
                res_nlrvreq += 1;
            }

            operators[op] = 1;
        }

        // if LRV is required, add frequency
        size_t nw = (req->resp_order == req->res_order)? 0: req->resp_order-req->res_order+1;
        for (size_t iw = 0; iw < nw; ++iw) {
            stdl_operator op = req->ops[iw];
            struct _w_list* elm = NULL;
            err = _w_list_new(req->w[iw], &elm);
            STDL_ERROR_CODE_HANDLE(err, return err);

            if(lrvs_w[op] == NULL) {
                lrvs_w[op] = elm;
            } else {
                struct _w_list* last = lrvs_w[op];
                int already_in = 0;
                while (last->next != NULL && !already_in) {
                    if(stdl_float_equals(req->w[iw], last->w, 1e-6f)) {
                        already_in = 1;
                    }

                    last = last->next;
                }

                if(!already_in && !stdl_float_equals(req->w[iw], last->w, 1e-6f)) {
                    last->next = elm;
                }
                else
                    _w_list_delete(elm);
            }
        }

        // check out if it requires amplitudes
        if(req->res_order > 0) {
            if(req->nroots < 0)
                res_nexci = ctx->ncsfs;
            else if((size_t) req->nroots > res_nexci) {
                if((size_t) req->nroots > ctx->ncsfs) {
                    STDL_WARN("%ld excited states requested, which is more than the number of CSFs", req->nroots);
                    req->nroots = (int) ctx->ncsfs;
                }
                res_nexci = (size_t) req->nroots;
            }
        }

        req = req->next;
    }

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | build requests ");

    err = stdl_responses_handler_new(res_nops, res_nlrvreq, res_nexci, inp->data_output, ctx, rh_ptr);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // copy operators
    int ioffset = 0;
    for (int iop = 0; iop < STDL_OP_COUNT; ++iop) {
        if(operators[iop]) {
            (*rh_ptr)->ops[ioffset] = iop;
            ioffset++;
        }
    }

    // create LRV requests
    stdl_lrv_request* lrvs[STDL_OP_COUNT] = {NULL};
    ioffset = 0;
    for (int iop = 0; iop < STDL_OP_COUNT; ++iop) {
        if(islrvs[iop]) {
            // count the number of frequencies
            size_t nw = 0;
            struct _w_list* last = lrvs_w[iop];
            while (last != NULL) {
                nw++;
                last = last->next;
            }

            totnw += nw;

            STDL_ERROR_HANDLE_AND_REPORT(nw == 0, return STDL_ERR_INPUT, "LRV but nw=0");

            // create LRV request
            (*rh_ptr)->lrvreqs[ioffset] = NULL;
            err = stdl_lrv_request_new(iop, nw, ctx->ncsfs, (*rh_ptr)->lrvreqs + ioffset);
            STDL_ERROR_CODE_HANDLE(err, return err);

            lrvs[iop] = (*rh_ptr)->lrvreqs[ioffset];

            // copy frequencies
            nw = 0;
            struct _w_list* curr = lrvs_w[iop];
            struct _w_list* prev = NULL;
            while (curr != NULL) {
                (*rh_ptr)->lrvreqs[ioffset]->w[nw] = curr->w;

                nw++;

                prev = curr;
                curr = curr->next;
                prev->next = NULL;
                _w_list_delete(prev);
            }

            ioffset++;
        }
    }

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | assign each response to its request ");

    req = inp->res_resreqs;
    while (req != NULL) {
        if(req->resp_order != req->res_order) {
            size_t nops = req->resp_order - req->res_order + 1;
            for (size_t iop = 0; iop < nops; ++iop) {
                stdl_lrv_request *lrvreq = lrvs[req->ops[iop]];
                req->lrvreqs[iop] = lrvreq;
                for (size_t jw = 0; jw < lrvreq->nw; ++jw) {
                    if (stdl_float_equals(req->w[iop], lrvreq->w[jw], 1e-6f)) {
                        req->wpos[iop] = jw;
                        break;
                    }
                }
            }
        }

        req = req->next;
    }

    stdl_log_msg(0, "< done\n");
    stdl_log_msg(0, "Will compute %ld matrix(ces) of MO integrals, %ld response vector(s), and %ld amplitude vector(s)\n", (*rh_ptr)->nops, totnw, (*rh_ptr)->nexci);

    return STDL_ERR_OK;
}

int stdl_user_input_handler_compute_properties(stdl_user_input_handler* inp, stdl_context* ctx, stdl_responses_handler* rh) {
    assert(inp != NULL && ctx != NULL && rh != NULL);

    int err;

    stdl_response_request* req = inp->res_resreqs;
    while (req != NULL) {

        if(req->resp_order == 1 && req->res_order == 0) { // linear
            size_t dim0 = STDL_OPERATOR_DIM[req->lrvreqs[0]->op], dim1 = STDL_OPERATOR_DIM[req->lrvreqs[1]->op];

            float* tensor = malloc(dim0 * dim1 * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(tensor == NULL, return STDL_ERR_MALLOC, "malloc");

            err = stdl_response_lr_tensor(
                    ctx,
                    (size_t[]) {dim0, dim1}, (int[]) {1, 1},
                    req->lrvreqs[0]->op_integrals,
                    req->lrvreqs[1]->X + req->wpos[1] * dim1 * ctx->ncsfs,
                    req->lrvreqs[1]->Y + req->wpos[1] * dim1 * ctx->ncsfs,
                    0, tensor);
            STDL_ERROR_CODE_HANDLE(err, free(tensor); return err);

            if(req->ops[0] == STDL_OP_DIPL && req->ops[1] == STDL_OP_DIPL)
                stdl_log_property_polarizability(req, tensor);

            free(tensor);

        } else if(req->resp_order == 2 && req->res_order == 0) { // quadratic
            size_t dim0 = STDL_OPERATOR_DIM[req->lrvreqs[0]->op],
                    dim1 = STDL_OPERATOR_DIM[req->lrvreqs[1]->op],
                    dim2 = STDL_OPERATOR_DIM[req->lrvreqs[2]->op];

            // TODO: it should be more general than that!
            float beta[3][3][3];
            stdl_property_first_hyperpolarizability(
                    ctx,
                    req->lrvreqs[1]->op_integrals,
                    (float*[]) {req->lrvreqs[0]->Y + req->wpos[0] * dim0 * ctx->ncsfs, req->lrvreqs[1]->X + req->wpos[1] * dim1 * ctx->ncsfs,req->lrvreqs[2]->X + req->wpos[2] * dim2 * ctx->ncsfs},
                    (float*[]) {req->lrvreqs[0]->X + req->wpos[0] * dim0 * ctx->ncsfs, req->lrvreqs[1]->Y + req->wpos[1] * dim1 * ctx->ncsfs,req->lrvreqs[2]->Y + req->wpos[2] * dim2 * ctx->ncsfs},
                    beta
                    );

            if(req->ops[0] == STDL_OP_DIPL && req->ops[1] == STDL_OP_DIPL && req->ops[2] == STDL_OP_DIPL)
                stdl_log_property_first_hyperpolarizability(req, beta);

        } else if(req->resp_order == 1 && req->res_order == 1) { // linear SR

            float* tdips = malloc(rh->nexci * 3 * sizeof(float ));
            STDL_ERROR_HANDLE_AND_REPORT(tdips == NULL, return STDL_ERR_MALLOC, "malloc");

            // TODO: it should be more general than that!
            stdl_property_transition_dipoles(ctx, rh->nexci, rh->ops_integrals[0] /* <- !!!! */, rh->Xamp, rh->Yamp, tdips);

            stdl_log_property_g2e_dipoles(rh, ctx, tdips, 5e-3f);

            STDL_FREE_IF_USED(tdips);
        }

        req = req->next;
    }

    return STDL_ERR_OK;
}

int stdl_user_input_handler_approximate_size(stdl_user_input_handler *inp, size_t *sz, size_t *respreq_sz) {
    assert(inp != NULL && sz != NULL);

    *respreq_sz = 0;

    if(inp->res_resreqs != NULL)
        stdl_response_request_approximate_size(inp->res_resreqs, respreq_sz);

    *sz = sizeof(stdl_user_input_handler)
            + *respreq_sz;

    return STDL_ERR_OK;
}
