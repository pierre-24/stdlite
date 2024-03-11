#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <toml.h>
#include <argtable3.h>

#include <stdlite/helpers.h>
#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/molden_parser.h>

#include "user_input.h"

int stdl_user_input_new(stdl_user_input** inp_ptr) {
    assert(inp_ptr != NULL);

    *inp_ptr = malloc(sizeof(stdl_user_input));
    STDL_ERROR_HANDLE_AND_REPORT(*inp_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    STDL_DEBUG("create user_input %p", *inp_ptr);

    (*inp_ptr)->title = NULL;

    // --- context:
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

    // -- response:
    (*inp_ptr)->res_nops = 0;
    (*inp_ptr)->res_ops = NULL;
    (*inp_ptr)->res_resreqs = NULL;
    (*inp_ptr)->res_nlrvreq = 0;
    (*inp_ptr)->res_lrvreqs = NULL;
    (*inp_ptr)->res_nexci = 0;
    (*inp_ptr)->res_eexci = NULL;
    (*inp_ptr)->res_Xamp = NULL;
    (*inp_ptr)->res_Yamp = NULL;

    return STDL_ERR_OK;
}

int stdl_user_input_delete(stdl_user_input* inp) {
    assert(inp != NULL);

    STDL_DEBUG("delete user_input %p", inp);

    if(inp->res_resreqs != NULL)
        stdl_response_request_delete(inp->res_resreqs);

    for (size_t i = 0; i < inp->res_nlrvreq; ++i) {
        stdl_lrv_request_delete(inp->res_lrvreqs[i]);
    }

    STDL_FREE_ALL(inp->title, inp->ctx_source, inp->ctx_output, inp->res_lrvreqs, inp->res_ops, inp->res_eexci, inp->res_Xamp, inp->res_Yamp, inp);

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
            int err = stdl_user_input_parse_frequency(energy.u.s, &val);
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

int stdl_user_input_fill_from_toml(stdl_user_input* inp, FILE *f) {
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
                err = stdl_response_request_new(1, 0, (stdl_operator[]) {opA, opB}, (float[]) {w}, 0, &req);
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
                err = stdl_response_request_new(2, 0, (stdl_operator[]) {opA, opB, opC}, (float[]) {wB, wC}, 0, &req);
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

                int nroot;
                toml_datum_t root = toml_int_in(t, "root");
                STDL_ERROR_HANDLE_AND_REPORT(!root.ok, err = STDL_ERR_INPUT; goto _end, "missing `root` in `linear_sr[%d]`", i);
                nroot = (int) root.u.i;

                stdl_response_request* req = NULL;
                err = stdl_response_request_new(1, 1, (stdl_operator[]) {opA}, NULL, nroot, &req);
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

int stdl_user_input_parse_frequency(char* input, double* result) {
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
        FILE* f = fopen(arg_input->filename[0], "r");
        STDL_ERROR_HANDLE_AND_REPORT(f == NULL, err = STDL_ERR_OPEN; goto _end, "cannot open `%s`", arg_input->filename[0]);

        err = stdl_user_input_fill_from_toml(inp, f);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        fclose(f);
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

int stdl_user_input_new_from_args(int argc, char* argv[], stdl_user_input** inp) {
    int err;

    err = stdl_user_input_new(inp);
    STDL_ERROR_CODE_HANDLE(err, return err);

    err = stdl_user_input_fill_from_args(*inp, argc, argv);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_delete(*inp); *inp = NULL; return err);

    err = stdl_user_input_check(*inp);
    STDL_ERROR_CODE_HANDLE(err, stdl_user_input_delete(*inp); *inp = NULL; return err);

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

        stdl_log_msg(0, "tda = %s\n", inp->ctx_tda ? "true" : "false");
        stdl_log_msg(0, "gammaJ = %f\ngammaK = %f\nax = %f\n", inp->ctx_gammaJ, inp->ctx_gammaK, inp->ctx_ax);
        stdl_log_msg(0, "ethr = %f # au\ne2thr = %e # au\n", inp->ctx_ethr, inp->ctx_e2thr);

        stdl_log_msg(0, "output = \"%s\"\n", inp->ctx_output);
    } else {
        stdl_log_msg(0, "# A' and B' are obtained from %s\n", inp->ctx_source);
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

int stdl_user_input_make_context(stdl_user_input* inp, stdl_context **ctx_ptr) {
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
        STDL_ERROR_CODE_HANDLE(err, stdl_basis_delete(bs); stdl_wavefunction_delete(wf); return err);

        // select and build A' and B'
        if(inp->ctx_method == STDL_METHOD_MONOPOLE)
            err = stdl_context_select_csfs_monopole(*ctx_ptr, !inp->ctx_tda);
        else if(inp->ctx_method == STDL_METHOD_MONOPOLE_DIRECT)
            err = stdl_context_select_csfs_monopole_direct(*ctx_ptr, !inp->ctx_tda);

        STDL_ERROR_CODE_HANDLE(err, return err);
    }

    if(strcmp(inp->ctx_source, inp->ctx_output) == 0)
        STDL_WARN("`context.source` and `contex.output` are the same, so the content of `%s` will be replaced", inp->ctx_output);

    hid_t file_id = H5Fcreate(inp->ctx_output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    STDL_ERROR_HANDLE_AND_REPORT(file_id == H5I_INVALID_HID, return STDL_ERR_OPEN, "cannot open %s", inp->ctx_output)

    err = stdl_context_dump_h5(*ctx_ptr, file_id);
    H5Fclose(file_id);

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

int stdl_user_input_prepare_responses(stdl_user_input* inp, stdl_context * ctx) {
    assert(inp != NULL && ctx != NULL);

    int err;

    stdl_log_msg(1, "+ ");
    stdl_log_msg(0, "Prepare responses >");
    stdl_log_msg(1, "\n  | Count requests ");

    // count the number of operators, LRV requests, amplitudes, and freqs.
    inp->res_nops = 0;
    inp->res_nlrvreq = 0;
    inp->res_nexci = 0;

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
                inp->res_nops += 1;

            if(!islrvs[op] && req->res_order < req->resp_order)
                inp->res_nlrvreq += 1;

            operators[op] = islrvs[op] = 1;
        }

        size_t nw = req->resp_order - req->res_order;
        for (size_t iw = 0; iw < nw; ++iw) {// TODO: in fact, a quadratic response requires 3 frequencies, not 2 :o
            stdl_operator op = req->ops[iw + 1];
            struct _w_list* elm = NULL;
            err = _w_list_new(req->w[iw], &elm);
            STDL_ERROR_CODE_HANDLE(err, return err);

            if(lrvs_w[op] == NULL) {
                lrvs_w[op] = elm;
            } else {
                struct _w_list* last = lrvs_w[op];
                int already_in = 0;
                while (last->next != NULL && !already_in) {
                    if(stdl_float_equals(req->w[iw], last->w, 1e-6f))
                        already_in = 1;

                    last = last->next;
                }

                if(!already_in)
                    last->next = elm;
            }
        }

        // check out if it requires amplitudes
        if(req->res_order > 0) {
            if(req->nroot < 0)
                inp->res_nexci = ctx->ncsfs;
            else if((size_t) req->nroot > inp->res_nexci)
                inp->res_nexci = (size_t) req->nroot;
        }

        req = req->next;
    }

    STDL_ERROR_HANDLE_AND_REPORT(inp->res_nops == 0, return STDL_ERR_INPUT, "No requests found, exiting");

    stdl_log_msg(0, "-");
    stdl_log_msg(1, "\n  | build requests ");

    inp->res_ops = malloc(inp->res_nops * sizeof(stdl_operator));
    inp->res_lrvreqs = malloc(inp->res_nlrvreq * sizeof(stdl_lrv_request*));
    STDL_ERROR_HANDLE_AND_REPORT(inp->res_ops == NULL || inp->res_lrvreqs == NULL, return STDL_ERR_MALLOC, "malloc");

    // find the number of excited states
    int ioffset = 0;
    for (int iop = 0; iop < STDL_OP_COUNT; ++iop) {
        if(operators[iop]) {
            inp->res_ops[ioffset] = iop;
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

            STDL_ERROR_HANDLE_AND_REPORT(nw == 0, return STDL_ERR_INPUT, "LRV but nw=0");

            // create LRV request
            inp->res_lrvreqs[ioffset] = NULL;
            err = stdl_lrv_request_new(iop, nw, (inp->res_lrvreqs) + ioffset);
            STDL_ERROR_CODE_HANDLE(err, return err);

            lrvs[iop] = inp->res_lrvreqs[ioffset];

            // copy frequencies
            nw = 0;
            struct _w_list* curr = lrvs_w[iop];
            struct _w_list* prev = NULL;
            while (curr != NULL) {
                inp->res_lrvreqs[ioffset]->w[nw] = curr->w;

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
    stdl_log_msg(1, "\n  | assign each response to its request");

    req = inp->res_resreqs;
    while (req != NULL) {
        size_t nw = req->resp_order - req->res_order;
        for (size_t iop = 0; iop < nw; ++iop) {
            stdl_lrv_request* lrvreq = lrvs[req->ops[iop + 1]];
            req->requests[iop] = lrvreq;
            for (size_t jw = 0; jw < lrvreq->nw; ++jw) {
                if(req->w[iop] == lrvreq->w[jw]) {
                    req->wpos[iop] = jw;
                    break;
                }
            }
        }

        req = req->next;
    }

    stdl_log_msg(0, "< done\n");


    return STDL_ERR_OK;
}
