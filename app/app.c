#include <stdlib.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/utils/molden_parser.h>
#include <string.h>

#include "app.h"

int stdl_app_set_debug_log_level() {
    char* end = NULL;
    char* loglevel = getenv("LOGLEVEL");
    if(loglevel != NULL) {
        long loglevel_d = strtol(loglevel, &end, 10);
        if(loglevel != end)
            stdl_set_log_level((int) loglevel_d);
    }

    char* debuglevel = getenv("DEBUGLEVEL");
    if(debuglevel != NULL) {
        long debuglevel_d = strtol(debuglevel, &end, 10);
        if(debuglevel != end)
            stdl_set_debug_level((int) debuglevel_d);
    }

    return STDL_ERR_OK;
}

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

int stdl_app_context(stdl_user_input* inp, stdl_context **ctx_ptr) {
    stdl_wavefunction* wf = NULL;
    stdl_basis* bs = NULL;
    int err;

    if(inp->ctx_source_type == STDL_SRC_CTX) { // directly read context from a previous calculation
        err = stdl_context_load_h5(inp->ctx_source, ctx_ptr);
        STDL_ERROR_CODE_HANDLE(err, return err);
    } else { // read file for wavefunction & basis, then select normally
        if(inp->ctx_source_type == STDL_SRC_CTX_WB) {
            err = _from_h5(inp->ctx_source, &wf, &bs);
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
        STDL_WARN("source and output are the same, so the content of `%s` is replaced", inp->ctx_source);

    err = stdl_context_dump_h5(*ctx_ptr, inp->ctx_output);

    return err;
}
