#include <assert.h>
#include <ctype.h>
#include <string.h>

#include "stdlite/utils/molden_parser.h"
#include "stdlite/logging.h"
#include "stdlite/helpers.h"

int stdl_molden_parser_read_section_title(stdl_lexer* lx, char** title) {

    assert(lx != NULL && title != NULL);

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != '[',
            return STDL_ERR_UTIL_MOLDEN,
            "expected MOLDEN section to start with `[`"
    );

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_pos_in_line != 1,
            return STDL_ERR_UTIL_MOLDEN,
            "expected MOLDEN section to start with `[` as the first character of a line"
    );

    stdl_lexer_advance(lx, 1);

    // read until `]`
    int err;
    char* str;
    int sz = 0, fac = 0;

    err = stdl_grow_string(&str, sz, &fac);
    STDL_ERROR_CODE_HANDLE(err, return err);

    while (lx->current_tk_type != STDL_TK_EOF  && lx->current_tk_type != STDL_TK_NL && lx->current_tk_value != ']') {
        err = stdl_parser_store_value_and_grow_string(lx, &str, &sz, &fac);
        STDL_ERROR_CODE_HANDLE(err, free(str); return err);
    }

    STDL_LEXER_ERROR_HAR(lx, lx->current_tk_value != ']', free(str); return STDL_ERR_UTIL_MOLDEN, "expected MOLDEN section to end with `]`");

    str[sz] = '\0';
    *title = str;

    return STDL_ERR_OK;
}

int _pred_not_begining_of_section(int c) {
    return c != '[';
}

int stdl_molden_parser_skip_section(stdl_lexer* lx) {
    assert(lx != NULL);

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to continue section"
    );

    return stdl_lexer_skip(lx, _pred_not_begining_of_section);
}

// just a linked-chain to temporary store atoms
struct _atom_info_ {
    double atm[4];
    struct _atom_info_* next;
};

void _atom_info_delete(struct _atom_info_* dt) {
    if(dt != NULL) {
        if(dt->next != NULL)
            _atom_info_delete(dt->next);

        free(dt);
    }
}

int _pred_notnl(int c) {
    return c != '\n';
}

int stdl_molden_parser_read_atoms_section(stdl_lexer* lx, size_t* natm, double** atm) {
    assert(lx != NULL && natm != NULL && atm != NULL);

    STDL_DEBUG("Read [Atoms] section");

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to continue section"
    );

    int err;

    err = stdl_lexer_eat(lx, STDL_TK_CHAR);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // Get unit
    int is_bohr = 0;
    err = stdl_lexer_skip(lx, isblank);
    STDL_ERROR_CODE_HANDLE(err, return err);

    char* unit = NULL;
    err = stdl_parser_get_literal(lx, isalpha, &unit);
    STDL_ERROR_CODE_HANDLE(err, return err);

    if(strcmp(unit, "AU") == 0)
        is_bohr = 1;

    free(unit);

    // go to next line
    err = stdl_lexer_skip_whitespace_and_nl(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // read atoms
    *natm = 0;
    struct _atom_info_* first = NULL;
    struct _atom_info_* current = NULL;
    while (lx->current_tk_type != STDL_TK_EOF && lx->current_tk_value != '[') {

        char* tmp;
        struct _atom_info_* inf = malloc(sizeof(struct _atom_info_));
        STDL_ERROR_HANDLE_AND_REPORT(inf == NULL, _atom_info_delete(first); return STDL_ERR_MALLOC, "malloc");

        inf->next = NULL;

        err = stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, _atom_info_delete(inf); _atom_info_delete(first); return err);
        free(tmp);

        long n;
        err = stdl_lexer_skip(lx, isblank)
                || stdl_parser_get_integer(lx, &n);
        STDL_ERROR_CODE_HANDLE(err, _atom_info_delete(inf); _atom_info_delete(first); return err);
        STDL_LEXER_ERROR_HAR(lx, (size_t) n != *natm + 1, _atom_info_delete(inf); _atom_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "expected atoms that are in consecutive order");

        err = stdl_lexer_skip(lx, isblank)
                || stdl_parser_get_integer(lx, &n);
        STDL_ERROR_CODE_HANDLE(err, _atom_info_delete(inf); _atom_info_delete(first); return err);
        inf->atm[0] = (double ) n;

        double d;
        for (int i = 0; i < 3; ++i) {
            err = stdl_lexer_skip(lx, isblank)
                    || stdl_parser_get_number(lx, &d);
            STDL_ERROR_CODE_HANDLE(err, _atom_info_delete(inf); _atom_info_delete(first); return err);
            inf->atm[i + 1] = d * (is_bohr ? 1: 1. / STDL_CONST_AU_TO_ANG);
        }

        err = stdl_lexer_skip_whitespace_and_nl(lx);
        STDL_ERROR_CODE_HANDLE(err, _atom_info_delete(inf); _atom_info_delete(first); return err);

        if(first == NULL) {
            first = inf;
            current = first;
        } else {
            current->next = inf;
            current = current->next;
        }

        *natm += 1;
    }

    STDL_DEBUG("Got %ld atoms", *natm);

    // actually create the `atm` array, now that we know its size
    *atm = malloc(*natm * 4 * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(*atm == NULL, return STDL_ERR_MALLOC, "malloc");

    size_t iatm = 0;
    current = first;
    struct _atom_info_* prev;

    while (current != NULL) {
        memcpy(*atm + iatm * 4, current->atm, 4 * sizeof (double ));
        prev = current;
        current = current->next;
        free(prev); // free linked list while we're at it

        iatm += 1;
    }

    return STDL_ERR_OK;
}

// linked chained list to store basis function
struct _basis_info_ {
    size_t iatm;
    size_t nprim;
    int btype;

    double* env;

    struct _basis_info_* next;
};

void _basis_info_delete(struct _basis_info_* inf) {
    if(inf != NULL) {
        if(inf->next != NULL)
            _basis_info_delete(inf->next);

        STDL_FREE_ALL(inf->env, inf);
    }
}

int _basis_info_new(size_t iatm, char* btype, size_t nprim, struct _basis_info_** inf_ptr) {
    assert(iatm > 0 && btype != NULL && nprim > 0 && inf_ptr != NULL);

    // find basis set type
    int angmom = -42;
    if(strcmp(btype, "s") == 0) {
        angmom = 0;
    } else if(strcmp(btype, "sp") == 0) {
        angmom = -1;
    } else if(strcmp(btype, "p") == 0) {
        angmom = 1;
    } else if(strcmp(btype, "d") == 0) {
        angmom = 2;
    } else if(strcmp(btype, "f") == 0) {
        angmom = 3;
    }

    STDL_ERROR_HANDLE_AND_REPORT(angmom == -42, return STDL_ERR_UTIL_MOLDEN, "Basis function of type `%s` is not handled yet in MOLDEN", btype);

    *inf_ptr = malloc(sizeof(struct _basis_info_));
    STDL_ERROR_HANDLE_AND_REPORT(*inf_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*inf_ptr)->iatm = iatm;
    (*inf_ptr)->nprim = nprim;
    (*inf_ptr)->btype = angmom;
    (*inf_ptr)->next = NULL;
    (*inf_ptr)->env = malloc(3 * nprim * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT((*inf_ptr)->env == NULL, _basis_info_delete(*inf_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_molden_parser_read_gto_section(stdl_lexer* lx, size_t natm, stdl_basis_data** dt_ptr) {
    assert(lx != NULL && natm > 0 && dt_ptr != NULL);

    STDL_DEBUG("Read [GTO] section");

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to continue section"
    );

    int err;
    size_t nbas = 0, nprim = 0;

    err = stdl_lexer_eat(lx, STDL_TK_CHAR) || stdl_lexer_skip_whitespace_and_nl(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    struct _basis_info_* first = NULL;
    struct _basis_info_* current;

    while (lx->current_tk_type != STDL_TK_EOF && lx->current_tk_value != '[') {
        long iatm;
        err = stdl_parser_get_integer(lx, &iatm);
        STDL_ERROR_CODE_HANDLE(err, return err);

        STDL_DEBUG("Read basis functions for iatm=%d", iatm);

        err = stdl_lexer_skip(lx, _pred_notnl) /* the zero at the end is not mandatory */
           || stdl_lexer_eat(lx, STDL_TK_NL)
           || stdl_lexer_skip(lx, isblank);
        STDL_ERROR_CODE_HANDLE(err, return err);

        // read basis
        while (lx->current_tk_type == STDL_TK_ALPHA) {
            char* btype;
            long bas_nprim;

            err = stdl_parser_get_literal(lx, isalpha, &btype)
                    || stdl_lexer_skip(lx, isblank)
                    || stdl_parser_get_integer(lx, &bas_nprim);
            STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(first); return err);

            STDL_DEBUG("basis function `%s` (with %ld primitives)", btype, bas_nprim);

            nprim += bas_nprim;

            struct _basis_info_* inf = NULL;
            err = _basis_info_new(iatm, btype, bas_nprim, &inf);
            STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(first); return err);
            free(btype);

            // the 1.0 at the end is not mandatory
            err = stdl_lexer_skip(lx, _pred_notnl);
            STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(inf); _basis_info_delete(first); return err);

            for (size_t iprim = 0; iprim < inf->nprim; ++iprim) {
                // exponent and coef
                err = stdl_lexer_skip(lx, isblank)
                        || stdl_lexer_eat(lx, STDL_TK_NL)
                        || stdl_lexer_skip(lx, isblank)
                        || stdl_parser_get_number(lx, inf->env + iprim)
                        || stdl_lexer_skip(lx, isblank)
                        || stdl_parser_get_number(lx, inf->env + bas_nprim + iprim);
                STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(inf); _basis_info_delete(first); return err);

                if(inf->btype == -1) { // extra sp coef
                    err = stdl_lexer_skip(lx, isblank)
                          || stdl_parser_get_number(lx, inf->env + 2 * bas_nprim + iprim);
                    STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(inf); _basis_info_delete(first); return err);
                }
            }

            err = stdl_lexer_skip(lx, isblank)
                    || stdl_lexer_eat(lx, STDL_TK_NL)
                    || stdl_lexer_skip(lx, isblank);
            STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(inf); _basis_info_delete(first); return err);

            nbas++;

            if(first == NULL) {
                first = inf;
                current = first;
            } else {
                current->next = inf;
                current = current->next;
            }
        }

        // skip empty line
        err = stdl_lexer_skip(lx, isblank)
              || stdl_lexer_eat(lx, STDL_TK_NL)
              || stdl_lexer_skip(lx, isblank);
        STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(first); return err);
    }

    STDL_DEBUG("found %ld primitives in %ld basis functions", nprim, nbas);

    // finally create data
    err = stdl_basis_data_new(nbas, nprim, dt_ptr);
    STDL_ERROR_CODE_HANDLE(err, _basis_info_delete(first); return err);

    size_t ibas = 0, ioffset = 0;
    current = first;
    struct _basis_info_* prev;
    while (current != NULL) {

        (*dt_ptr)->bas_types[ibas] = current->btype;
        (*dt_ptr)->prims_per_bas[ibas] = (long) current->nprim;
        (*dt_ptr)->bastoatm[ibas] = (long) current->iatm;

        memcpy((*dt_ptr)->benv + ioffset, current->env, current->nprim * sizeof(double ));
        memcpy((*dt_ptr)->benv + nprim + ioffset, current->env + current->nprim, current->nprim * sizeof(double ));
        if(current->btype == -1)
            memcpy((*dt_ptr)->benv + 2 *nprim + ioffset, current->env + 2 * current->nprim, current->nprim * sizeof(double ));

        ioffset += current->nprim;
        ibas++;

        prev = current;
        current = current->next;

        prev->next = NULL;
        _basis_info_delete(prev); // free linked list while we're at it
    }

    return STDL_ERR_OK;
}

// linked chained list to store MO
struct _mo_info_ {
    double e;
    double* c;

    struct _mo_info_* next;
};

void _mo_info_delete(struct _mo_info_* inf) {
    if(inf != NULL) {
        if(inf->next != NULL)
            _mo_info_delete(inf->next);

        STDL_FREE_ALL(inf->c, inf);
    }
}

int _mo_info_new(size_t nao, struct _mo_info_** inf_ptr) {
    assert(nao > 0 && inf_ptr != NULL);

    *inf_ptr = malloc(sizeof(struct _mo_info_));
    STDL_ERROR_HANDLE_AND_REPORT(*inf_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*inf_ptr)->e = .0;
    (*inf_ptr)->next = NULL;
    (*inf_ptr)->c = calloc(nao, sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT((*inf_ptr)->c == NULL, _mo_info_delete(*inf_ptr); return STDL_ERR_MALLOC, "malloc");

    return STDL_ERR_OK;
}

int stdl_molden_parser_read_mo_section(stdl_lexer *lx, size_t nao, size_t *nmo, size_t *nocc, double **e, double **C) {
    assert(lx != NULL && nao > 0 && nmo != NULL && e != NULL && C != NULL);

    STDL_DEBUG("Read [MO] section");

    STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_value != ']',
            return STDL_ERR_UTIL_MOLDEN,
            "expected `]` to continue section"
    );

    int err;

    err = stdl_lexer_eat(lx, STDL_TK_CHAR) || stdl_lexer_skip_whitespace_and_nl(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // read each MO
    *nmo = 0;
    *nocc = 0;

    struct _mo_info_* first = NULL;
    struct _mo_info_* current;
    while (lx->current_tk_type != STDL_TK_EOF && lx->current_tk_value != '[') {
        struct _mo_info_* inf = NULL;
        err = _mo_info_new(nao, &inf);
        STDL_ERROR_CODE_HANDLE(err, _mo_info_delete(first); return err);

        // Sym=occ
        char* tmp;
        err = stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);
        STDL_LEXER_ERROR_HAR(lx, strcmp("Sym", tmp) != 0,  free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "expected Sym");
        free(tmp);

        // Ene=occ
        err = stdl_lexer_skip(lx, _pred_notnl)
           || stdl_lexer_eat(lx, STDL_TK_NL)
           || stdl_lexer_skip(lx, isblank)
           || stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);
        STDL_LEXER_ERROR_HAR(lx, strcmp("Ene", tmp) != 0,  free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "expected Ene");
        free(tmp);

        err = stdl_lexer_eat(lx, STDL_TK_EQ)
                || stdl_lexer_skip(lx, isblank)
                || stdl_parser_get_number(lx, &(inf->e));
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

        // Spin=occ
        err = stdl_lexer_skip(lx, _pred_notnl)
              || stdl_lexer_eat(lx, STDL_TK_NL)
              || stdl_lexer_skip(lx, isblank)
              || stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);
        STDL_LEXER_ERROR_HAR(lx, strcmp("Spin", tmp) != 0,  free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "expected Spin");
        free(tmp);

        // Occ=occ
        err = stdl_lexer_skip(lx, _pred_notnl)
              || stdl_lexer_eat(lx, STDL_TK_NL)
              || stdl_lexer_skip(lx, isblank)
              || stdl_parser_get_literal(lx, isalpha, &tmp);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);
        STDL_LEXER_ERROR_HAR(lx, strcmp("Occup", tmp) != 0,  free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "expected Occup");
        free(tmp);

        double occ;
        err = stdl_lexer_eat(lx, STDL_TK_EQ)
              || stdl_lexer_skip(lx, isblank)
              || stdl_parser_get_number(lx, &occ);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

        if(occ > .01)
            *nocc += 1;

        // skip to coefficients (finally!)
        err = stdl_lexer_skip(lx, isblank)
              || stdl_lexer_eat(lx, STDL_TK_NL)
              || stdl_lexer_skip(lx, isblank);
        STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

        size_t icoef = 0;
        while(lx->current_tk_type == STDL_TK_DIGIT) {
            // read coef
            long n;
            err = stdl_parser_get_integer(lx, &n);
            STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

            err = stdl_lexer_skip(lx, isblank)
                  || stdl_parser_get_number(lx, inf->c + n - 1 /* index in 1-based */);
            STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

            // skip to next coef
            err = stdl_lexer_skip(lx, isblank)
                  || stdl_lexer_eat(lx, STDL_TK_NL)
                  || stdl_lexer_skip(lx, isblank);
            STDL_ERROR_CODE_HANDLE(err, free(tmp); _mo_info_delete(inf); _mo_info_delete(first); return err);

            icoef++;
        }

        STDL_DEBUG("read MO %ld (%ld coefficients set)", *nmo + 1, icoef);

        (*nmo)++;

        if(first == NULL) {
            first = inf;
            current = first;
        } else {
            current->next = inf;
            current = current->next;
        }
    }

    STDL_ERROR_HANDLE_AND_REPORT(*nmo > nao, _mo_info_delete(first); return STDL_ERR_UTIL_MOLDEN, "found %ld MOs, which is larger than the number of AO (%ld)", *nmo, nao);
    STDL_DEBUG("Found %ld MOs", *nmo);

    *e = malloc(*nmo * sizeof(double ));
    *C = malloc(*nmo * nao * sizeof(double ));
    STDL_ERROR_HANDLE_AND_REPORT(*e == NULL || *C == NULL, _mo_info_delete(first); return STDL_ERR_MALLOC, "malloc");

    current = first;
    struct _mo_info_* prev;
    size_t imo = 0;
    while (current != NULL) {
        (*e)[imo] = current->e;
        memcpy((*C) + imo * nao, current->c, nao * sizeof(double ));
        imo++;

        prev = current;
        current = current->next;

        prev->next = NULL;
        _mo_info_delete(prev); // free linked list while we're at it
    }

    return STDL_ERR_OK;
}


int stdl_molden_parser_extract(stdl_lexer* lx, stdl_wavefunction** wf_ptr, stdl_basis** bs_ptr) {
    assert(lx != NULL && wf_ptr != NULL && bs_ptr != NULL);

    stdl_log_msg(0, "Extract wavefunction and basis set from MOLDEN >");

    char* title = NULL;
    int err;

    err = stdl_molden_parser_read_section_title(lx, &title);
    STDL_ERROR_CODE_HANDLE(err, return err);

    STDL_LEXER_ERROR_HAR(lx, strcmp(title, "Molden Format") != 0, free(title); return STDL_ERR_UTIL_MOLDEN, "expected [Molden Format] on the first line of MOLDEN file");
    free(title);

    err = stdl_molden_parser_skip_section(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    size_t natm, nao = 0, nmo, nocc;
    double* atm = NULL, *e = NULL, *C = NULL;
    stdl_basis_data* dt = NULL;

    int all_read = 0;
    while (lx->current_tk_type != STDL_TK_EOF && err == STDL_ERR_OK) {
        err = stdl_molden_parser_read_section_title(lx, &title);
        STDL_ERROR_CODE_HANDLE(err, goto _end);

        if(strcmp(title, "Atoms") == 0) {
            err = stdl_molden_parser_read_atoms_section(lx, &natm, &atm);
        } else if(strcmp(title, "GTO") == 0) {
            err = stdl_molden_parser_read_gto_section(lx, natm, &dt) || stdl_basis_data_count_nao(dt, &nao);
        } else if(strcmp(title, "MO") == 0) {
            err = stdl_molden_parser_read_mo_section(lx, nao, &nmo, &nocc, &e, &C);
        } else {
            STDL_DEBUG("Skip section [%s]", title);
            err = stdl_molden_parser_skip_section(lx);
        }

        free(title);
    }

    all_read = atm != NULL && dt != NULL && C != NULL && e != NULL;
    STDL_ERROR_HANDLE_AND_REPORT(!all_read, err = STDL_ERR_UTIL_MOLDEN; goto _end, "MOLDEN file was missing certain sections");

    stdl_log_msg(0, "-");

    // create wavefunction
    stdl_wavefunction_new(natm, nocc, nao, nmo, wf_ptr);

    memcpy((*wf_ptr)->e, e, nmo * sizeof(double ));
    memcpy((*wf_ptr)->C, C, nao * nmo * sizeof(double ));

    stdl_log_msg(0, "-");

    // create basis
    err = stdl_basis_data_to_basis(dt, natm, atm, bs_ptr);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    stdl_basis_data_delete(dt);
    dt = NULL;

    stdl_log_msg(0, "-");

    // map each AO to its atom
    int si, center, shift = 0, nprim = 0;
    for (int ibas = 0; ibas < (*bs_ptr)->nbas; ++ibas) {
        center = (*bs_ptr)->bas[ibas * 8 + 0];
        nprim += (*bs_ptr)->bas[ibas * 8 + 2];

        if((*bs_ptr)->use_spherical)
            si = CINTcgto_spheric(ibas, (*bs_ptr)->bas);
        else
            si = CINTcgtos_cart(ibas, (*bs_ptr)->bas);

        for(int j = 0; j < si; j++) {
            (*wf_ptr)->aotoatm[shift + j] = center;
        }

        shift += si;
    }

    stdl_log_msg(0, "-");

    // TODO: reorder AOs

    // create the S matrix
    err = stdl_basis_dsp_ovlp((*bs_ptr), (*wf_ptr)->S);
    STDL_ERROR_CODE_HANDLE(err, goto _end);

    stdl_log_msg(0, "< done\n");
    stdl_log_msg(0, "Got %d atoms, %d AOs (%d primitives in %d basis functions), and %d MOs\n", natm, nao, nprim, (*bs_ptr)->nbas, nmo);

    _end:
    if(dt != NULL)
        stdl_basis_data_delete(dt);

    STDL_FREE_ALL(C, e, atm);

    return err;
}
