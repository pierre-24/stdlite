#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

#include "stdlite/logging.h"
#include "stdlite/utils/fchk_parser.h"
#include "stdlite.h"
#include "stdlite/basis.h"

int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar) {
    assert(lx != NULL && name != NULL && type != NULL && is_scalar != NULL);

    STDL_LEXER_ERROR_HAR(
        lx,
        lx->current_tk_type != STDL_TK_ALPHA,
        return STDL_ERR_UTIL_FCHK,
        "expected FCHK section to start with ALPHA"
    );

    // read up the 40 next characters (name)
    int err;
    char buff[41]; // 40 + '\0'
    int i = 0, last_alnum = 0;

    while (i < 40) {
        STDL_LEXER_ERROR_HAR(
            lx,
            lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF,
            return STDL_ERR_UTIL_FCHK,
            "unexpected end in section name"
        );

        buff[i] = lx->current_tk_value;

        if(isalnum(lx->current_tk_value))
            last_alnum = i;

        err = stdl_lexer_advance(lx, 1);
        STDL_ERROR_CODE_HANDLE(err, return err);

        i++;
    }

    // read out the next 7 character (only the fourth is relevant)
    i = 0;
    while (i < 7) {
        if (i == 3) {
            STDL_LEXER_ERROR_HAR(
                lx,
                lx->current_tk_type != STDL_TK_ALPHA || (lx->current_tk_value != 'I' && lx->current_tk_value != 'R' && lx->current_tk_value != 'C'),
                return STDL_ERR_UTIL_FCHK,
                "expecting data type to be I/R/C"
            );

            *type = lx->current_tk_value;
        }

        err = stdl_lexer_advance(lx, 1);
        STDL_ERROR_CODE_HANDLE(err, return err);
        i++;
    }


    // time to know if it is a scalar or a vector
    STDL_LEXER_ERROR_HAR(
        lx,
        lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF,
        return STDL_ERR_UTIL_FCHK,
        "line is too short to check for kind (scalar/vector)"
    );

    *is_scalar = lx->current_tk_type != STDL_TK_ALPHA || lx->current_tk_value != 'N';

    // just advance to the value/size
    if(!(*is_scalar)) {
        // skip `N`
        stdl_lexer_eat(lx, STDL_TK_ALPHA);

        // skip `=`
        err = stdl_lexer_eat(lx, STDL_TK_EQ);
        STDL_ERROR_CODE_HANDLE(err, return err);
    }

    err = stdl_lexer_skip(lx, isspace);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // ok, normally we are good. Time to copy the name
    *name = malloc((last_alnum + 2) * sizeof(char));
    STDL_ERROR_HANDLE_AND_REPORT(*name == NULL, return STDL_ERR_MALLOC, "malloc");

    memcpy(*name, buff, last_alnum + 1);
    (*name)[last_alnum + 1] = '\0';

    return STDL_ERR_OK;
}

// Skips a NL. Stops right after.
int _skip_NL(stdl_lexer* lx) {
    if(lx->current_tk_type != STDL_TK_EOF) {
        int err = stdl_lexer_eat(lx, STDL_TK_NL);
        STDL_ERROR_CODE_HANDLE(err, return err);
    }

    return STDL_ERR_OK;
}


int stdl_fchk_parser_get_scalar_integer(stdl_lexer* lx, long *value) {
    return stdl_parser_get_integer(lx, value) || _skip_NL(lx);
}


int stdl_fchk_parser_get_scalar_number(stdl_lexer* lx, double* value) {
    return stdl_parser_get_number(lx, value) || _skip_NL(lx);
}

// Reads the size. Stops after NL.
int _get_vec_sz(stdl_lexer* lx, size_t* sz) {
    assert(lx != NULL && sz != NULL);

    STDL_LEXER_ERROR_HAR(
        lx,
        lx->current_tk_type != STDL_TK_DIGIT,
        return STDL_ERR_UTIL_FCHK,
        "expected FCHK vector to start with DIGIT"
    );

    long sz_read = -1;
    int err;

    err = stdl_parser_get_integer(lx, &sz_read);
    STDL_ERROR_CODE_HANDLE(err, return err);

    STDL_LEXER_ERROR_HAR(
        lx,
        sz_read < 0,
        return STDL_ERR_UTIL_FCHK,
        "size of the vector is <0"
    );

    err = _skip_NL(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    *sz = (size_t) sz_read;
    return err;
}

int _get_vector_ints(stdl_lexer* lx, size_t sz, long **vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    // get vector size
    int err;

    size_t i = 0;
    long x;
    while (i < sz) {
        err = stdl_lexer_skip(lx, isspace) || stdl_parser_get_integer(lx, &x);
        STDL_ERROR_CODE_HANDLE(err, return err);

        (*vector)[i] = x;

        if(i % 6 == 5 && i != sz - 1) {
            err = _skip_NL(lx);
            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // ok, we should be at the next section
    return STDL_ERR_OK;
}


int stdl_fchk_parser_get_vector_integers(stdl_lexer* lx, size_t* sz, long **vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(long));
    STDL_ERROR_HANDLE_AND_REPORT(*vector == NULL, return STDL_ERR_MALLOC, "malloc");

    // read
    err = _get_vector_ints(lx, *sz, vector);
    STDL_ERROR_CODE_HANDLE(err, free(*vector); return err);

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_integers_immediate(stdl_lexer* lx, size_t sz, long** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    // get vector size
    size_t actual_size;
    int err = _get_vec_sz(lx, &actual_size);
    STDL_ERROR_CODE_HANDLE(err, return err);

    STDL_LEXER_ERROR_HAR(
        lx,
        sz != actual_size,
        return STDL_ERR_UTIL_FCHK,
        "expected %d data, got %d", sz, actual_size
    );

    return _get_vector_ints(lx, sz, vector);
}

int _get_vector_numbers(stdl_lexer* lx, size_t sz, double** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    int err;
    size_t i = 0;
    double x;
    while (i < sz) {
        err = stdl_lexer_skip(lx, isspace) || stdl_parser_get_number(lx, &x);
        STDL_ERROR_CODE_HANDLE(err, return err);

        (*vector)[i] = x;

        if(i % 5 == 4 && i != sz - 1) {
            err = _skip_NL(lx);
            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // now, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_numbers(stdl_lexer* lx, size_t* sz, double** vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(double));
    STDL_ERROR_HANDLE_AND_REPORT(*vector == NULL, return STDL_ERR_MALLOC, "malloc");

    // read
    err = _get_vector_numbers(lx, *sz, vector);
    STDL_ERROR_CODE_HANDLE(err, free(*vector); return err);

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_numbers_immediate(stdl_lexer* lx, size_t sz, double** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    // get vector size
    size_t actual_size;
    int err = _get_vec_sz(lx, &actual_size);
    STDL_ERROR_CODE_HANDLE(err, return err);

    STDL_LEXER_ERROR_HAR(
        lx,
        sz != actual_size,
        return STDL_ERR_UTIL_FCHK,
        "expected %d data, got %d", sz, actual_size
    );

    return _get_vector_numbers(lx, sz, vector);
}

int _get_vector_string(stdl_lexer* lx, size_t sz, char **out) {
    assert(lx != NULL && sz > 0 && out != NULL && *out != NULL);
    int err;

    size_t i = 0, j;
    while (i < sz) {
        // read the 12 characters of a pack
        j = 0;
        while (j < 12) {
            (*out)[i * 12 + j] = lx->current_tk_value;
            err = stdl_lexer_advance(lx, 1);
            STDL_ERROR_CODE_HANDLE(err, return err);
            j++;
        }

        if(i % 5 == 4 && i != sz - 1) { // ensure NL
            err = _skip_NL(lx);
            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        i++;
    }

    // skip last NL if any
    err = _skip_NL(lx);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // safeguard
    (*out)[(sz) * 12] = '\0';

    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, char **out) {
    assert(lx != NULL && sz != NULL && out != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_ERROR_CODE_HANDLE(err, return err);

    // allocate string
    *out = malloc(((*sz) * 12 + 1) * sizeof(char));
    STDL_ERROR_HANDLE_AND_REPORT(*out == NULL, return STDL_ERR_MALLOC, "malloc");

    err = _get_vector_string(lx, *sz, out);
    STDL_ERROR_CODE_HANDLE(err, free(*out); return err);

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_skip_section(stdl_lexer* lx, char type, int is_scalar) {
    assert(lx != NULL);

    int err, nl_to_skip = 1;
    if(!is_scalar) {
        // fetch vector size
        size_t sz = -1;
        err = _get_vec_sz(lx, &sz);
        STDL_ERROR_CODE_HANDLE(err, return err);

        switch (type) {
            case 'R':
            case 'C':
                nl_to_skip = (int) sz / 5 + ((sz % 5 == 0)? 0 : 1);
                break;
            case 'I':
                nl_to_skip = ((int) sz / 6) + ((sz % 6 == 0)? 0 : 1);
                break;
            default:
                STDL_LEXER_ERROR_HAR(lx, 1, return STDL_ERR_UTIL_FCHK, "unknown label %c", lx->current_tk_value);
        }
    }

    int i = 0;
    while (i < nl_to_skip) {
        while (lx->current_tk_type != STDL_TK_NL && lx->current_tk_type != STDL_TK_EOF) {
            err = stdl_lexer_advance(lx, 1);
            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        err = stdl_lexer_eat(lx, STDL_TK_NL);
        STDL_ERROR_CODE_HANDLE(err, return err);

        i++;
    }

    return STDL_ERR_OK;
}

int stdl_fchk_parser_skip_intro(stdl_lexer* lx) {
    assert(lx != NULL);

    int i = 0, err;
    while (i < 2) {
        while (lx->current_tk_type != STDL_TK_NL && lx->current_tk_type != STDL_TK_EOF) {
            err = stdl_lexer_advance(lx, 1);
            STDL_ERROR_CODE_HANDLE(err, return err);
        }

        STDL_LEXER_ERROR_HAR(lx, lx->current_tk_type != STDL_TK_NL, return STDL_ERR_UTIL_FCHK, "FCHK is too short!");

        err = stdl_lexer_advance(lx, 1);
        STDL_ERROR_CODE_HANDLE(err, return err);

        i++;
    }

    return STDL_ERR_OK;
}


// what is required to build a basis set
struct _fchk_data_basis {
    size_t
        nbas, // number of basis functions
        nprim // number of primitives, nprim >= nbas
        ;
    long
        *shell_types,  // [nbas]
        *prims_per_shell, // [nbas]
        *bastoatm // [nbas], 1-based list
        ;
    double
        *benv // [3 * nprim]
        ;
};


int _fchk_data_delete(struct _fchk_data_basis* dt) {
    STDL_FREE_ALL(dt->shell_types, dt->prims_per_shell, dt->bastoatm, dt->benv,dt);

    return STDL_ERR_OK;
}

// create the structure. If it fails, it is an issue with `malloc()`
int _fchk_data_new(struct _fchk_data_basis **dt_ptr, size_t nbas, size_t nprims) {
    assert(nbas > 0 && nprims > 0 && nbas <= nprims);

    *dt_ptr = malloc(sizeof(struct _fchk_data_basis));
    STDL_ERROR_HANDLE_AND_REPORT(*dt_ptr == NULL, return STDL_ERR_MALLOC, "malloc");

    (*dt_ptr)->nbas = nbas;
    (*dt_ptr)->nprim = nprims;

    // ints
    (*dt_ptr)->shell_types = NULL;
    (*dt_ptr)->prims_per_shell = NULL;

    // doubles
    (*dt_ptr)->benv = NULL;

    // ints
    (*dt_ptr)->shell_types = malloc(nbas * sizeof(long));
    (*dt_ptr)->prims_per_shell = malloc(nbas * sizeof(long));
    (*dt_ptr)->bastoatm = malloc(nbas * sizeof(long));

    STDL_ERROR_HANDLE_AND_REPORT(
        (*dt_ptr)->shell_types == NULL || (*dt_ptr)->prims_per_shell == NULL || (*dt_ptr)->bastoatm == NULL,
        _fchk_data_delete(*dt_ptr); return STDL_ERR_MALLOC,
        "malloc"
    );


    // doubles
    (*dt_ptr)->benv = malloc(3 * nprims * sizeof(double));

    STDL_ERROR_HANDLE_AND_REPORT(
        (*dt_ptr)->benv == NULL,
        _fchk_data_delete(*dt_ptr); return STDL_ERR_MALLOC,
        "malloc"
    );

    return STDL_ERR_OK;
}

// make the basis set, so shuffle things around to convert the format of Gaussian to the one of libcint.
// If this fails, it is most probably from `malloc()`.
int _make_basis_set(stdl_basis** bs_ptr, stdl_wavefunction* wf, struct _fchk_data_basis* dt) {
    assert(bs_ptr != NULL && wf != NULL && dt != NULL);

    STDL_DEBUG("creating the basis set");

    int nbas = (int) dt->nbas, extra_coefs = 0, fct_type = 0 /* 1 = cartesian, -1 = spherical */;

    for(int i=0; i < (int) dt->nbas; i++) {
        if(dt->shell_types[i] == -1) { // count extra p function due to sp
            nbas += 1;
            extra_coefs += (int) dt->prims_per_shell[i];
        }

        if(dt->shell_types[i] < -1) {
            STDL_ERROR_HANDLE_AND_REPORT(fct_type == 1, return STDL_ERR_UTIL_FCHK, "mixing cartesian and spherical basis functions is not supported");
            fct_type = -1;
        } else if(dt->shell_types[i] > 0) {
            STDL_ERROR_HANDLE_AND_REPORT(fct_type == -1, return STDL_ERR_UTIL_FCHK, "mixing cartesian and spherical basis functions is not supported");
            fct_type = 1;
        }
    }

    if(fct_type == 0)
        fct_type = 1;

    size_t env_size = 20 + wf->natm * 3 /* coordinates */ + 2 * dt->nprim /* exp + contractions */ + extra_coefs /* sp functions */;
    stdl_debug_msg(__FILE__, __LINE__, "%d atoms and %d basis functions (including sp→s,p) = %ld bytes of env", wf->natm, nbas, env_size);

    int err = stdl_basis_new((int) wf->natm, nbas, env_size, fct_type == -1, bs_ptr);
    STDL_ERROR_CODE_HANDLE(err, return err);

    size_t base_offset = 20;

    // atoms
    for(size_t i=0; i < wf->natm; i++) {
        (*bs_ptr)->atm[i * 6 + 0] = (int) wf->atm[i * 4 + 0];
        (*bs_ptr)->atm[i * 6 + 1] = (int) (base_offset + i * 3);

        (*bs_ptr)->env[base_offset + i * 3 + 0] = wf->atm[i * 4 + 1];
        (*bs_ptr)->env[base_offset + i * 3 + 1] = wf->atm[i * 4 + 2];
        (*bs_ptr)->env[base_offset + i * 3 + 2] = wf->atm[i * 4 + 3];
    }

    // copy exps and coefs
    size_t offset_exps = base_offset + 3 * wf->natm, offset_coefs = base_offset +3 * wf->natm + dt->nprim, offset_coefs_sp = base_offset +3 * wf->natm + 2 * dt->nprim;
    memcpy(&((*bs_ptr)->env[offset_exps]), dt->benv, 2 * dt->nprim * sizeof(double));

    size_t offset_bas = 0;
    int iprim = 0, ipprim=0;

    // basis
    for(size_t i=0; i < dt->nbas; i++) {
        int angular = (dt->shell_types[i] == -1) ? 0 : abs((int) dt->shell_types[i]);

        (*bs_ptr)->bas[(offset_bas + i) * 8 + 0] = (int) dt->bastoatm[i] - 1; // Gaussian gives a 1-based list
        (*bs_ptr)->bas[(offset_bas + i) * 8 + 1] = angular;
        (*bs_ptr)->bas[(offset_bas + i) * 8 + 2] = (int) dt->prims_per_shell[i];
        (*bs_ptr)->bas[(offset_bas + i) * 8 + 3] = 1;
        (*bs_ptr)->bas[(offset_bas + i) * 8 + 5] = (int) offset_exps + iprim;
        (*bs_ptr)->bas[(offset_bas + i) * 8 + 6] = (int) offset_coefs + iprim;

        // normalize coefs
        for(int j=0; j < (int) dt->prims_per_shell[i]; j++)
            (*bs_ptr)->env[(int) offset_coefs + iprim + j] *= CINTgto_norm(angular, (*bs_ptr)->env[(int) offset_exps + iprim + j]);

        if(dt->shell_types[i] == -1) { // sp
            offset_bas += 1;

            (*bs_ptr)->bas[(offset_bas + i) * 8 + 0] = (int) dt->bastoatm[i] - 1; // Gaussian gives a 1-based list
            (*bs_ptr)->bas[(offset_bas + i) * 8 + 1] = 1;
            (*bs_ptr)->bas[(offset_bas + i) * 8 + 2] = (int) dt->prims_per_shell[i];
            (*bs_ptr)->bas[(offset_bas + i) * 8 + 3] = 1;
            (*bs_ptr)->bas[(offset_bas + i) * 8 + 5] = (int) offset_exps + iprim;
            (*bs_ptr)->bas[(offset_bas + i) * 8 + 6] = (int) offset_coefs_sp + ipprim;

            // copy p-coefs and normalize them
            memcpy(&((*bs_ptr)->env[offset_coefs_sp + ipprim]), &(dt->benv[2 * dt->nprim + iprim]), ((int) dt->prims_per_shell[i]) * sizeof(double));
            for(int j=0; j < (int) dt->prims_per_shell[i]; j++)
                (*bs_ptr)->env[(int) offset_coefs_sp + ipprim + j] *= CINTgto_norm(1, (*bs_ptr)->env[(int) offset_exps + iprim + j]);

            ipprim += (int) dt->prims_per_shell[i];
        }

        iprim += (int) dt->prims_per_shell[i];
    }

    return STDL_ERR_OK;
}

// Get the number of ao, given the type of each basis function.
// According to the Gaussian documentation, `0=s, 1=p, -1=sp, 2=6d, -2=5d, 3=10f, -3=7f` (and so all).
// In most of the case, `nao=nmo`.
size_t _count_nao(size_t nbas, long* shell_types) {
    assert(nbas > 0);

    size_t total = 0;

    for (size_t i = 0; i < nbas; ++i) {
        switch (shell_types[i]) {
            case 0:
                total += 1;
                break;
            // cartesian's
            case 1: // 3p
                total += 3;
                break;
            case 2: // 6d
                total += 6;
                break;
            case 3: // 10f
                total += 10;
                break;
            case 4: // 15g
                total += 15;
                break;
            case 5: // 21h
                total += 21;
                break;
            case 6: // 28i
                total += 28;
                break;
            // spherical's
            case -1: // sp
                total += 4;
                break;
            case -2: // 5d
                total += 5;
                break;
            case -3: // 7f
                total += 7;
                break;
            case -4: // 9g
                total += 9;
                break;
            case -5: // 11h
                total += 11;
                break;
            case -6: // 13i
                total += 13;
                break;
            default:
                STDL_ERROR_HANDLE_AND_REPORT(1, return 0, "encountered shell type=%d, but it was not taken into account. This will cause issues.", shell_types[i]);
        }
    }

    return total;
}

int stdl_fchk_parser_extract(stdl_lexer *lx, stdl_wavefunction **wf_ptr, stdl_basis **bs_ptr) {
    assert(wf_ptr != NULL && bs_ptr != NULL && lx != NULL);

    stdl_log_msg(0, "Extract wavefunction and basis set from FCHK >");

    STDL_DEBUG("reading FCHK file");

    // useful variables
    int error = STDL_ERR_OK;
    char* name = NULL;
    char type;
    int is_scalar, finished = 0;

    // extracted from file
    long an_integer;
    double* a_vector = NULL;
    size_t a_size, nocc = 0, natm = 0, nao = 0, nmo = 0, nbas = 0, nprim = 0;

    struct _fchk_data_basis* dt = NULL;
    double* atm = NULL; // [natm*4]

    // read the file and extract what is required
    while (lx->current_tk_type != STDL_TK_EOF && error == STDL_ERR_OK  && !finished) {
        error = stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar);
        STDL_ERROR_CODE_HANDLE(error, goto _end);

        STDL_DEBUG("- reading section `%s`", name);

        if(strcmp("Number of electrons", name) == 0) {
            error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
            nocc = (size_t) an_integer / 2;
        } /* --- GEOMETRY: --- */
        else if(strcmp("Nuclear charges", name) == 0) {
            error = stdl_fchk_parser_get_vector_numbers(lx, &natm, &a_vector);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            atm = malloc( natm * 4 * sizeof(double));
            STDL_ERROR_HANDLE_AND_REPORT(atm == NULL, error = STDL_ERR_MALLOC; free(name); goto _end, "malloc");

            for(size_t i = 0; i < natm; i++)
                atm[i * 4 + 0] = a_vector[i];

            free(a_vector);
        } else if(strcmp("Current cartesian coordinates", name) == 0) {
            /* Note: assume that "Nuclear charges" was read before. */
            error = stdl_fchk_parser_get_vector_numbers(lx, &a_size, &a_vector);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            for(size_t i = 0; i < natm; i++) {
                atm[i * 4 + 1] = a_vector[i * 3 + 0];
                atm[i * 4 + 2] = a_vector[i * 3 + 1];
                atm[i * 4 + 3] = a_vector[i * 3 + 2];
            }

            free(a_vector);
        } /* --- BASIS SET: --- */
        else if(strcmp("Number of contracted shells", name) == 0) {
            error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            nbas = (size_t) an_integer;
        } else if(strcmp("Number of primitive shells", name) == 0) {
            /* Note: assume that "Number of contracted shells" was read before. */
            error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            nprim = (size_t) an_integer;
            error = _fchk_data_new(&dt, nbas, nprim);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("Shell types", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->shell_types));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            nao = _count_nao(nbas, dt->shell_types);
            if(nao != nmo)
                stdl_warning_msg(__FILE__, __LINE__, "number of MO (%d) and AO (%d) does not match", nmo, nao);
        } else if(strcmp("Number of primitives per shell", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->prims_per_shell));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("Shell to atom map", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->bastoatm));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("Primitive exponents", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &(dt->benv));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("Contraction coefficients", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            double* ptr = &(dt->benv[nprim]);
            error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &ptr);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("P(S=P) Contraction coefficients", name) == 0) {
            /* Note: assume that "Number of primitive shells" was read before. */
            double* ptr = &(dt->benv[2 * nprim]);
            error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &ptr);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } /* --- nmo, C, e --- */
        else if(strcmp("Number of independent functions", name) == 0) {
            error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            nmo = (size_t) an_integer;
        } else if(strcmp("Alpha Orbital Energies", name) == 0) {
            /* Note: assume that "Number of independent functions" was read before. */
            error = stdl_wavefunction_new(natm, nocc, nao, nmo, wf_ptr);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nmo, &((*wf_ptr)->e));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        } else if(strcmp("Alpha MO coefficients", name) == 0) {
            /* Note: assume that "Number of independent functions" was read before. */
            error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nmo * nao, &((*wf_ptr)->C));
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);

            // Normally, there is nothing else to read!
            finished = 1;
        } /* --- OTHERS: --- */
        else { // not interesting, skip section
            error = stdl_fchk_parser_skip_section(lx, type, is_scalar);
            STDL_ERROR_CODE_HANDLE(error, free(name); goto _end);
        }

        free(name);
    }

    stdl_log_msg(0, "-");

    STDL_DEBUG("reading done");

    STDL_ERROR_HANDLE_AND_REPORT(!finished, error = STDL_ERR_UTIL_FCHK; goto _end, "FCHK was missing certain sections (dt=0x%x, wf_ptr=0x%x)", dt, wf_ptr);

    // at that point, we should have read everything
    // copy remaining stuffs
    memcpy((*wf_ptr)->atm, atm, 4 * natm * sizeof(double));
    free(atm);
    atm = NULL;

    stdl_log_msg(0, "-");

    // create the basis set
    error = _make_basis_set(bs_ptr, *wf_ptr, dt);
    STDL_ERROR_CODE_HANDLE(error, goto _end);

    _fchk_data_delete(dt);
    dt = NULL;

    stdl_log_msg(0, "-");

    // map each AO to its atom
    int si, center, shift = 0;
    for (int i = 0; i < (*bs_ptr)->nbas; ++i) {
        center = (*bs_ptr)->bas[i * 8 + 0];
        if((*bs_ptr)->use_spherical)
            si = CINTcgto_spheric(i, (*bs_ptr)->bas);
        else
            si = CINTcgtos_cart(i, (*bs_ptr)->bas);

        for(int j = 0; j < si; j++) {
            (*wf_ptr)->aotoatm[shift + j] = center;
        }

        shift += si;
    }

    stdl_log_msg(0, "-");

    // create the S matrix
    error = stdl_basis_dsp_ovlp((*bs_ptr), (*wf_ptr)->S);

    stdl_log_msg(0, "< done\n");

    // clean up stuffs
    _end:
        STDL_FREE_IF_USED(atm);
        if(dt != NULL)
            _fchk_data_delete(dt);
        if(error != STDL_ERR_OK && *wf_ptr != NULL)
            stdl_wavefunction_delete(*wf_ptr);

    return error;
}
