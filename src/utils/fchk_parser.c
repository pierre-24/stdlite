#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "stdlite/errors.h"
#include "stdlite/utils/fchk_parser.h"
#include "stdlite.h"
#include "stdlite/basis.h"

int stdl_fchk_parser_get_section_info(stdl_lexer* lx, char** name, char* type, int* is_scalar) {
    assert(lx != NULL && name != NULL && type != NULL && is_scalar != NULL);

    if(lx->current_tk_type != STDL_TK_ALPHA) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected FCHK section to start with ALPHA");
        return STDL_ERR_UTIL_FCHK;
    }

    // read up the 40 next characters (name)
    int err = STDL_ERR_OK;
    char buff[41]; // 40 + '\0'
    int i = 0, last_alnum = 0;

    while (i < 40) {
        if(lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF)  {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "unexpected end in section name");
            err = STDL_ERR_UTIL_FCHK;
        }

        if(err == STDL_ERR_OK) {
            buff[i] = lx->current_tk_value;

            if(isalnum(lx->current_tk_value))
                last_alnum = i;

            err = stdl_lexer_advance(lx, 1);
        }

        STDL_RETURN_ON_ERROR(err);
        i++;
    }

    // read out the next 7 character (only the fourth is relevant)
    i = 0;
    while (i < 7) {
        if (i == 3) {
            if(lx->current_tk_type != STDL_TK_ALPHA || (lx->current_tk_value != 'I' && lx->current_tk_value != 'R' && lx->current_tk_value != 'C')) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expecting data type to be I/R/C");
                err = STDL_ERR_UTIL_FCHK;
            } else
                *type = lx->current_tk_value;
        }

        if(err == STDL_ERR_OK)
            err = stdl_lexer_advance(lx, 1);

        STDL_RETURN_ON_ERROR(err);
        i++;
    }

    // time to know if it is a scalar or a vector
    if (lx->current_tk_type == STDL_TK_NL || lx->current_tk_type == STDL_TK_EOF) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "line is too short to check for kind (scalar/vector)");
        err = STDL_ERR_UTIL_FCHK;
    } else {
        *is_scalar = lx->current_tk_type != STDL_TK_ALPHA || lx->current_tk_value != 'N';

        // just advance to the value/size
        if(!(*is_scalar)) { // skip '='
            err = stdl_lexer_eat(lx, STDL_TK_ALPHA);

            if(err == STDL_ERR_OK)
                err = stdl_lexer_eat(lx, STDL_TK_EQ);

            if(err != STDL_ERR_OK) {
                stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected EQ after 'N'");
                err = STDL_ERR_UTIL_FCHK;
            }
        }
    }

    if(err == STDL_ERR_OK)
        err = stdl_lexer_skip(lx, isspace);

    STDL_RETURN_ON_ERROR(err);

    // ok, normally we are good. Time to copy the name
    *name = malloc((last_alnum + 2) * sizeof(char));
    if(*name == NULL)
        return STDL_ERR_MALLOC;

    memcpy(*name, buff, last_alnum + 1);
    (*name)[last_alnum + 1] = '\0';

    return STDL_ERR_OK;
}

int _skip_NL(stdl_lexer* lx) {
    /* Skips a NL. Stops right after.
     */
    int err = STDL_ERR_OK;

    if(lx->current_tk_type != STDL_TK_EOF) {
        err = stdl_lexer_eat(lx, STDL_TK_NL);
        if(err != STDL_ERR_OK) {
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected NL");
            err = STDL_ERR_UTIL_FCHK;
        }
    }

    return err;
}


int stdl_fchk_parser_get_scalar_integer(stdl_lexer* lx, long *value) {
    int err = stdl_parser_get_integer(lx, value);
    if(err == STDL_ERR_OK)
        err = _skip_NL(lx);

    return err;
}


int stdl_fchk_parser_get_scalar_number(stdl_lexer* lx, double* value) {
    int err = stdl_parser_get_number(lx, value);
    if(err == STDL_ERR_OK)
        err = _skip_NL(lx);

    return err;
}

int _get_vec_sz(stdl_lexer* lx, size_t* sz) {
    /* Reads the size. Stops after NL.
     */
    assert(lx != NULL && sz != NULL);

    if (lx->current_tk_type != STDL_TK_DIGIT) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected FCHK vector to start with DIGIT");
        return STDL_ERR_UTIL_FCHK;
    }

    long sz_read = -1;
    int err;

    err = stdl_parser_get_integer(lx, &sz_read);
    STDL_RETURN_ON_ERROR(err);

    if (sz_read < 0) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "size of the vector is <0");
        return STDL_ERR_UTIL_FCHK;
    }

    err = _skip_NL(lx);

    if(err == STDL_ERR_OK)
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
        err = stdl_lexer_skip(lx, isspace);

        if(err == STDL_ERR_OK)
            err = stdl_parser_get_integer(lx, &x);

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 6 == 5 && i != sz - 1)
                err = _skip_NL(lx);
        }

        STDL_RETURN_ON_ERROR(err);
        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);
    STDL_RETURN_ON_ERROR(err);

    // ok, we should be at the next section
    return STDL_ERR_OK;
}


int stdl_fchk_parser_get_vector_integers(stdl_lexer* lx, size_t* sz, long **vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_RETURN_ON_ERROR(err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(long));
    if(*vector == NULL)
        return STDL_ERR_MALLOC;

    err = _get_vector_ints(lx, *sz, vector);
    if(err != STDL_ERR_OK) {
        free(*vector);
        return err;
    }

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_integers_immediate(stdl_lexer* lx, size_t sz, long** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    // get vector size
    size_t actual_size;
    int err = _get_vec_sz(lx, &actual_size);
    STDL_RETURN_ON_ERROR(err);

    if(sz != actual_size) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected %d data, got %d", sz, actual_size);
        return STDL_ERR_UTIL_FCHK;
    }

    return _get_vector_ints(lx, sz, vector);
}

int _get_vector_numbers(stdl_lexer* lx, size_t sz, double** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    int err;
    size_t i = 0;
    double x;
    while (i < sz) {
        err = stdl_lexer_skip(lx, isspace);
        if(err == STDL_ERR_OK) {
            err = stdl_parser_get_number(lx, &x);
        }

        if(err == STDL_ERR_OK)  {
            (*vector)[i] = x;

            if(i % 5 == 4 && i != sz - 1)
                err = _skip_NL(lx);
        }

        STDL_RETURN_ON_ERROR(err);
        i++;
    }

    // skip last NL, if any
    err = _skip_NL(lx);
    STDL_RETURN_ON_ERROR(err);

    // now, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_numbers(stdl_lexer* lx, size_t* sz, double** vector) {
    assert(lx != NULL && sz != NULL && vector != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_RETURN_ON_ERROR(err);

    // allocate vector
    *vector = malloc((*sz) * sizeof(double));
    if(*vector == NULL)
        return STDL_ERR_MALLOC;

    err = _get_vector_numbers(lx, *sz, vector);
    if(err != STDL_ERR_OK) {
        free(*vector);
        return err;
    }

    // ok, we should be at the next section
    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_numbers_immediate(stdl_lexer* lx, size_t sz, double** vector) {
    assert(lx != NULL && sz > 0 && vector != NULL && *vector != NULL);

    // get vector size
    size_t actual_size;
    int err = _get_vec_sz(lx, &actual_size);
    STDL_RETURN_ON_ERROR(err);

    if(sz != actual_size) {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "expected %d data, got %d", sz, actual_size);
        return STDL_ERR_UTIL_FCHK;
    }

    return _get_vector_numbers(lx, sz, vector);
}

int _get_vector_string(stdl_lexer* lx, size_t sz, char **out) {
    assert(lx != NULL && sz > 0 && out != NULL && *out != NULL);
    int err = STDL_ERR_OK;

    size_t i = 0, j;
    while (i < sz) {
        // read the 12 characters of a pack
        j = 0;
        while (j < 12 && err == STDL_ERR_OK) {
            (*out)[i * 12 + j] = lx->current_tk_value;
            err = stdl_lexer_advance(lx, 1);
            j++;
        }

        if(err == STDL_ERR_OK && i % 5 == 4 && i != sz - 1)
            err = _skip_NL(lx);

        STDL_RETURN_ON_ERROR(err);
        i++;
    }

    // skip last NL if any
    err = _skip_NL(lx);
    STDL_RETURN_ON_ERROR(err);

    // safeguard
    (*out)[(sz) * 12] = '\0';

    return STDL_ERR_OK;
}

int stdl_fchk_parser_get_vector_string(stdl_lexer* lx, size_t* sz, char **out) {
    assert(lx != NULL && sz != NULL && out != NULL);

    // get vector size
    int err = _get_vec_sz(lx, sz);
    STDL_RETURN_ON_ERROR(err);

    // allocate string
    *out = malloc(((*sz) * 12 + 1) * sizeof(char));
    if(*out == NULL)
        return STDL_ERR_MALLOC;

    err = _get_vector_string(lx, *sz, out);
    if(err != STDL_ERR_OK) {
        free(*out);
        return err;
    }

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
        STDL_RETURN_ON_ERROR(err);

        switch (type) {
            case 'R':
            case 'C':
                nl_to_skip = (int) sz / 5 + ((sz % 5 == 0)? 0 : 1);
                break;
            case 'I':
                nl_to_skip = ((int) sz / 6) + ((sz % 6 == 0)? 0 : 1);
                break;
        }
    }

    int i = 0;
    while (i < nl_to_skip) {
        while (lx->current_tk_type != STDL_TK_NL && lx->current_tk_type != STDL_TK_EOF) {
            err = stdl_lexer_advance(lx, 1);
            STDL_RETURN_ON_ERROR(err);
        }

        err = stdl_lexer_eat(lx, STDL_TK_NL);
        STDL_RETURN_ON_ERROR(err);

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
            STDL_RETURN_ON_ERROR(err);
        }

        if(lx->current_tk_type == STDL_TK_NL) {
            err = stdl_lexer_advance(lx, 1);
            STDL_RETURN_ON_ERROR(err);
        } else { // met EOF, so stop there!
            stdl_error_msg_parser(__FILE__, __LINE__, lx, "FCHK is too short!");
            return STDL_ERR_UTIL_FCHK;
        }

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


void _fchk_data_delete(struct _fchk_data_basis* dt) {
    STDL_FREE_IF_USED(dt->shell_types);
    STDL_FREE_IF_USED(dt->prims_per_shell);
    STDL_FREE_IF_USED(dt->bastoatm);

    STDL_FREE_IF_USED(dt->benv);

    free(dt);
}

// create the structure. If it fails, it is an issue with `malloc()`
struct _fchk_data_basis* _fchk_data_new(size_t nbas, size_t nprims) {
    assert(nbas > 0 && nprims > 0 && nbas <= nprims);

    struct _fchk_data_basis* dt = malloc(sizeof(struct _fchk_data_basis));

    if(dt != NULL) {
        dt->nbas = nbas;
        dt->nprim = nprims;

        // ints
        dt->shell_types = NULL;
        dt->prims_per_shell = NULL;

        // doubles
        dt->benv = NULL;

        // ints
        dt->shell_types = malloc(nbas * sizeof(long));
        dt->prims_per_shell = malloc(nbas * sizeof(long));
        dt->bastoatm = malloc(nbas * sizeof(long));

        if(dt->shell_types == NULL || dt->prims_per_shell == NULL || dt->bastoatm == NULL) {
            _fchk_data_delete(dt);
            return NULL;
        }

        // doubles
        dt->benv = malloc(3 * nprims * sizeof(double));

        if(dt->benv == NULL) {
            _fchk_data_delete(dt);
            return NULL;
        }
    }

    return dt;
}

// make the basis set
int _make_basis_set(stdl_basis** bs_ptr, stdl_wavefunction* wf, struct _fchk_data_basis* dt) {

    int nbas = (int) dt->nbas, extra_coefs = 0, fct_type = 0 /* 1 = cartesian, -1 = spherical */;

    for(int i=0; i < (int) dt->nbas; i++) {
        if(dt->shell_types[i] == -1) { // count extra p function due to sp
            nbas += 1;
            extra_coefs += (int) dt->prims_per_shell[i];
        }

        if(dt->shell_types[i] < -1) {
            if(fct_type == 1) {
                stdl_error_msg(__FILE__, __LINE__, "mixing cartesian and spherical basis functions is not supported");
                return STDL_ERR_UTIL_FCHK;
            }
            fct_type = -1;
        } else if(dt->shell_types[i] > 0) {
            if(fct_type == -1) {
                stdl_error_msg(__FILE__, __LINE__, "mixing cartesian and spherical basis functions is not supported");
                return STDL_ERR_UTIL_FCHK;
            }
            fct_type = 1;
        }
    }

    if(fct_type == 0)
        fct_type = 1;

    size_t env_size = wf->natm * 3 /* coordinates */ + 2 * dt->nprim /* exp + contractions */ + extra_coefs /* sp functions */;
    stdl_debug_msg(__FILE__, __LINE__, "%d atoms and %d basis functions (including sp→s,p) = %ld bytes of env", wf->natm, nbas, env_size);

    STDL_RETURN_ON_ERROR(stdl_basis_new(bs_ptr, (int) wf->natm, nbas, env_size, fct_type == -1));

    // atoms
    for(size_t i=0; i < wf->natm; i++) {
        (*bs_ptr)->atm[i * 6 + 0] = (int) wf->atm[i * 4 + 0];
        (*bs_ptr)->atm[i * 6 + 1] = (int) (i * 3);

        (*bs_ptr)->env[i * 3 + 0] = wf->atm[i * 4 + 1];
        (*bs_ptr)->env[i * 3 + 1] = wf->atm[i * 4 + 2];
        (*bs_ptr)->env[i * 3 + 2] = wf->atm[i * 4 + 3];
    }

    // copy exps and coefs
    size_t offset_exps = 3 * wf->natm, offset_coefs = 3 * wf->natm + dt->nprim, offset_coefs_sp = 3 * wf->natm + 2 * dt->nprim;
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
size_t _get_nao(size_t nbas, long* shell_types) {
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
                stdl_error_msg(__FILE__, __LINE__, "encountered shell type=%d, but it was not taken into account. This will cause issues.", shell_types[i]);
                return 0;
        }
    }

    return total;
}

int stdl_fchk_parser_extract(stdl_wavefunction **wf_ptr, stdl_basis **bs_ptr, stdl_lexer *lx) {
    assert(wf_ptr != NULL && bs_ptr != NULL && lx != NULL);

    // useful variables
    int error = STDL_ERR_OK;
    char* name = NULL;
    char type;
    int is_scalar, finished = 0;

    // extracted from file
    long an_integer;
    double* a_vector = NULL;
    size_t a_size, nelec = 0, natm = 0, nao = 0, nmo = 0, nbas = 0, nprim = 0;

    struct _fchk_data_basis* dt = NULL;
    double* atm = NULL; // [natm*4]

    // read the file and extract what is required
    while (lx->current_tk_type != STDL_TK_EOF && error == STDL_ERR_OK  && !finished) {
        error = stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar);

        if(error == STDL_ERR_OK) {
            stdl_debug_msg(__FILE__, __LINE__, "FCHK: section `%s`", name);

            if(strcmp("Number of electrons", name) == 0) {
                error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
                if(error == STDL_ERR_OK)
                    nelec = (size_t) an_integer;
            } /* --- GEOMETRY: --- */
            else if(strcmp("Nuclear charges", name) == 0) {
                error = stdl_fchk_parser_get_vector_numbers(lx, &natm, &a_vector);
                if(error == STDL_ERR_OK) {
                    atm = malloc( natm * 4 * sizeof(double));
                    if(atm != NULL) {
                        for(size_t i = 0; i < natm; i++)
                            atm[i * 4 + 0] = a_vector[i];
                    } else
                        error = STDL_ERR_MALLOC;

                    free(a_vector);
                }
            } else if(strcmp("Current cartesian coordinates", name) == 0) {
                /* Note: assume that "Nuclear charges" was read before. */
                error = stdl_fchk_parser_get_vector_numbers(lx, &a_size, &a_vector);
                if(error == STDL_ERR_OK) {
                    for(size_t i = 0; i < natm; i++) {
                        atm[i * 4 + 1] = a_vector[i * 3 + 0];
                        atm[i * 4 + 2] = a_vector[i * 3 + 1];
                        atm[i * 4 + 3] = a_vector[i * 3 + 2];
                    }

                    free(a_vector);
                }
            } /* --- BASIS SET: --- */
            else if(strcmp("Number of contracted shells", name) == 0) {
                error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
                if(error == STDL_ERR_OK)
                    nbas = (size_t) an_integer;
            } else if(strcmp("Number of primitive shells", name) == 0) {
                /* Note: assume that "Number of contracted shells" was read before. */
                error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
                if(error == STDL_ERR_OK) {
                    nprim = (size_t) an_integer;
                    dt = _fchk_data_new(nbas, nprim);
                    if(dt == NULL)
                        error = STDL_ERR_MALLOC;
                }
            } else if(strcmp("Shell types", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->shell_types));
                if(error == STDL_ERR_OK) {
                    nao = _get_nao(nbas, dt->shell_types);
                    if(nao != nmo)
                        stdl_warning_msg(__FILE__, __LINE__, "number of MO (%d) and AO (%d) does not match", nmo, nao);
                }
            } else if(strcmp("Number of primitives per shell", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->prims_per_shell));
            } else if(strcmp("Shell to atom map", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                error = stdl_fchk_parser_get_vector_integers_immediate(lx, nbas, &(dt->bastoatm));
            } else if(strcmp("Primitive exponents", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &(dt->benv));
            } else if(strcmp("Contraction coefficients", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                double* ptr = &(dt->benv[nprim]);
                error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &ptr);
            } else if(strcmp("P(S=P) Contraction coefficients", name) == 0) {
                /* Note: assume that "Number of primitive shells" was read before. */
                double* ptr = &(dt->benv[2 * nprim]);
                error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nprim, &ptr);
            } /* --- nmo, C, e --- */
            else if(strcmp("Number of independent functions", name) == 0) {
                error = stdl_fchk_parser_get_scalar_integer(lx, &an_integer);
                if(error == STDL_ERR_OK)
                    nmo = (size_t) an_integer;
            } else if(strcmp("Alpha Orbital Energies", name) == 0) {
                /* Note: assume that "Number of independent functions" was read before. */
                error = stdl_wavefunction_new(wf_ptr, natm, nelec, nao, nmo);
                if(error == STDL_ERR_OK)
                    error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nmo, &((*wf_ptr)->e));
                else
                    error = STDL_ERR_MALLOC;
            } else if(strcmp("Alpha MO coefficients", name) == 0) {
                /* Note: assume that "Number of independent functions" was read before. */
                error = stdl_fchk_parser_get_vector_numbers_immediate(lx, nmo * nao, &((*wf_ptr)->C));

                // Normally, there is nothing else to read!
                finished = 1;
            } /* --- OTHERS: --- */
            else // not interesting, skip section
                stdl_fchk_parser_skip_section(lx, type, is_scalar);

            free(name);
        }
    }

   // at that point, we should have read everything
    if(finished) {
        // copy geometry where it belongs
        memcpy((*wf_ptr)->atm, atm, 4 * natm * sizeof(double));
        free(atm);
        atm = NULL;

        // create basis set
        error = _make_basis_set(bs_ptr, *wf_ptr, dt);
    } else {
        stdl_error_msg_parser(__FILE__, __LINE__, lx, "FCHK was missing certain sections (dt=0x%x, wf_ptr=0x%x)", dt, wf_ptr);
        error = STDL_ERR_READ;
    }

    // clean up stuffs
    STDL_FREE_IF_USED(atm);
    if(dt != NULL)
        _fchk_data_delete(dt);
    if(error != STDL_ERR_OK && *wf_ptr != NULL)
        stdl_wavefunction_delete(*wf_ptr);

    return error;
}
