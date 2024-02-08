#include <stdlib.h>
#include <string.h>

#include <stdlite/utils/fchk_parser.h>

#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}


void test_parse_section_info_ok() {
    // scalar
    char* scalar = "Multiplicity                               I                1";
    char* title = NULL;
    int is_scalar = -1;
    char type;
    stdl_lexer* lx = NULL;

    rewind(stream);
    fputs(scalar, stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

    TEST_ASSERT_EQUAL_STRING("Multiplicity", title);
    TEST_ASSERT_EQUAL_CHAR('I', type);
    TEST_ASSERT_EQUAL_INT(1, is_scalar);
    TEST_ASSERT_EQUAL_INT(STDL_TK_DIGIT, lx->current_tk_type);

    free(title);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    // vector
    char* vector = "Atomic numbers                             I   N=           3";

    rewind(stream);
    fputs(vector, stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

    TEST_ASSERT_EQUAL_STRING("Atomic numbers", title);
    TEST_ASSERT_EQUAL_CHAR('I', type);
    TEST_ASSERT_EQUAL_INT(0, is_scalar);
    TEST_ASSERT_EQUAL_INT(STDL_TK_DIGIT, lx->current_tk_type);

    free(title);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}

void test_parse_incorrect_section_info_ko() {
    char* incorrect[] = {
            "Multi", // too short
            "Multiplicity                               ", // too short
            "Multiplicity                               I   ", // too short
            " Multiplicity                               I                1", // begin with space
            "Multiplicity                               Q                1", // incorrect type
            "Multiplicity                                I                1", // incorrect position for type
            "Atomic numbers                             I   N            3", // missing "="
    };

    char* title = NULL;
    int is_scalar = -1;
    char type;
    stdl_lexer* lx = NULL;

    for(int i=0; i < 7; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
        STDL_NOK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

        free(title);

        ASSERT_STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_parser_vector_ints_ok() {
    char* fmt =
            "%d\n"
            "           0           1           2           3           4           5\n"
            "           6           7";
    int to_test[] = {6, 8};

    long* values;
    size_t sz;
    stdl_lexer* lx = NULL;

    for(int i = 0; i < 2; i++) {
        int mx = to_test[i];
        rewind(stream);
        fprintf(stream, fmt, mx);
        rewind(stream);

        ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
        ASSERT_STDL_OK(stdl_fchk_parser_get_vector_integers(lx, &sz, &values));

        TEST_ASSERT_EQUAL_INT(mx, sz);
        for (long j = 0; j < mx; ++j) {
            TEST_ASSERT_EQUAL_INT(j, values[j]);
        }

        free(values);

        ASSERT_STDL_OK(stdl_lexer_delete(lx));
    }

}

void test_parser_incorrect_vector_ints_ko() {
    char* incorrect[] = {
            "2", // just a number
            "2\n", // no value
            "2\n1", // missing one value
            "2\n1  x", // value is not integer
            "2\n 1 1x", // should end with \n
    };

    long* values;
    size_t sz;
    stdl_lexer* lx = NULL;

    for(int i = 0; i < 5; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
        STDL_NOK(stdl_fchk_parser_get_vector_integers(lx, &sz, &values));
        ASSERT_STDL_OK(stdl_lexer_delete(lx));
    }

}

void test_parser_vector_numbers_ok() {
    char* sec =
            "9\n"
            "  1.38394898E-17  1.38394898E-17  2.26016022E-01  2.09095447E-16  1.43961729E+00\n"
            " -9.04064083E-01 -4.96113637E-16 -1.43961729E+00 -9.04064083E-01";

    double expected_values[] = {
            1.38394898E-17,  1.38394898E-17, 2.26016022E-01 ,
            2.09095447E-16, 1.43961729E+00, -9.04064083E-01,
            -4.96113637E-16, -1.43961729E+00, -9.04064083E-01
    };

    rewind(stream);
    fputs(sec, stream);
    rewind(stream);

    double* values;
    size_t sz;
    stdl_lexer* lx = NULL;

    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &values));

    TEST_ASSERT_EQUAL_INT(9, sz);
    for (int j = 0; j < 9; ++j) {
        TEST_ASSERT_EQUAL_DOUBLE(expected_values[j], values[j]);
    }

    free(values);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}

void test_parser_incorrect_vector_numbers_ko() {
    char* incorrect[] = {
            "2", // just a number
            "2\n", // no value
            "2\n1.", // missing one value
            "2\n1.  x", // value is not integer
            "2\n 1. 1.e", // incorrect number
            "2\n 1. 1.x", // should end with \n
    };

    double* values;
    size_t sz;
    stdl_lexer* lx = NULL;

    for(int i = 0; i < 6; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
        STDL_NOK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &values));
        ASSERT_STDL_OK(stdl_lexer_delete(lx));
    }

}


void test_parser_vector_string_ok() {
    char* sec =
            "6\n"
            "#p wB97XD/6-311+G(d) Opt=tight freq scrf=(SMD,solvent=aceton\n"
            "itrile)     ";

    char* expected = "#p wB97XD/6-311+G(d) Opt=tight freq scrf=(SMD,solvent=acetonitrile)     ";

    rewind(stream);
    fputs(sec, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));

    char* actual;
    size_t sz;
    ASSERT_STDL_OK(stdl_fchk_parser_get_vector_string(lx, &sz, &actual));

    TEST_ASSERT_EQUAL_INT(6, sz);
    TEST_ASSERT_EQUAL_STRING(expected, actual);
    TEST_ASSERT_EQUAL_INT(6 * 12, strlen(actual));

    free(actual);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}

void test_read_fchk_skip_all_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    char* name = NULL;
    char type;
    int is_scalar = -1, nsections = 0;

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    while (lx->current_tk_type != STDL_TK_EOF) {
        ASSERT_STDL_OK(stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar));
        free(name);

        ASSERT_STDL_OK(stdl_fchk_parser_skip_section(lx, type, is_scalar));
        nsections++;
    }

    TEST_ASSERT_EQUAL_INT(86, nsections);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_read_fchk_read_all_section_ok() {

    char* fchk_path = "../tests/test_files/water_631g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    char* name = NULL;
    char type;
    int is_scalar = -1, nsections = 0;
    long an_integer;
    long* vector_integers = NULL;
    double a_number;
    double* vector_numbers = NULL;
    char* a_string;
    size_t sz;

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    while (lx->current_tk_type != STDL_TK_EOF) {
        ASSERT_STDL_OK(stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar));
        free(name);

        if(is_scalar) {
            switch (type) {
                case 'I':
                    ASSERT_STDL_OK(stdl_fchk_parser_get_scalar_integer(lx, &an_integer));
                    break;
                case 'R':
                    ASSERT_STDL_OK(stdl_fchk_parser_get_scalar_number(lx, &a_number));
                    break;
                default:
                    TEST_FAIL_MESSAGE("Unknown type for scalar");
                    break;
            }
        } else {
            switch (type) {
                case 'I':
                    ASSERT_STDL_OK(stdl_fchk_parser_get_vector_integers(lx, &sz, &vector_integers));
                    free(vector_integers);
                    break;
                case 'R':
                    ASSERT_STDL_OK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &vector_numbers));
                    free(vector_numbers);
                    break;
                case 'C':
                    ASSERT_STDL_OK(stdl_fchk_parser_get_vector_string(lx, &sz, &a_string));
                    free(a_string);
                    break;
                default:
                    TEST_FAIL_MESSAGE("Unknown type for vector");
                    break;
            }
        }

        nsections++;
    }

    TEST_ASSERT_EQUAL_INT(86, nsections);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}


void test_extract_wavefunction_and_basis_ok() {

    char* fchk_path = "../tests/test_files/water_sto3g.fchk";

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer *lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction *wf = NULL;
    stdl_basis *bs = NULL;
    ASSERT_STDL_OK(stdl_fchk_parser_extract(lx, &wf, &bs));
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    // check basis and geometry
    TEST_ASSERT_EQUAL_INT(3, wf->natm);
    TEST_ASSERT_EQUAL_INT(7, wf->nao); // O[1s,2s,2px,2py,2pz] + H[1s] + H[1s]
    TEST_ASSERT_EQUAL_INT(7, wf->nmo);
    TEST_ASSERT_EQUAL_INT(5, wf->nocc);

    ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_extract_wf_too_short_ko() {
    char* vector = "whatever\nwhatever\nAtomic numbers                             I   N=           3";

    rewind(stream);
    fputs(vector, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));

    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction *wf = NULL;
    stdl_basis *bs = NULL;
    STDL_NOK(stdl_fchk_parser_extract(lx, &wf, &bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}
