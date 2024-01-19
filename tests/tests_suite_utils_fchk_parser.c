#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>

#include "tests_suite.h"

FILE* stream;

void setUp(void) {
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
    stdl_lexer* lx;

    rewind(stream);
    fputs(scalar, stream);
    rewind(stream);

    lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

    TEST_ASSERT_EQUAL_STRING("Multiplicity", title);
    TEST_ASSERT_EQUAL_CHAR('I', type);
    TEST_ASSERT_EQUAL_INT(1, is_scalar);

    TEST_ASSERT_EQUAL_INT(STDL_TK_DIGIT, lx->current_tk_type);

    free(title);

    STDL_OK(stdl_lexer_delete(lx));

    // vector
    char* vector = "Atomic numbers                             I   N=           3";

    rewind(stream);
    fputs(vector, stream);
    rewind(stream);

    lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

    TEST_ASSERT_EQUAL_STRING("Atomic numbers", title);
    TEST_ASSERT_EQUAL_CHAR('I', type);
    TEST_ASSERT_EQUAL_INT(0, is_scalar);

    TEST_ASSERT_EQUAL_INT(STDL_TK_DIGIT, lx->current_tk_type);

    free(title);

    STDL_OK(stdl_lexer_delete(lx));
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
    stdl_lexer* lx;

    for(int i=0; i < 7; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);

        STDL_NOK(stdl_fchk_parser_get_section_info(lx, &title, &type, &is_scalar));

        free(title);

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_parser_vector_ints_ok() {
    char* fmt =
            "%d\n"
            "           0           1           2           3           4           5\n"
            "           6           7";
    int to_test[] = {6, 8};

    stdl_lexer* lx;

    for(int i = 0; i < 2; i++) {
        int mx = to_test[i];
        rewind(stream);
        fprintf(stream, fmt, mx);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);

        long* values;
        size_t sz;
        STDL_OK(stdl_fchk_parser_get_vector_ints(lx, &sz, &values));

        TEST_ASSERT_EQUAL_INT(mx, sz);
        for (long j = 0; j < mx; ++j) {
            TEST_ASSERT_EQUAL_INT(j, values[j]);
        }

        free(values);

        STDL_OK(stdl_lexer_delete(lx));
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

    stdl_lexer* lx;
    long* values;
    size_t sz;

    for(int i = 0; i < 5; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);

        STDL_NOK(stdl_fchk_parser_get_vector_ints(lx, &sz, &values));

        STDL_OK(stdl_lexer_delete(lx));
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

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    double* values;
    size_t sz;
    STDL_OK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &values));

    TEST_ASSERT_EQUAL_INT(9, sz);
    for (int j = 0; j < 9; ++j) {
        TEST_ASSERT_EQUAL_DOUBLE(expected_values[j], values[j]);
    }

    free(values);

    STDL_OK(stdl_lexer_delete(lx));
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

    stdl_lexer* lx;
    double* values;
    size_t sz;

    for(int i = 0; i < 6; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);

        STDL_NOK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &values));

        STDL_OK(stdl_lexer_delete(lx));
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

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    char* actual;
    size_t sz;
    STDL_OK(stdl_fchk_parser_get_vector_string(lx, &sz, &actual));

    TEST_ASSERT_EQUAL_INT(6, sz);
    TEST_ASSERT_EQUAL_STRING(expected, actual);
    TEST_ASSERT_EQUAL_INT(6 * 12, strlen(actual));

    free(actual);

    STDL_OK(stdl_lexer_delete(lx));
}

void test_read_fchk_skip_all_ok() {
    char cwd[512], fchk_path[1024];
    TEST_ASSERT_NOT_NULL(getcwd(cwd, 512));

    sprintf(fchk_path, "%s/../../tests/test_files/water.fchk", cwd);

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    char* name = NULL;
    char type;
    int is_scalar = -1, nsections = 0;

    stdl_lexer* lx = stdl_lexer_new(f);
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_begin(lx));

    while (lx->current_tk_type != STDL_TK_EOF) {
        STDL_OK(stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar));
        free(name);

        STDL_OK(stdl_fchk_parser_skip_section(lx, type, is_scalar));
        nsections++;
    }

    TEST_ASSERT_EQUAL_INT(86, nsections);

    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}

void test_read_fchk_read_all_section_ok() {
    char cwd[512], fchk_path[1024];
    TEST_ASSERT_NOT_NULL(getcwd(cwd, 512));

    sprintf(fchk_path, "%s/../../tests/test_files/water.fchk", cwd);

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

    stdl_lexer* lx = stdl_lexer_new(f);
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_begin(lx));

    while (lx->current_tk_type != STDL_TK_EOF) {
        STDL_OK(stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar));
        free(name);

        if(is_scalar) {
            switch (type) {
                case 'I':
                    STDL_OK(stdl_fchk_parser_get_scalar_int(lx, &an_integer));
                    break;
                case 'R':
                    STDL_OK(stdl_fchk_parser_get_scalar_number(lx, &a_number));
                    break;
                default:
                    TEST_FAIL_MESSAGE("Unknown type for scalar");
                    break;
            }
        } else {
            switch (type) {
                case 'I':
                    STDL_OK(stdl_fchk_parser_get_vector_ints(lx, &sz, &vector_integers));
                    free(vector_integers);
                    break;
                case 'R':
                    STDL_OK(stdl_fchk_parser_get_vector_numbers(lx, &sz, &vector_numbers));
                    free(vector_numbers);
                    break;
                case 'C':
                    STDL_OK(stdl_fchk_parser_get_vector_string(lx, &sz, &a_string));
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

    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
