#include <stdlib.h>

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

void test_parser_vector_ints() {
    char* fmt =
            "%d\n"
            "           0           1           2           3           4           5\n"
            "           6           7";
    int to_test[] = {1, 2, 5, 7};

    stdl_lexer* lx;

    for(int i = 0; i < 4; i++) {
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

void test_parser_vector_numbers() {
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
