#include <ctype.h>
#include <stdlib.h>

#include <stdlite/utils/base_parser.h>

#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}

void wipe_and_write_integer(FILE* f, const long x) {
    rewind(f);
    fprintf(f, "%ld", x);
    fputc(' ', f);
    rewind(f);
}

void test_read_integer_ok() {
    long test_integers[] = {1, -5, 17, 64555};
    long expected, actual;
    stdl_lexer* lx;

    for(int i=0; i < 4; i++) {
        expected = test_integers[i];
        wipe_and_write_integer(stream, expected);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_OK(stdl_parser_get_integer(lx, &actual));

        TEST_ASSERT_EQUAL_INT(expected, actual);

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_read_non_integer_ko() {
    char* wrong[] = {
        ".",
        "x",
        "x1",
    };

    stdl_lexer* lx;
    long actual;

    for (int i=0; i < 3; i++) {
        rewind(stream);
        fputs(wrong[i], stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_NOK(stdl_parser_get_integer(lx, &actual)); // that does not works ;)

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_read_integers_ok() {
    long test_integers[] = {1, -5, 17};
    rewind(stream);
    fprintf(stream , "%ld %ld %ld ", test_integers[0], test_integers[1], test_integers[2]);
    rewind(stream);

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    long expected, actual;

    for(int i=0; i < 3; i++) {
        expected = test_integers[i];

        // read number
        STDL_OK(stdl_parser_get_integer(lx, &actual));
        TEST_ASSERT_EQUAL_INT(expected, actual);

        // eat up next space
        STDL_OK(stdl_lexer_eat(lx, STDL_TK_WHITESPACE));
    }

    STDL_OK(stdl_lexer_delete(lx));
}

void wipe_and_write_real(FILE* f, char* fmt, const double x) {
    rewind(f);
    fprintf(f, fmt, x);
    fputc(' ', f);
    rewind(f);
}

void test_read_real_ok() {
    double test_reals[] = {1.25, -5.3, 17.1, 6.5e3};
    double expected, actual;
    stdl_lexer* lx;

    for(int i=0; i < 4; i++) {
        expected = test_reals[i];

        // pure float
        wipe_and_write_real(stream, "%f", expected);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_OK(stdl_parser_get_real(lx, &actual));

        TEST_ASSERT_EQUAL_DOUBLE(expected, actual);

        STDL_OK(stdl_lexer_delete(lx));

        // scientific notation
        wipe_and_write_real(stream, "%e", expected);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_OK(stdl_parser_get_real(lx, &actual));

        TEST_ASSERT_EQUAL_DOUBLE(expected, actual);

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_read_real_many_notations_ok() {
    double expected = 1.25, actual;
    stdl_lexer* lx;

    char* nt[] = { // many forms of 1.25
        "1.25",
        "+1.25",
        "0.125e1",
        ".125e1",
        ".125e+1",
        "12.5e-1",
        // also, stops before it gets weird:
        "1.25x",
        "1.25.",
        ".125e1."
    };

    for(int i=0; i < 9; i++) {
        rewind(stream);
        fputs(nt[i], stream);
        fputc(' ', stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_OK(stdl_parser_get_real(lx, &actual));

        TEST_ASSERT_EQUAL_DOUBLE(expected, actual);

        STDL_OK(stdl_lexer_delete(lx));
    }

    expected = -0.6;
    char* nt2[] = { // many forms of -0.6
        "-0.6",
        "-.6",
        "-.06e1",
        "-.06e+1",
        "-6e-1",
        "-6.e-1",
    };

    for(int i=0; i < 6; i++) {
        rewind(stream);
        fputs(nt2[i], stream);
        fputc(' ', stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_OK(stdl_parser_get_real(lx, &actual));

        TEST_ASSERT_EQUAL_DOUBLE(expected, actual);

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_read_non_real_ko() {
    char* wrong[] = {
        ".",
        "x",
        ".x",
        "1e",
        "1e."
    };

    stdl_lexer* lx;
    double actual;

    for (int i=0; i < 5; i++) {
        rewind(stream);
        fputs(wrong[i], stream);
        rewind(stream);

        lx = stdl_lexer_new(stream);
        TEST_ASSERT_NOT_NULL(lx);
        STDL_NOK(stdl_parser_get_real(lx, &actual)); // that does not works!

        STDL_OK(stdl_lexer_delete(lx));
    }
}

void test_read_litteral_ok() {
    stdl_lexer* lx;
    char* expected = "keyword";

    rewind(stream);
    fputs(expected, stream);
    rewind(stream);

    lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    char* actual;
    STDL_OK(stdl_parser_get_literal(lx, isalpha, &actual));
    TEST_ASSERT_EQUAL_STRING(expected, actual);

    free(actual);

    STDL_OK(stdl_lexer_delete(lx));
}

void test_read_litterals_ok() {
    char* test_litterals[] = {"this", "is", "nice"};
    rewind(stream);
    fprintf(stream , "%s %s %s ", test_litterals[0], test_litterals[1], test_litterals[2]);
    rewind(stream);

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    char* expected;
    char* actual = NULL;

    for(int i=0; i < 3; i++) {
        expected = test_litterals[i];

        // read
        STDL_OK(stdl_parser_get_literal(lx, isalpha, &actual));
        TEST_ASSERT_EQUAL_STRING(expected, actual);

        free(actual);

        // eat up next space
        STDL_OK(stdl_lexer_eat(lx, STDL_TK_WHITESPACE));
    }

    STDL_OK(stdl_lexer_delete(lx));
}
