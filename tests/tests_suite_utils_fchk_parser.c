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

    free(title);

    STDL_OK(stdl_lexer_delete(lx));
}
