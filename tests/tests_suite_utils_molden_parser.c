#include "tests_suite.h"

#include <stdlite/utils/matrix.h>
#include <stdlite/utils/molden_parser.h>

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}

void test_read_title_ok() {
    stdl_lexer* lx = NULL;

    char* title = "test", *rtitle = NULL;
    char buff[512];
    sprintf(buff, "[%s]", title);

    rewind(stream);
    fputs(buff, stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));

    ASSERT_STDL_OK(stdl_molden_parser_read_section_title(lx, &rtitle));

    TEST_ASSERT_EQUAL_STRING(title, rtitle);

    free(rtitle);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}

void test_parse_incorrect_title_ko() {
    char* incorrect[] = {
            "[x", // too short
            "[x\ny]", // contains \n
            " [x", // not the first character of the line
    };

    char* title = NULL;
    stdl_lexer* lx = NULL;

    for(int i=0; i < 3; i++) {
        rewind(stream);
        fputs(incorrect[i], stream);
        rewind(stream);

        ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));
        STDL_NOK(stdl_molden_parser_read_section_title(lx, &title));

        free(title);

        ASSERT_STDL_OK(stdl_lexer_delete(lx));
    }
}


void test_skip_section_ok() {
    stdl_lexer* lx = NULL;

    char* title = "test", *rtitle = NULL;
    char buff[512];
    sprintf(buff, "[whatever]\nwhatever\n[%s]", title);

    rewind(stream);
    fputs(buff, stream);
    rewind(stream);

    ASSERT_STDL_OK(stdl_lexer_new(stream, &lx));

    ASSERT_STDL_OK(stdl_molden_parser_read_section_title(lx, &rtitle));
    TEST_ASSERT_EQUAL_STRING("whatever", rtitle);
    free(rtitle);

    ASSERT_STDL_OK(stdl_molden_parser_skip_section(lx));

    ASSERT_STDL_OK(stdl_molden_parser_read_section_title(lx, &rtitle));
    TEST_ASSERT_EQUAL_STRING(title, rtitle);
    free(rtitle);

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
}

void test_read_file_ok() {
    char* molden_path = "../tests/test_files/water_631g.molden";

    FILE* f = fopen(molden_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));

    stdl_wavefunction* wf = NULL;
    stdl_basis* bs = NULL;

    ASSERT_STDL_OK(stdl_molden_parser_extract(lx, &wf, &bs));

    // ASSERT_STDL_OK(stdl_wavefunction_delete(wf));
    ASSERT_STDL_OK(stdl_basis_delete(bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));
    fclose(f);
}
