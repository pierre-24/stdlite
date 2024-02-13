#include <stdlite/utils/fchk_parser.h>
#include <string.h>

#include "tests_suite.h"


void read_fchk(char* fchk_path, stdl_wavefunction** wf, stdl_basis** bs) {
    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    ASSERT_STDL_OK(stdl_lexer_new(f, &lx));
    ASSERT_STDL_OK(stdl_fchk_parser_skip_intro(lx));

    ASSERT_STDL_OK(stdl_fchk_parser_extract(lx, wf, bs));

    ASSERT_STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
