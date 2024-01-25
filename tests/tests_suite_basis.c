#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>

#include "tests_suite.h"

void test_basis_functions_ovlp_ok() {
    char cwd[512], fchk_path[1024];
    TEST_ASSERT_NOT_NULL(getcwd(cwd, 512));

    sprintf(fchk_path, "%s/../tests/test_files/water_631g.fchk", cwd);

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, f));
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = NULL;
    stdl_basis* bs = NULL;
    STDL_OK(stdl_fchk_parser_wavefunction_new(&wf, &bs, lx));

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
