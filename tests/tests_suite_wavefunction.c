#include <string.h>
#include <unistd.h>

#include <stdlite/utils/fchk_parser.h>
#include <stdlite/matrix.h>

#include "tests_suite.h"

void test_extract_wavefunction_from_fchk_ok() {
    char cwd[512], fchk_path[1024];
    TEST_ASSERT_NOT_NULL(getcwd(cwd, 512));

    sprintf(fchk_path, "%s/../tests/test_files/water_sto3g.fchk", cwd);

    FILE* f = fopen(fchk_path, "r");
    TEST_ASSERT_NOT_NULL(f);

    stdl_lexer* lx = stdl_lexer_new(f);
    TEST_ASSERT_NOT_NULL(lx);

    STDL_OK(stdl_fchk_parser_skip_intro(lx));

    stdl_wavefunction * wf = stdl_fchk_parser_wavefunction_new(lx);
    TEST_ASSERT_NOT_NULL(wf);

    TEST_ASSERT_EQUAL_INT(3, wf->natm);
    TEST_ASSERT_EQUAL_INT(7, wf->nao); // O[1s,2s,2px,2py,2pz] + H[1s] + H[1s]
    TEST_ASSERT_EQUAL_INT(7, wf->nmo);
    TEST_ASSERT_EQUAL_INT(10, wf->nelec);

    stdl_matrix_numbers_print(wf->nmo, wf->nao, wf->C, 0);

    STDL_OK(stdl_wavefunction_delete(wf));
    STDL_OK(stdl_lexer_delete(lx));

    fclose(f);
}
