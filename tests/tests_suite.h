#ifndef TESTS_SUITE_H
#define TESTS_SUITE_H

#define UNITY_INCLUDE_DOUBLE
#include <unity.h>

#include <stdlite/logging.h>
#include <stdlite/wavefunction.h>
#include <stdlite/basis.h>

#define ASSERT_STDL_OK(v) TEST_ASSERT_EQUAL_INT(STDL_ERR_OK, v)
#define STDL_NOK(v) TEST_ASSERT_NOT_EQUAL_INT(STDL_ERR_OK, v)

void read_fchk(char* fchk_path, stdl_wavefunction** wf, stdl_basis** bs);


#endif //TESTS_SUITE_H
