#ifndef STDLITE_TESTS_SUITE_H
#define STDLITE_TESTS_SUITE_H

#include <unity.h>
#include <stdlite/logging.h>

#define ASSERT_STDL_OK(v) TEST_ASSERT_EQUAL_INT(STDL_ERR_OK, v)

#endif //STDLITE_TESTS_SUITE_H
