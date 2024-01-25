#ifndef TESTS_SUITE_H
#define TESTS_SUITE_H

#define UNITY_INCLUDE_DOUBLE
#include <unity.h>

#include <stdlite/errors.h>

#define STDL_OK(v) TEST_ASSERT_EQUAL_INT(STDL_ERR_OK, v)
#define STDL_NOK(v) TEST_ASSERT_NOT_EQUAL_INT(STDL_ERR_OK, v)


#endif //TESTS_SUITE_H
