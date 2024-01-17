#ifndef TESTS_SUITE_H
#define TESTS_SUITE_H

#include <unity.h>

#include <stdlite/errors.h>

#define _OK(v) TEST_ASSERT_EQUAL_INT(v, STDL_ERR_OK)
#define _NOK(v) TEST_ASSERT_NOT_EQUAL_INT(v, STDL_ERR_OK)


#endif //TESTS_SUITE_H
