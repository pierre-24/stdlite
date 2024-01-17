//
// Created by pierre on 1/17/24.
//

#include <unity.h>
#include "dummy.h"


void setUp(void) {
}

void tearDown(void) {
}

void test_x() {
    TEST_ASSERT_EQUAL_INT(8, test_func(4));
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_x);
    return UNITY_END();
}

