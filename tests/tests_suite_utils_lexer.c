#include <string.h>

#include "tests_suite.h"

#include <stdlite/utils/lexer.h>


void setUp(void) {
}

void tearDown(void) {
}

void test_lexer_ok() {
    char* str = "ab[-9]=3!";
    int l = strlen(str);

    FILE* stream = tmpfile();
    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    stdl_token_type tab[] = {
        STDL_TK_ALPHA, STDL_TK_ALPHA, STDL_TK_CHAR, STDL_TK_DASH, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EQ, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EOS};

    for(int i=0; i <= l; i++) {
        TEST_ASSERT_EQUAL_INT(lx->current_tk_type, tab[i]);
        stdl_lexer_advance(lx, 1);
    }

    STDL_OK(stdl_lexer_delete(lx));
    fclose(stream);
}

void test_lexer_line_ok() {
    char* str = "a\nb1\ncde\nf2";
    int l = strlen(str);

    FILE* stream = tmpfile();
    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = stdl_lexer_new(stream);
    TEST_ASSERT_NOT_NULL(lx);

    int line = 1, pos_in_line = 0;

    for(int i=0; i < l; i++) {
        TEST_ASSERT_EQUAL_INT(lx->current_line, line);
        TEST_ASSERT_EQUAL_INT(lx->current_pos_in_line, pos_in_line);

        stdl_lexer_advance(lx, 1);
        if(lx->current_tk_value == '\n') {
            line += 1;
            pos_in_line = 0;
        } else {
            pos_in_line += 1;
        }
    }

    STDL_OK(stdl_lexer_delete(lx));
    fclose(stream);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_lexer_ok);
    RUN_TEST(test_lexer_line_ok);
    return UNITY_END();
}