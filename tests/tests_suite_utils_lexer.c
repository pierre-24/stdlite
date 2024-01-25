#include <string.h>

#include <stdlite/utils/lexer.h>

#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}

void test_lexer_ok() {
    char* str = "ab[-9]=3!";
    int l = strlen(str);
    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, stream));

    stdl_token_type tab[] = {
        STDL_TK_ALPHA, STDL_TK_ALPHA, STDL_TK_CHAR, STDL_TK_DASH, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EQ, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EOF};

    for(int i=0; i <= l; i++) {
        TEST_ASSERT_EQUAL_INT(tab[i], lx->current_tk_type);
        stdl_lexer_advance(lx, 1);
    }

    // if one continue to advance, it results in an error
    STDL_NOK(stdl_lexer_advance(lx, 1));
    TEST_ASSERT_EQUAL_INT(STDL_TK_EOF, lx->current_tk_type);

    STDL_OK(stdl_lexer_delete(lx));
}

void test_lexer_line_ok() {
    char* str = "a\nb1\ncde\nf2\0";
    int l = strlen(str);

    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(&lx, stream));

    int line = 1, pos_in_line = 0;

    for(int i=0; i < l; i++) {
        TEST_ASSERT_EQUAL_INT(line, lx->current_line);
        TEST_ASSERT_EQUAL_INT(pos_in_line, lx->current_pos_in_line);

        stdl_lexer_advance(lx, 1);
        if(lx->current_tk_value == '\n') {
            line += 1;
            pos_in_line = 0;
        } else {
            pos_in_line += 1;
        }
    }

    STDL_OK(stdl_lexer_delete(lx));
}