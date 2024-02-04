#include <string.h>

#include <stdlite/utils/lexer.h>
#include <ctype.h>

#include "tests_suite.h"

FILE* stream;

void setUp(void) {
    stdl_set_debug_level(-1);

    stream = tmpfile();
}

void tearDown(void) {
    fclose(stream);
}

void test_lexer_ok() {
    char* str = "ab[-9]=3!";
    size_t l = strlen(str);
    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(stream, &lx));

    stdl_token_type tab[] = {
        STDL_TK_ALPHA, STDL_TK_ALPHA, STDL_TK_CHAR, STDL_TK_DASH, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EQ, STDL_TK_DIGIT, STDL_TK_CHAR, STDL_TK_EOF};

    for(size_t i=0; i <= l; i++) {
        TEST_ASSERT_EQUAL_INT(tab[i], lx->current_tk_type);
        if(lx->current_tk_type != STDL_TK_EOF)
            STDL_OK(stdl_lexer_advance(lx, 1));
    }

    // if one continue to advance, it results in an error!
    STDL_NOK(stdl_lexer_advance(lx, 1));
    TEST_ASSERT_EQUAL_INT(STDL_TK_EOF, lx->current_tk_type);

    STDL_OK(stdl_lexer_delete(lx));
}

void test_lexer_line_ok() {
    char* str = "a\nb1\ncde\nf2\0";
    size_t l = strlen(str);

    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(stream, &lx));

    int line = 1, pos_in_line = 0;

    for(size_t i=0; i < l; i++) {
        TEST_ASSERT_EQUAL_INT(line, lx->current_line);
        TEST_ASSERT_EQUAL_INT(pos_in_line, lx->current_pos_in_line);

        STDL_OK(stdl_lexer_advance(lx, 1));
        if(lx->current_tk_value == '\n') {
            line += 1;
            pos_in_line = 0;
        } else {
            pos_in_line += 1;
        }
    }

    STDL_OK(stdl_lexer_delete(lx));
}

void test_lexer_eat_ok() {
    char* str = "a\n";

    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(stream, &lx));

    STDL_OK(stdl_lexer_eat(lx, STDL_TK_ALPHA)); // works on `a`, advance
    STDL_NOK(stdl_lexer_eat(lx, STDL_TK_ALPHA)); // does not work on `\n`

    STDL_OK(stdl_lexer_delete(lx));
}


void test_lexer_skip_ok() {
    char* str = "abc1\n";

    fputs(str, stream);
    rewind(stream);

    stdl_lexer* lx = NULL;
    STDL_OK(stdl_lexer_new(stream, &lx));

    STDL_OK(stdl_lexer_skip(lx, isalpha)); // advance up to `1`

    TEST_ASSERT_EQUAL_INT(3, lx->pos_in_stream);
    TEST_ASSERT_EQUAL_INT(3, lx->current_pos_in_line);
    TEST_ASSERT_EQUAL_CHAR(lx->current_tk_value, '1');

    STDL_OK(stdl_lexer_delete(lx));
}

