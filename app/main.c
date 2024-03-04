#include <stdio.h>
#include <stdlib.h>

#include <argtable3.h>

#include "app.h"

/**
 * Parse the arguments of the program.
 *
 * @param argc number of arguments
 * @param argv value of arguments
 * @return the exit code, so `EXIT_SUCCESS` if everything went well.
 */
int parse_arguments(int argc, char* argv[]) {
    char* self = argv[0];
    struct arg_lit* arg_help;
    struct arg_file* arg_input;
    struct arg_end* arg_end_;

    int exit_code = EXIT_SUCCESS;

    void* args_table[] = {
            arg_help = arg_litn("h", "help", 0, 1, "display this help and exit"),
            arg_input = arg_file1(NULL, NULL, "<input>", "input file (TOML format)"),
            arg_end_ = arg_end(20)
    };

    int nerrors = arg_parse(argc, argv, args_table);

    /* `--help` takes precedence over the rest */
    if (arg_help->count > 0) {
        printf("Usage: %s", self);
        arg_print_syntax(stdout, args_table, "\n");
        printf(STDL_APP_ARG_DOC);
        arg_print_glossary(stdout, args_table, "  %-25s %s\n");

        exit_code = EXIT_SUCCESS;
    }  else if (nerrors > 0) {
        arg_print_errors(stderr, arg_end_, self);
        exit_code = EXIT_FAILURE;
    }

    arg_freetable(args_table, sizeof(args_table) / sizeof(args_table[0]));

    return exit_code;
}

int main(int argc, char* argv[]) {

    int exit_code = parse_arguments(argc, argv);
    if(exit_code != EXIT_SUCCESS)
        return exit_code;

    return EXIT_SUCCESS;
}
