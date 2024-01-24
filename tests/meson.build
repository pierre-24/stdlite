# dependencies
tests_dep = [stdlite_lib_dep]

unity_dep = cc.find_library('unity', required: false)
if not unity_dep.found()
    unity_proj = subproject('unity', default_options: ['support_double=true'])
    unity_dep = unity_proj.get_variable('unity_dep')
    runner_gen = unity_proj.get_variable('gen_test_runner')  # generate a runner corresponding to test suite
else
    runner_gen = generator(
        find_program('generate_test_runner.rb'),
        output: '@BASENAME@_Runner.c',
        arguments: ['@INPUT@', '@OUTPUT@']
    )
endif

tests_dep += unity_dep

# create tests
tests_to_run = [
    'utils_lexer',
    'utils_base_parser',
    'utils_fchk_parser',
]

foreach t: tests_to_run
    source_file = files('tests_suite_' +  t + '.c')

    test('tests_suite_' + t,
        executable(
            'tests_suite_' +  t,
            sources: [source_file, runner_gen.process(source_file)],
            dependencies: tests_dep,
            install: false,
            build_by_default: false,
        )
    )
endforeach

