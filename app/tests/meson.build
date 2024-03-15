# -- code tests
# dependencies
tests_dep = [stdlite_app_lib_dep]

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
    'user_input_handler',
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
             c_args : app_build_args,
         ),
        suite: 'app_code',
    )
endforeach

# -- app test
test('help', stdlite_exe, args: ['-h'], suite: 'app')
test('valid input (TOML)', stdlite_exe, args: ['../app/tests/test_files/test.toml'], suite: 'app')
test('valid input (context)', stdlite_exe, args: ['--ctx_source=../app/tests/test_files/context_water_sto3g.h5', '--ctx_source_type=STDL_CTX'], suite: 'app')

test('noargs', stdlite_exe, suite: 'app', should_fail: true)
test('unknown arg', stdlite_exe, args: ['-k'], suite: 'app', should_fail: true)
test('invalid TOML', stdlite_exe, args: ['../app/tests/test_files/wrong.toml'], suite: 'app', should_fail: true)
test('invalid TOML 2', stdlite_exe, args: ['../app/tests/test_files/wrong2.toml'], suite: 'app', should_fail: true)