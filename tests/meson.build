unity_proj = subproject('Unity')
unity_dep = unity_proj.get_variable('unity_dep')

tests_to_run = [
    'utils_lexer'
]

foreach t: tests_to_run
    test('tests_suite_' + t,
         executable(
             'tests_suite_' +  t,
             files('tests_suite_' +  t + '.c'),
             dependencies: [stdlite_lib_dep, unity_dep],
             install: false
         )
    )
endforeach

