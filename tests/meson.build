unity_proj = subproject('Unity')
unity_dep = unity_proj.get_variable('unity_dep')

test('a_test',
     executable(
        'a_test',
        files('test_dummy.c'),
        dependencies: [stdlite_lib_dep, unity_dep],
        install: false
     )
)