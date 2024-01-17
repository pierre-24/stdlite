exe_args = []

# executable
stdlite_exe = executable(
    meson.project_name() + '_run',
    sources: files(
        'main.c'
    ),
    dependencies: stdlite_lib_deps,
    c_args: exe_args,
    link_language: 'c',
    install: true,
)

# tests
test('noargs', stdlite_exe, suite: 'app')