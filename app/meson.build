exe_args = []

# dependencies
app_dep = [stdlite_lib_dep]

argtable_proj = subproject('argtable')
argtable_dep = argtable_proj.get_variable('argtable_dep')
app_dep += argtable_dep

tomlc99_proj = subproject('tomlc99')
tomlc99_dep = tomlc99_proj.get_variable('tomlc99_dep')
app_dep += tomlc99_dep

# sources
srcs = files('main.c')

# executable
stdlite_exe = executable(
    meson.project_name() + '_run',
    sources: srcs,
    dependencies: app_dep,
    install: true,
)

# tests
test('help', stdlite_exe, args: ['-h'], suite: 'app')

test('noargs', stdlite_exe, suite: 'app', should_fail: true)
test('unknown arg', stdlite_exe, args: ['-k'], suite: 'app', should_fail: true)
test('test input', stdlite_exe, args: ['../tests/test_files/test.toml'], suite: 'app')