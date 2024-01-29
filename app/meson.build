exe_args = []

app_dep = [stdlite_lib_dep]

argtable_proj = subproject('argtable')
argtable_dep = argtable_proj.get_variable('argtable_dep')

app_dep += argtable_dep

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