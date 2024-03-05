exe_args = []

# build args
app_build_args = [
    '-DAPP_NAME="' + meson.project_name() + '_run"',
    '-DAPP_VERSION="' + meson.project_version() + '"',
]

# dependencies
app_dep = [stdlite_lib_dep]

argtable_dep = cc.find_library('argtable3', required: false)
if not argtable_dep.found()
    argtable_proj = subproject('argtable')
    argtable_dep = argtable_proj.get_variable('argtable_dep')
endif
app_dep += argtable_dep

tomlc99_dep = cc.find_library('toml', required: false)
if not tomlc99_dep.found()
    tomlc99_proj = subproject('tomlc99')
    tomlc99_dep = tomlc99_proj.get_variable('tomlc99_dep')
endif
app_dep += tomlc99_dep

# sources
srcs = files(
    'user_input.c'
)

includes = include_directories('.')

# library
stdlite_app_lib = static_library(
    meson.project_name() + '_app',
    srcs,
    include_directories : includes,
    dependencies: app_dep,
)

stdlite_app_lib_dep = declare_dependency(
    link_with: stdlite_app_lib,
    include_directories: includes,
    dependencies: app_dep,
)

# executable
stdlite_exe = executable(
    meson.project_name() + '_run',
    sources: ['main.c'],
    dependencies: app_dep,
    link_with: stdlite_app_lib,
    install: true,
    c_args : app_build_args,
)

subdir('tests')