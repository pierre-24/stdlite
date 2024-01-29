# this is a minimal project file for argtable3 (https://github.com/argtable/argtable3/).
# It uses the amalgamed files.

project(
    'argtable', 'c',
    version : '3.2.2',
)

cc = meson.get_compiler('c')

# sources & includes & dependencies
srcs = ['dist/argtable3.c']

project_dep = [
    cc.find_library('m', required : true)
]

includes = include_directories('dist')

# make a library
argtable_lib = library(
    meson.project_name(),
    srcs,
    include_directories: includes,
    dependencies: project_dep,
    install: not meson.is_subproject(),
)

argtable_dep = declare_dependency(
    link_with: argtable_lib,
    include_directories: includes,
    dependencies: project_dep,
)