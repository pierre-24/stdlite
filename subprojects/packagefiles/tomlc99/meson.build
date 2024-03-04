# this is a minimal project file for tomlc99 (https://github.com/cktan/tomlc99).

project(
    'tomlc99', 'c',
    version : '0.1', # ?!?
)

cc = meson.get_compiler('c')

# sources & includes
srcs = ['toml.c']
include_dir = ['.']

# make a library
tomlc99_lib = library(
    meson.project_name(),
    srcs,
    include_directories: include_dir,
    install: not meson.is_subproject(),
)

tomlc99_dep = declare_dependency(
    link_with: tomlc99_lib,
    include_directories: include_dir
)