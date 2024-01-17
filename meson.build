project('stdlite', 'c',
  version : '0.1',
  default_options : ['warning_level=3'],
)

# configure dependencies
project_deps = []

# fetch sources and includes
srcs = []
subdir('src')
includes = include_directories('include')

# library
build_args = [
  '-DPROJECT_NAME=' + meson.project_name(),
  '-DPROJECT_VERSION=' + meson.project_version(),
]

stdlite_lib = library(
  meson.project_name(),
  srcs,
  install : true,
  c_args : build_args,
  include_directories : includes,
)

stdlite_lib_deps = declare_dependency(
  link_with: stdlite_lib,
  include_directories: includes,
  dependencies: project_deps,
)

# app
subdir('app')
