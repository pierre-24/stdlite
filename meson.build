project('stdlite', 'c',
  version : '0.1',
  default_options : ['warning_level=3'],
)

cc = meson.get_compiler('c')

# configure dependencies
project_dep = []

build_args = [
  '-DPROJECT_NAME=' + meson.project_name(),
  '-DPROJECT_VERSION=' + meson.project_version(),
]

libcint_dep = cc.find_library('libcint', required: false)
if not libcint_dep.found()
  libcint_proj = subproject('libcint', default_options: [])
  libcint_dep = libcint_proj.get_variable('libcint_dep')
endif
project_dep += libcint_dep

la_backend = get_option('la_backend')
if la_backend == 'openblas' or la_backend == 'default'
  project_dep += cc.find_library('openblas', required : true)
  project_dep += cc.find_library('lapack', required : true)
endif

# fetch sources and includes
srcs = []
subdir('src')
includes = include_directories('include')

# library
stdlite_lib = library(
  meson.project_name(),
  srcs,
  install : true,
  c_args : build_args,
  include_directories : includes,
  dependencies: project_dep,
)

stdlite_lib_dep = declare_dependency(
  link_with: stdlite_lib,
  include_directories: includes,
  dependencies: project_dep,
)

# app
subdir('app')

# tests
subdir('tests')
