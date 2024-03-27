project('stdlite', 'c',
  version : '0.2.0',
  default_options : ['warning_level=3'],
)

cc = meson.get_compiler('c')

# configure dependencies
project_dep = []

build_args = [
  '-DPROJECT_NAME="lib' + meson.project_name() + '"',
  '-DPROJECT_VERSION="' + meson.project_version() + '"',
]

libcint_dep = cc.find_library('libcint', required: false)
if not libcint_dep.found()
  libcint_proj = subproject('libcint', default_options: [])
  libcint_dep = libcint_proj.get_variable('libcint_dep')
endif
project_dep += libcint_dep

if get_option('openmp')
  build_args += '-DUSE_OPENMP'
  project_dep += dependency('openmp')
endif

# note: cannot use flexiblas at the moment
# see https://github.com/mpimd-csc/flexiblas/issues/2
la_backend = get_option('la_backend')
if la_backend == 'openblas' or la_backend == 'default'
  project_dep += cc.find_library('openblas', required : true)
  project_dep += cc.find_library('lapack', required : true)
  project_dep += cc.find_library('lapacke', required : true)
  build_args += '-DUSE_OPENBLAS'
elif la_backend == 'mkl'
  # found using `pkg-config --list-all | grep "mkl"`
  if get_option('openmp')
    project_dep += dependency('mkl-dynamic-ilp64-iomp', required: true)
  else
    project_dep += dependency('mkl-dynamic-ilp64-seq', required: true)
  endif
  build_args += '-DUSE_MKL'
endif

h5noserial = cc.find_library('hdf5', required: false)
if not h5noserial.found()
  project_dep += [cc.find_library('hdf5_serial', required: true), cc.find_library('hdf5_serial_hl', required: true)]
  build_args += ['-DUSE_HDF5_SERIAL']
else
  project_dep += [cc.find_library('hdf5', required: true), cc.find_library('hdf5_hl', required: true)]
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
