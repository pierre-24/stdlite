project('stdlite', 'c',
  version : '0.2.1',
  default_options : ['warning_level=3'],
)

cc = meson.get_compiler('c')

# configure dependencies
project_dep = []

build_args = [
  '-DPROJECT_NAME="lib' + meson.project_name() + '"',
  '-DPROJECT_VERSION="' + meson.project_version() + '"',
]

project_info = ''
date_program = find_program('date', required: false)
if date_program.found()
  build_args += '-DPROJECT_BUILD_DATE="' + run_command(date_program).stdout().replace('\n', '') + '"'
endif

git_program = find_program('git', required: false)
if git_program.found()
  build_args += '-DPROJECT_BUILD_COMMIT="' + run_command(git_program, 'rev-parse', 'HEAD').stdout().substring(0, 7) + '"'
endif


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
    dep_mkl = dependency('mkl-dynamic-ilp64-iomp', required: false)
    if dep_mkl.found()
      project_dep += dep_mkl
    else
      project_dep += cc.find_library('mkl_core')
      project_dep += cc.find_library('mkl_intel_ilp64')
      project_dep += cc.find_library('mkl_intel_thread')
      project_dep += cc.find_library('iomp5')
    endif
  else
    dep_mkl = dependency('mkl-dynamic-ilp64-seq', required: false)
    if dep_mkl.found()
      project_dep += dep_mkl
    else
      project_dep += cc.find_library('mkl_core')
      project_dep += cc.find_library('mkl_intel_ilp64')
      project_dep += cc.find_library('mkl_sequential')
    endif
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
if get_option('tests')
  subdir('tests')
endif