# Translation of most of https://github.com/sunqm/libcint/blob/master/CMakeLists.txt (no test or example)
# To be used with a wrap file as a Meson subproject.

project(
    'cint', 'c',
    version : '6.1',
)

project_config = {
    'cint_VERSION': '6.1.0',
    'cint_SOVERSION': '6'
}

# configure files
configure_file(
    format: 'cmake@',
    input: 'src/config.h.in',
    output: 'config.h',
    configuration: project_config
)

configure_file(
    format: 'cmake@',
    input: 'include/cint.h.in',
    output: 'cint.h',
    configuration: project_config
)

# configure dependencies
cc = meson.get_compiler('c')

project_dep = [
    cc.find_library('m', required : true)
]

# sources & includes
srcs = files(
    'src/c2f.c', 'src/cart2sph.c', 'src/cint1e.c', 'src/cint2e.c', 'src/cint_bas.c',
    'src/fblas.c', 'src/g1e.c', 'src/g2e.c', 'src/misc.c', 'src/optimizer.c',
    'src/fmt.c', 'src/rys_wheeler.c', 'src/eigh.c', 'src/rys_roots.c', 'src/find_roots.c',
    'src/cint2c2e.c', 'src/g2c2e.c', 'src/cint3c2e.c', 'src/g3c2e.c',
    'src/cint3c1e.c', 'src/g3c1e.c',  'src/breit.c',
    'src/cint1e_a.c', 'src/cint3c1e_a.c',
    'src/cint1e_grids.c', 'src/g1e_grids.c',
    'src/autocode/breit1.c', 'src/autocode/dkb.c', 'src/autocode/gaunt1.c',
    'src/autocode/grad1.c', 'src/autocode/grad2.c', 'src/autocode/hess.c',
    'src/autocode/int3c1e.c', 'src/autocode/int3c2e.c', 'src/autocode/intor1.c',
    'src/autocode/intor2.c', 'src/autocode/intor3.c', 'src/autocode/intor4.c',
    'src/autocode/deriv3.c', 'src/autocode/int1e_grids1.c',
    'src/autocode/deriv4.c', 'src/autocode/lresc.c'
)

includes = include_directories('src', 'include')

# build args
build_args = []

if cc.has_argument('-march=native')
    build_args += ['-march=native', '-fno-math-errno']
    if cc.get_id() == 'intel'
        build_args += ['-unroll-aggressive', '-ipo']
    else
        build_args += ['-funroll-loops', '-ftree-vectorize']
    endif
    message('added -march=native and loop unroling flags')
endif

quadmath_dep = cc.find_library('quadmath', required: false)
if quadmath_dep.found()
    project_dep += quadmath_dep
    build_args += '-DHAVE_QUADMATH_H'
    message('added quadmath')
endif

# deal with options:
if get_option('with_polynomial_fit')
    srcs += files('src/polyfits.c', 'src/sr_rys_polyfits.c')
    build_args += '-DWITH_POLYNOMIAL_FIT'
endif

if get_option('with_f12')
    srcs += files('src/g2e_f12.c', 'src/stg_roots.c', 'src/cint2e_f12.c')
    build_args += '-DWITH_F12'
endif

if get_option('with_4c1e')
    srcs += files('src/cint4c1e.c', 'src/g4c1e.c')
    message('Enabled WITH_4C1E. Note there are bugs in 4c1e integral functions')
endif

if get_option('pypzpx')
    build_args += '-DPYPZPX'
endif

if get_option('min_expcutoff') != ''
    build_args += '-DMIN_EXPCUTOFF=' + get_option('min_expcutoff')
endif

if get_option('keep_going')
    build_args += '-DKEEP_GOING'
endif

if get_option('with_fortran')
    build_args += '-DWITH_FORTRAN'
endif

if get_option('with_cint2_interface')
    build_args += '-DWITH_CINT2_INTERFACE'
endif

# library
libcint = library(
    meson.project_name(),
    srcs,
    install : not meson.is_subproject(),
    c_args : build_args,
    include_directories : includes,
    dependencies: project_dep
)

libcint_dep = declare_dependency(
    link_with: libcint,
    include_directories: includes,
    dependencies: project_dep,
)