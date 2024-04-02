title: Build & install

## Preparation

To build `stdlite` from its sources, you'll need:

1. A C compiler. Currently, GCC and [clang](https://clang.llvm.org/) are supported.
2. The [Meson build system](https://github.com/mesonbuild/meson), with a backend (generally [ninja](https://github.com/ninja-build/ninja)). This is probably available in you package manager.
3. A linear algebra backend which provides CBLAS and [LAPACKe](https://netlib.org/lapack/lapacke.html) interfaces for C. Currently, either [openblas64](https://www.openblas.net/)+Lapack(e) or [Intel MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library) are supported. To install MKL, which seems to provide better performances, see, *e.g.*, [this link](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html). Since [LAPACKe is not yet suported by flexiblas](https://github.com/mpimd-csc/flexiblas/issues/2), this option is not available at the moment. While not supported yet, Scalapack might also be an option, [see below](#custom-build).
4. The [HDF5 library](https://github.com/HDFGroup/hdf5) (and its headers), which is most probably available in your favorite distribution package manager.

??? note "Optional dependencies"

    Optionally, there are some (light-weight) dependencies that are installed by the project if not found, but that you can also install yourself if you prefer.

    For the library:
    
    + [`libcint`](https://github.com/sunqm/libcint/).

    For the tests:

    + [`Unity`](https://github.com/ThrowTheSwitch/Unity) (which is a test framework, not the game engine).

    For the application:

    + [`argtable3`](https://github.com/argtable/argtable3),
    + [`tomlc99`](https://github.com/cktan/tomlc99).

    Then, make sure that they are available in the corresponding `PATH`-like variables if needed.

Then, download the sources using `git`:

```bash
git clone https://github.com/pierre-24/stdlite.git
```

or directly on the [GitHub repository](https://github.com/pierre-24/stdlite), as a ZIP.

## Compilation

If you want to use `clang`, set the `CC` variable to do so:

```bash
export CC=clang
```

Use the `meson` subcommands to set up the project.

```bash
# go in folder
cd stdlite

# Setup the build.
# It downloads extra dependencies if required. 
meson setup _build --buildtype=release
```

The default build instruction use OpenBLAS+LAPACK (with `ILP64` enabled) and OpenMP. To change this, use:

```bash
# to use MKL (this will probably improve performances)
meson configure _build -Dla_backend=mkl

# to disable OpenMP (not recommended)
meson configure _build -Dopenmp=false
```

Then, compile the project:

```bash
meson compile -C _build
```

If you want, you can also run the test suite to check that everything is ok:

```bash
meson test -C _build
```

## Installation

To install the libraries and executables, use:

```bash
# (optional) if you want to change the default installation prefix:
meson configure _build --prefix=$HOME/.local

# install
meson install -C _build
```

Depending on the installation prefix and your user rights, root access might be required to perform this last action.

## Custom build

Successful build and run have also been achieved with:

+ OpenBLAS+ScaLAPACK: `-Dla_backend=custom -Dblas_vendor=openblas -Dlapack_vendor=scalpack`. Since [version 2.2.0](https://netlib.org/scalapack/scalapack-2.2.0.html), `ILP64` is supported, so `-Dc_args="-DInt=long -DSTDL_LA_INT=long"` should use 8-bytes integers.
+ OpenBLAS+LAPACK (no `ILP64`, so default integers): `-Dla_backend=openblas+lapack`

... But this has not been thoroughly tested yet.