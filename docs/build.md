title: Build & install

## Preparation

To build `stdlite` from its sources, you'll need:

1. The [Meson build system](https://github.com/mesonbuild/meson), with a backend (generally [ninja](https://github.com/ninja-build/ninja)).
2. A linear algebra backend. Currently, only [openblas](https://www.openblas.net/)+lapack are supported.
3. The [HDF5 library](https://github.com/HDFGroup/hdf5) (and its headers), which is most probably available in your favorite distribution package manager.


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

Use the `meson` subcommands to set up and compile the project:

```bash
# go in folder
cd stdlite

# Setup the build.
# It downloads extra dependencies if required. 
meson setup _build

# Compile dependencies (if any) and stdlite.
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