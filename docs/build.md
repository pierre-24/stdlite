title: Build & install

## Preparations

To build `stdlite` from its sources, you need:

1. The [Meson build system](https://github.com/mesonbuild/meson), with a backend (generally [ninja](https://github.com/ninja-build/ninja)).
2. A linear algebra backend. Currently, only [openblas](https://www.openblas.net/) is supported.


!!! note

    Optionally, there are some dependencies that are installed by the project if not found, but that you can also install yourself if you prefer:
    
    + [`libcint`](https://github.com/sunqm/libcint/),
    + [`Unity`](https://github.com/ThrowTheSwitch/Unity) (which is a test framework, not the game engine).

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