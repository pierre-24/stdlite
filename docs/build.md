title: Build & install `stdlite`

To build `stdlite` from its sources, you need:

1. The [Meson build system](https://github.com/mesonbuild/meson), with a backend (generally [ninja](https://github.com/ninja-build/ninja)).
2. A linear algebra backend. Currently, only [openblas](https://www.openblas.net/) is supported.

Then, download the sources:

```bash
# download
git clone https://github.com/pierre-24/stdlite.git

# go in folder
cd stdlite
```

Use the `meson` subcommands to set up and compile the project:

```bash
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

To install the libraries and executables, use:

```bash
# (optional) if you want to change the default installation prefix:
meson configure _build --prefix=$HOME/.local

# install
meson install -C _build
```

Depending on the installation prefix and your user rights, root access might be required to perform this last action.