# Using the application

It is possible to compute responses and properties through `stdlite_run`.

!!! note "Preamble"

    To run a calculation, one needs a wavefunction from another quantum chemistry program as a *source*.
    Currently, two formats are supported: [FHCK](https://gaussian.com/interfacing/) and [MOLDEN](https://www.theochem.ru.nl/molden/molden_format.html).
    Utilities such as [`iodata`](https://github.com/theochem/iodata) provide ways to convert different output to one of these format (preferably MOLDEN).

## Example

The utility `stdlite_run` primarily read an input file written in the [TOML format](https://toml.io/).
For example, the following input performs a polarizability calculation, which is a linear response:

```toml
--8<-- "docs/assets/example.toml"
```

To perform the calculation, create a directory and download the inputs:

```bash
mkdir test
cd tests

# download source (water, computed at the HF/STO-3G level)
wget https://github.com/pierre-24/stdlite/raw/master/docs/assets/water.molden

# download TOML input
wget https://github.com/pierre-24/stdlite/raw/master/docs/assets/example.toml
```

and run the program:

```bash
stdlite_run example.toml > output.log
```

(...)

It is also possible to pass arguments to the program, which takes precedence over the input file (see below). 
For example, to compute the same response without the Tamm-Dancoff approximation:

```bash
stdlite_run example.toml --ctx_tda=0 > output.log
```

## Description of the input

(...)