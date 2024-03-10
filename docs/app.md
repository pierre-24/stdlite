# Description of `stdlite_run` inputs

This page describe in details the inputs that `stdlite_run` takes.

## Preamble

In the following, "Keyword" refers to keywords that can be put in the TOML input, while "command line option" refers to option you can directly give to the program.
Note that a command line option, when it exists, have precedence over the TOML input file.

When a keyword type is "`str`/`float` (energy)", the energy is to be given in the format `"xxxYY"`, where `xxx` is a number, and `YY` is a unit.
Currently, 3 units are supported: `au` ([atomic units](https://en.wikipedia.org/wiki/Atomic_units)), `eV` ([electronvolts](https://en.wikipedia.org/wiki/Electronvolt)), and `nm` (nanometers).
If no unit are given, atomic units are assumed.
Valid inputs are, *e.g.*, `"8.5eV"`, `"1200nm"`, `"0.25au"`, `1.25`, etc.

## General 

!!! abstract "Title"

    **Type**: `str`
    **Keyword**: `title`
    **Default**: None

    Title of the calculation, in free format.
    It is there for documentation/archiving purposes, and will never be parsed.

## Context (`[context]`)

These are the keywords related to the creation of the context, *i.e.*, the selection of the CSFs and the creation of the $\mathbf A'$ and $\mathbf B'$ super-matrices. 

!!! abstract "Source"

    **Type**: `str` (path)
    **Keyword**: `source`
    **Command line option**: `--ctx_source`

    Path to a QM file containing at least a wavefunction and a corresponding basis set. 
    This is the only **mandatory** keyword.

!!! abstract "Source type"

    **Type**: `str`
    **Keyword**: `source_type`
    **Command line option**: `--ctx_source_type`
    **Default**: `"FCHK"`

    Type of source. Currently the possible values are:

    + `"MOLDEN"` (Molden file), 
    + `"FCHK"` (Gaussian FCHK), 
    + `"STDL_CTX"` (context from a previous calculation), and 
    + `"STDL_CTX_WB"` (context from a previous calculation, but only the wavefunction and basis are used).

    Note that `"STDL_CTX"` implies that exactly the same context, *i.e.*, the same CSFs and the same $\mathbf A'$ and $\mathbf B'$ super-matrices are used, so it totally skips the context creation.
    If you want to re-create a new context, use `"STDL_CTX_WB"`.

!!! abstract "Energy threshold"

    **Type**: `str`/`float` (energy)
    **Keyword**: `ethr`
    **Command line option**: `--ctx_ethr`
    **Default**: `"7eV"`

    Energy threshold for the [truncation of the CI space](theory.md#truncation-of-the-active-space).

!!! abstract "Perturbation energy threshold"

    **Type**: `str`/`float` (energy)
    **Keyword**: `e2thr`
    **Command line option**: `--ctx_e2thr`
    **Default**: `"1e-4au"`

    Perturbation energy threshold for the selection of secondary CSFs in the [truncation of the CI space](theory.md#truncation-of-the-active-space).

!!! abstract "Method for the creation of context"

    **Type**: `str`
    **Keyword**: `method`
    **Command line option**: `--ctx_method`
    **Default**: `"monopole"`

    Method for the calculation of the elements of the  $\mathbf A'$ and $\mathbf B'$ super-matrices. 
    Currently the possible values are:

    + `"monopole"` ([default sTD-DFT procedure](theory.md#std-dft-or-the-monopole-approximation), *i.e.*, the monopole approximation), 
    + `"monopole_direct"` (same as default, but the integrals are evaluated on demand).

    The `"xxx_direct"` version requires less memory.

!!! abstract "Tamm-Dancoff approximation"

    **Type**: `int`
    **Keyword**: `tda`
    **Command line option**: `--ctx_tda`
    **Default**: `1`

    If different from 0, use the [Tamm-Dancoff approximation](theory.md#application-to-dft-td-dft).

!!! abstract "Parameter for Coulomb integrals"

    **Type**: `float`
    **Keyword**: `gammaJ`
    **Command line option**: `--ctx_gammaJ`
    **Default**: `4.0`

    Parameters for the Coulomb integral, $\gamma_J$, in the [monopole approximation](theory.md#std-dft-or-the-monopole-approximation).

!!! abstract "Parameter for exchange integrals"

    **Type**: `float`
    **Keyword**: `gammaK`
    **Command line option**: `--ctx_gammaK`
    **Default**: `2.0`

    Parameters for the exchange integral, $\gamma_K$, in the [monopole approximation](theory.md#std-dft-or-the-monopole-approximation).

!!! abstract "Hartree-Fock exchange percentage"

    **Type**: `float`
    **Keyword**: `ax`
    **Command line option**: `--ctx_ax`
    **Default**: `0.5`

    Assuming a global hybrid, amount of Hartree-Fock exchange, between 0 and 1.

!!! abstract "Output"

    **Type**: `str` (path)
    **Keyword**: `output`
    **Command line option**: `--ctx_output`
    **Default**: `"context.h5"`

    Path to the place where the context will be saved, after it has been completed.

## Responses (`[responses]`)

These are the keywords related to the calculation of responses, their residues, and the related properties.

In the following, `wX` is a frequency, thus following the syntax for energy mentioned in the preamble.  
`opX` is an operator, which should be one of:

+ `"dipl"`: diple length operator;

Other operators will be added in the future.

!!! abstract "Linear reponses"

    **Type**: `list`
    **Keyword**: `linear`
    **Default**: `[]`

    List the linear reponses to compute. 
    Each linear response $\braket{\braket{\hat A;\hat B}}_\omega$ is to be given as: `{opA = "A", opB = "B", wB="w"}`.

    For example, the following input will compute the electric polarizability at 512 and 1064nm:

    ```toml
    [responses]
    linear = [
        {opA = 'dipl', opB = 'dipl', wB = '512nm'}
        {opA = 'dipl', opB = 'dipl', wB = '1064nm'}, 
    ]
    ```

!!! abstract "Single residue of the linear response"

    **Type**: `list`
    **Keyword**: `linear_sr`
    **Default**: `[]`

    List the ground-to-excited states properties to compute.
    Each request for $\braket{0|\hat A|m}$ is to be given as: `{op = "A", root = N}`.
    `N` is the number of excited states, `\ket{m}`, to consider.
    If `N` < 0, all possible excited states (*i.e.*, corresponding to the number of CSFs) are computed.

    For example, the following input will compute the transition dipole moments for the 15 first excited states:

    ```toml
    [responses]
    linear_sr = [
        {op = 'dipl', root = 15}
    ]
    ```

!!! abstract "Quadratic responses"

    **Type**: `list`
    **Keyword**: `quadratic`
    **Default**: `[]`

    List the quadratic reponses to compute. 
    Each quadratic response $\braket{\braket{\hat A;\hat B, \hat C}}_{\omega_1,\omega_2}$ is to be given as: `{opA = "A", opB = "B", opC = "C", wB="w1", wC="w2"}`.

    For example, the following input  will compute the SHG first hyperpolarizability at 1064nm:

    ```toml
    [responses]
    quadratic = [
        {opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = '1064nm', wC = '1064nm"}, 
    ]
    ```