# Tutorial

It is possible to compute responses and properties within the sTD-DFT framework by using `libstdlite` through the `stdlite_run` program.
This is obviously not restricted to `stdlite`, and also possible, *e.g.*, with the original implementation called [`stda`](https://github.com/grimme-lab/stda). 
All implementations should normally give the same answer.

!!! note "Preamble"

    To run a calculation, one needs a wavefunction from another quantum chemistry program as a *source*.
    Currently, two formats are supported: [FHCK](https://gaussian.com/interfacing/) and [MOLDEN](https://www.theochem.ru.nl/molden/molden_format.html).
    Utilities such as [`iodata`](https://github.com/theochem/iodata) provide ways to convert different outputs to one of these format (preferably MOLDEN).


## A first calculation: polarizabilities and beyond

The utility `stdlite_run` primarily read an input file written in the [TOML format](https://toml.io/).
For example, the following input set up a polarizability calculation, which is a linear response, given by $\alpha(-\omega;\omega) = -\braket{\braket{\mu;\mu}}_\omega$, followed by a first hyperpolarizability calculation, a quadratic response, $\beta(-\omega_\varsigma;\omega_1,\omega_2) = -\braket{\braket{\mu;\mu,\mu}}_{\omega_1,\omega_2}$ :

```toml
--8<-- "docs/assets/example.toml"
```

To perform the calculation, create a directory and download the inputs:

```bash
mkdir test
cd tests

# download source (water, computed at the HF/6-31G(d,f) level)
wget https://github.com/pierre-24/stdlite/raw/master/docs/assets/water.molden

# download TOML input
wget https://github.com/pierre-24/stdlite/raw/master/docs/assets/example.toml
```

and run the program:

```bash
stdlite_run example.toml > output.log
```

### Analysis of the output

As seen in the output file, a `stdlite` calculation is divided in different parts:

+ Read the user input and check that everything is ok;
+ Create the "context", *i.e.*, a) read the wavefuntion in the MOLDEN file, b) select the CSFs, and c) build the corresponding $\mathbf A'$ and $\mathbf B'$ super-matrices;
+ Compute the linear response and amplitude vectors;
+ Use said vectors to compute actual properties (here, the polarizability and hyperpolarizabilities).

Here, the result is:

```text 
--- Properties -----------------------------------------------------------------
Compute polarizability tensor >------< done
** alpha(-w;w), w=0.000000 (inf nm)
         x            y            z
x      1.62024      0.00000     -0.00000
y      0.00000      7.02043      0.00000
z     -0.00000      0.00000      3.69417
iso =        4.11162
aniso =      4.71843
Compute polarizability tensor >------< done
** alpha(-w;w), w=0.042823 (1064.00 nm)
         x            y            z
x      1.62570      0.00000     -0.00000
y      0.00000      7.05151      0.00000
z     -0.00000      0.00000      3.71297
iso =        4.13006
aniso =      4.74035
Compute polarizability tensor >------< done
** alpha(-w;w), w=0.085645 (532.00 nm)
         x            y            z
x      1.64285      0.00000     -0.00000
y      0.00000      7.14678      0.00000
z     -0.00000      0.00000      3.77115
iso =        4.18692
aniso =      4.80717
Compute first hyperpolarizability tensor >------------------< done
** beta(-w1-w2;w1,w2), w1=0.000000 (inf nm), w2=0.000000 (inf nm)
          x            y            z
xx      0.00000      0.00000      1.17469
xy      0.00000      0.00000      0.00000
xz      1.17469      0.00000     -0.00000
yx      0.00000      0.00000      0.00000
yy      0.00000      0.00000     12.28710
yz      0.00000     12.28710      0.00000
zx      1.17469      0.00000     -0.00000
zy      0.00000     12.28710      0.00000
zz     -0.00000      0.00000      9.31467
<B2ZZZ> =     75.54126
<B2ZXX> =     15.77636
BHRS    =      9.55603
DR      =      4.78826
Compute first hyperpolarizability tensor >------------------< done
** beta(-w1-w2;w1,w2), w1=0.042823 (1064.00 nm), w2=0.042823 (1064.00 nm)
          x            y            z
xx      0.00000      0.00000      1.25188
xy      0.00000      0.00000      0.00000
xz      1.25188      0.00000     -0.00000
yx      0.00000      0.00000      0.00000
yy      0.00000      0.00000     12.72614
yz      0.00000     12.72614      0.00000
zx      1.21095      0.00000     -0.00000
zy      0.00000     12.72475      0.00000
zz     -0.00000      0.00000      9.73270
<B2ZZZ> =     81.56844
<B2ZXX> =     16.89920
BHRS    =      9.92309
DR      =      4.82676
```

### To go further

It is also possible to pass arguments to the program, which takes precedence over the input file.
For example, to compute the same response with the Tamm-Dancoff approximation:

```bash
# this is equivalent to `tda = 1` in example.toml
stdlite_run example.toml --ctx_tda=1 > output.log
```

In this case, the results are similar with or without the Tamm-Dancoff approximation. This might not be the case for other systems.

You can also change the amount of molecular orbitals that are included:

```bash
# this is equivalent to `ethr = '25eV'` in example.toml
stdlite_run example.toml --ctx_ethr=25eV > output.log
```

In this case, 16 MOs and 28 CSFs are considered, rather than "only" 15 MOs and 23 CSFs when the threshold energy was of 20 eV.

## A second calculation: excitation energies

Excitation energies are obtained by solving the [Casida equation](theory.md#excitations), while ground to excited dipole moments are computed as the first residue of the linear response function, $\lim_{\omega\to\omega_m}\braket{\braket{\mu;\mu}}_\omega$.

The input file now looks like this (`linear_sr` request the single residue of linear response):

```toml
--8<-- "docs/assets/example2.toml"
```

This time, the Tamm-Dancoff approximation was requested.

Run the program:

```bash
# download this TOML input
wget https://github.com/pierre-24/stdlite/raw/master/docs/assets/example.toml

# run the program
stdlite_run example2.toml > output2.log
```

This time, once the wavefunction is extracted and $\mathbf A'$ is built, the amplitude vectors are computed.
 
### Analysis of the output

The result list each excitation energies with its corresponding transition dipole moment (in the dipole-length formalism), and the resulting oscillator strength (`fL`).

```text
--- Properties -----------------------------------------------------------------
Compute ground to excited transition dipole moments >< done
**    -------- Energy -------- ------ Transition dipole ---------
       (Eh)     (eV)    (nm)      X        Y        Z      fL   
   1  0.36143   9.835  126.06  0.23181  0.00000  0.00000 0.01295
   2  0.45262  12.317  100.66  0.00000 -0.00000  0.62895 0.11937
   3  0.45320  12.332  100.54 -0.00000 -0.00000 -0.00001 0.00000
   4  0.53684  14.608   84.87 -0.00000 -0.74083 -0.00000 0.19642
   5  0.60510  16.466   75.30 -0.00000 -0.98007 -0.00000 0.38748
   6  0.77077  20.974   59.11  0.00000 -0.00000  0.64883 0.21632
   7  1.08574  29.544   41.97 -0.00000 -0.00000 -0.00000 0.00000
   8  1.11337  30.296   40.92  0.11547  0.00000  0.00000 0.00990
   9  1.17348  31.932   38.83 -0.00000  0.73807  0.00000 0.42617
  10  1.21420  33.040   37.53  0.00000  0.00000 -0.16143 0.02110
```

### To go further

You can totally combine linear/quadratic response calculations and excitations with the following input:

```toml
# compute the static polarizability
linear = [
    {opA = 'dipl', opB = 'dipl', wB = 0}, 
]

# compute the static hyperpolarizability
quadratic = [
    {opA = 'dipl', opB = 'dipl', opC = 'dipl', wB = 0, wC = 0},
]

# compute the 10 first excitation energies
linear_sr = [
    {opA = 'dipl', nroots = -1},
]
```

Note that `nroots=-1` request the computation of all possible excitations (in this case, 23), which might take longer.