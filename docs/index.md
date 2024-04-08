title: Home

# `stdlite`

An attempt at creating a lightweight library to perform sTD-DFT calculations.

If you are not familiar with the sTD-DFT theory of Grimme and co-workers, it is explained [in the next page](theory.md).

Then, if you want to run a sTD-DFT calculation, you can [compile and install](build.md) `stdlite_run`, a standalone program that uses the [`stdlite` library](API/index.md) developed for this project.
A [tutorial](tutorial.md) as well as a [complete description of its inputs](app.md) is provided.

## What?

My Ph.D. was mainly dedicated to nonlinear optics (in particular the first and second hyperpolarizability) and during a collaboration with [Prof. M. de Wergifossse](https://uclouvain.be/en/research-institutes/imcn/most/prof-marc-de-wergifosse.html), I discovered the simplified approaches developed by [Prof S. Grimme](https://www.uni-bonn.de/en/research-and-teaching/research-profile/transdisciplinary-research-areas/tra-matter/members-directory/stefan-grimme).
Among other, we computed the first hyperpolarizability of a whole protein with almost 4000 atoms, in good agreement with experimental data (see [10.1021/acs.jpclett.1c02911](https://dx.doi.org/10.1021/acs.jpclett.1c02911)).

The sTDA/sTD-DFT approaches are mainly implement in [`stda`](https://github.com/grimme-lab/stda).
However, its integration with other quantum chemistry programs is not straightforward.
Following the example of [`tblite`](https://tblite.readthedocs.io/en/latest/), I thus decided to develop a standalone library.

Contributions are welcomed, as described on [this page](contributing.md).

## How

Rather than Fortran, this library is developed in C, with [libcint](https://github.com/sunqm/libcint) (used by [pycsf](https://pyscf.org/)) to evaluate the ops_integrals. 

## Who?

My name is [Pierre Beaujean](https://pierrebeaujean.net), and I'm a Ph.D. in quantum chemistry from the [University of Namur](https://unamur.be) (Belgium).
While now working on batteries, I still continue to look into nonlinear optics in collaboration with other people of my lab... And beyond.
It is also a good opportunity to improve my skills.
