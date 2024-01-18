# `stdlite`

An attempt at creating a lightweight library to perform sTD-DFT calculations (if possible, providing ways to use OpenMP/MPI).

## Note to self

+ There is a submodule, so one should also clone that:
  ```bash
  git clone https://github.com/pierre-24/stdlite --recurse-submodules
  ```
+ Testing framework is [Unity](https://github.com/ThrowTheSwitch/Unity).
+ Code documentation is generated through [doxide](https://doxide.org/), which itself generates a [mkdocs](https://squidfunk.github.io/mkdocs-material/) documentation. Better install that one in a *virtualenv*.