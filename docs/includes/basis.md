Follows the structure of [`libcint`](https://github.com/sunqm/libcint/blob/master/doc/program_ref.txt), so that it can be used to compute extra integrals.

!!! note
  
    Most of the documentation (and all examples) of libcint actually refers to the old (version < 3.0) API.
    The new API actually includes two extra parameters [described there](https://github.com/sunqm/libcint/blob/master/doc/program_ref.tex#L228).
    Furthermore, according to [this](https://github.com/sunqm/libcint/issues/76), the first 20 values of `env` are reserved.

!!! warning

    Per [this comment](https://github.com/sunqm/libcint/issues/22#issuecomment-333268020), cartesian Gaussians are not (correctly) normalized by libcint :(