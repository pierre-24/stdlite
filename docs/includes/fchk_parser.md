While this parser contains a few safeguards, it works better if the FCHK is correctly formatted.

!!! warning

    The logical type, `L`, is not implemented. However, it does not seem to be used in production FCHK.
    Scalar `C` is also not implemented nor used.

The ordering of the orbitals for a given angular momentum is specific to Gaussian, but `stdlite` expects inputs and produce outputs following the order of `libcint`.
The following transposition table is used for LCAO coefficients:

| libcint ([source](https://github.com/sunqm/libcint/blob/master/doc/program_ref.txt))                                                                 | Gaussian ([source](https://gaussian.com/interfacing/))                                                                                               |
|------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------|
| $s$                                                                                                                                                  | $s$                                                                                                                                                  |
| $p_x$, $p_y$, $p_z$                                                                                                                                  | $p_x$, $p_y$, $p_z$                                                                                                                                  |
| $d_{xx}$ (0), $d_{xy}$ (1), $d_{xz}$ (2), $d_{yy}$ (3), $d_{yz}$ (4), $d_{zz}$ (5)                                                                   | $d_{xx}$ (0), $d_{yy}$ (3), $d_{zz}$ (5), $d_{xy}$ (1), $d_{xz}$ (2), $d_{yz}$ (4)                                                                   |
| $d_{xy}$ (0), $d_{yz}$ (1), $d_{z^2}$ (2), $d_{xz}$ (3), $d_{x^2-y^2}$ (4)                                                                           | $d_{z^2}$ (2), $d_{xz}$ (3), $d_{yz}$ (1), $d_{x^2-y^2}$ (4), $d_{xy}$ (0)                                                                           |
| $f_{xxx}$ (0), $f_{xxy}$ (1), $f_{xxz}$ (2), $f_{xyy}$ (3), $f_{xyz}$ (4), $f_{xzz}$ (5), $f_{yyy}$ (6), $f_{yyz}$ (7), $f_{yzz}$ (8), $f_{zzz}$ (9) | $f_{xxx}$ (0), $f_{yyy}$ (6), $f_{zzz}$ (9), $f_{xyy}$ (3), $f_{xxy}$ (1), $f_{xxz}$ (2), $f_{yzz}$ (8), $f_{xzz}$ (5), $f_{yyz}$ (7), $f_{xyz}$ (4) |
| $f_{y(3x^2-y^2)}$ (0), $f_{xyz}$ (1), $f_{yz^2}$ (2), $f_{z^3}$ (3), $f_{xz^2}$ (4), $f_{z(x^2-y^2)}$ (5), $f_{x(x^2-3y^2)}$ (6)                     | $f_{z^3}$ (3), $f_{xz^2}$ (4), $f_{yz^2}$ (2), $f_{z(x^2-y^2)}$ (5), $f_{xyz}$ (1), $f_{x(x^2-3y^2)}$ (6), $f_{y(3x^2-y^2)}$ (0)                     |


!!! warning

    For the moment, only functions up to $f$ are thus handled.

In order to parse a FCHK, use something of the form:

```c
char* name = NULL;
char type;
int error, is_scalar;

// 1. Open file, create lexer
FILE* f = fopen("file.fchk", "r");
stdl_lexer* lx = stdl_lexer_new(f);

// 2. Skip intro
stdl_fchk_parser_skip_intro(lx);

// 3. Read sections
while (lx->current_tk_type != STDL_TK_EOF) {
  // read section info, giving:
  // a) the name of the section, b) its type, and c) if it is a scalar
  error = stdl_fchk_parser_get_section_info(lx, &name, &type, &is_scalar);
  
  if(error == STDL_ERR_OK) {
    if(strcmp("an interesting section", name) == 0) {
      /* Read the content of this section by using either:
       * - `stdl_fchk_parser_get_scalar_*()`, or
       * - `stdl_fchk_parser_get_vector_*()`.
       * Don't forget to free the result after use.
       */
    } else if(/* ... */) {
      /* ... */
    } else {
      // not interesting, skip section
      stdl_fchk_parser_skip_section(lx, type, is_scalar);
    }
    free(name);
  }
}

stdl_lexer_delete(lx);
fclose(f);
```