While this parser contains a few safeguards, it works better if the FCHK is correctly formatted.

!!! warning

    The logical type, `L`, is not implemented. However, it does not seem to be used in production FCHK.
    Scalar `C` is also not implemented nor used.

The ordering of the orbitals for a given angular momentum is specific to Gaussian, but `stdlite` expects inputs and produce outputs following the order of `libcint`.
The following transposition table is used for LCAO coefficients.

For cartesian orbitals (`6d`, `10f`, `15g`, etc)

| libcint ([source](https://github.com/sunqm/libcint/blob/master/doc/program_ref.txt))                                                                                                                                                                | Gaussian ([source](https://gaussian.com/interfacing/))                                                                                                                                                                                               |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| $s$                                                                                                                                                                                                                                                 | $s$                                                                                                                                                                                                                                                  |
| $p_x$, $p_y$, $p_z$                                                                                                                                                                                                                                 | $p_x$, $p_y$, $p_z$                                                                                                                                                                                                                                  |
| $d_{xx}$ (0), $d_{xy}$ (1), $d_{xz}$ (2), $d_{yy}$ (3), $d_{yz}$ (4), $d_{zz}$ (5)                                                                                                                                                                  | $d_{xx}$ (0), $d_{yy}$ (3), $d_{zz}$ (5), $d_{xy}$ (1), $d_{xz}$ (2), $d_{yz}$ (4)                                                                                                                                                                   |
| $f_{xxx}$ (0), $f_{xxy}$ (1), $f_{xxz}$ (2), $f_{xyy}$ (3), $f_{xyz}$ (4), $f_{xzz}$ (5), $f_{yyy}$ (6), $f_{yyz}$ (7), $f_{yzz}$ (8), $f_{zzz}$ (9)                                                                                                | $f_{xxx}$ (0), $f_{yyy}$ (6), $f_{zzz}$ (9), $f_{xyy}$ (3), $f_{xxy}$ (1), $f_{xxz}$ (2), $f_{xzz}$ (5), $f_{yzz}$ (8), $f_{yyz}$ (7), $f_{xyz}$ (4)                                                                                                 |
| $g_{xxxx}$ (0), $g_{xxxy}$ (1), $g_{xxxz}$ (2), $g_{xxyy}$ (3), $g_{xxyz}$ (4), $g_{xxzz}$ (5), $g_{yyxz}$ (6), $g_{yyyx}$ (7), $g_{yyyy}$ (8), $g_{yyyz}$ (9), $g_{yyzz}$ (10), $g_{zzxy}$ (11), $g_{zzzx}$ (12), $g_{zzzy}$ (13), $g_{zzzz}$ (14) | $g_{xxxx}$ (0), $g_{yyyy}$ (8), $g_{zzzz}$ (14), $g_{xxxy}$ (1), $g_{xxxz}$ (2), $g_{yyyx}$ (7), $g_{yyyz}$ (9), $g_{zzzx}$ (12), $g_{zzzy}$ (13), $g_{xxyy}$ (3), $g_{xxzz}$ (5), $g_{yyzz}$ (10), $g_{xxyz}$ (4), $g_{yyxz}$ (6), $g_{zzxy}$ (11)  |

For spherical orbitals (`5d`, `7f`, etc., see [there for definition](https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics)):

| libcint ([source](https://github.com/sunqm/libcint/blob/master/doc/program_ref.txt))                          | Gaussian ([source](https://gaussian.com/interfacing/))                                                        |
|---------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| $d_{-2}$ (0), $d_{-1}$ (1), $d_0$ (2), $d_1$ (3), $d_2$ (4)                                                   | $d_0$ (2), $d_1$ (3), $d_{-1}$ (1), $d_2$ (4), $d_{-2}$ (0)                                                   |
| $f_{-3}$ (0), $f_{-2}$ (1), $f_{-1}$ (2), $f_0$ (3), $f_1$ (4), $f_2$ (5), $f_3$ (6)                          | $f_0$ (3), $f_{1}$ (4), $f_{-1}$ (2), $f_{2}$ (5), $f_{-2}$ (1), $f_{3}$ (6), $f_{-3}$ (0)                    |
| $g_{-4}$ (0), $g_{-3}$ (1), $g_{-2}$ (2), $g_{-1}$ (3), $g_0$ (4), $g_1$ (5), $g_2$ (6), $g_3$ (7), $g_4$ (8) | $g_0$ (4), $g_1$ (5), $g_{-1}$ (3), $g_2$ (6), $g_{-2}$ (2), $g_3$ (7), $g_{-3}$ (1), $g_4$ (8), $g_{-4}$ (0) |


!!! warning

    For the moment, only functions up to $g$ are thus handled.

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