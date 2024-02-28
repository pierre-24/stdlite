#ifndef STDLITE_MOLDEN_PARSER_H
#define STDLITE_MOLDEN_PARSER_H

#include <stdlite/utils/base_parser.h>
#include <stdlite/utils/fchk_parser.h>

/**
 * Read section title.
 * Expects `[` (as the first character of a line!), ends after `]`.
 *
 * @param lx a valid lexer
 * @param title title to be read.
 * @return error code
 * @ingroup molden_parser
 */
int stdl_molden_parser_read_section_title(stdl_lexer* lx, char** title);

/**
 * Skip everything 'til next section. Expects `]`. Ends on `[`.
 *
 * @param lx a valid lexer
 * @return error code
 * @ingroup molden_parser
 */
int stdl_molden_parser_skip_section(stdl_lexer* lx);

/**
 * Read `[Atoms]` section content. Expects `]`. Ends on `[`.
 *
 * @param lx a valid lexer
 * @param[out] natm number of atoms, should be >0.
 * @param[out] atm `double[4*natm]` list of atoms with nuclear charge (0) and coordinates (1:3) to be created
 * @return error code
 * @ingroup molden_parser
 */
int stdl_molden_parser_read_atoms_section(stdl_lexer* lx, size_t* natm, double** atm);

/**
 * Read `[GTO]` section content. Expects `]`. Ends on `[`.
 *
 * @param lx a valid lexer
 * @param atm number of atoms, must be >0.
 * @param use_spherical assume spherical function instead of cartesian
 * @param[out] dt_ptr basis data to be created
 * @return error code
 * @ingroup molden_parser
 */
int stdl_molden_parser_read_gto_section(stdl_lexer *lx, size_t natm, int use_spherical, stdl_basis_data **dt_ptr);

/**
 * Read `[MO]` section content. Expects `]`. Ends on `[`.
 *
 * @param lx a valid lexer
 * @param nao number of ao, must be >0.
 * @param[out] nmo number of MO
 * @param e `double[nmo]` energies of each MO
 * @param C `double[nmo, nao]` coefficients for each MO
 * @return error code
 * @ingroup molden_parser
 */
int stdl_molden_parser_read_mo_section(stdl_lexer *lx, size_t nao, size_t *nmo, size_t *nocc, double **e, double **C);

/**
 * Extract a wavefunction (`stdl_wavefunction`) and a basis set (`stdl_basis`) from a MOLDEN file.
 * Expects the file to starts with `[Molden Format]`.
 *
 * @param lx a valid lexer
 * @param wf_ptr wavefunction to be created
 * @param bs_ptr basis set to be created
 * @return error code
 * @ingroup molden_parser
 *
 */
int stdl_molden_parser_extract(stdl_lexer* lx, stdl_wavefunction** wf_ptr, stdl_basis** bs_ptr);

#endif //STDLITE_MOLDEN_PARSER_H
