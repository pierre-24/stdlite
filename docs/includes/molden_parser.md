By default, cartesian orbitals are assumed, unless `[5D]` or `[5D7F]` is found.
The ordering of the orbitals for a given angular momentum [follows the one of Gaussian](../fchk_parser/) (see [there](https://www.theochem.ru.nl/molden/molden_format.html)).

Program such as Dalton sort the MO according to symmetry, so the MO are re-sorted based on their energy, to correctly separate occupied and unoccupied MOs.

Dalton also name its symbol `X_nnn`, where `X` is a chemical symbol and `nnn` is a number. 
The latter is skipped.

!!! warning
    
    Mixing of cartesian and spherical functions is not supported. Thus, `[5D10F]` and `[7F]` will result in an error.