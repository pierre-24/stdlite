<!-- to be included in API/index.md -->

## Basic principles

+ Unless otherwise mentioned, everything inputs and outputs are (to be) given in [atomic units](https://en.wikipedia.org/wiki/Atomic_units) (of energy [hartree], of length [bohr], etc).

+ Every object, enum, macro, or function that belongs to the public library API starts with `stdl_` (or `STDL_`).
  Everything that does not is internal, and subject to change without notice.

+ Every array is stored in the usual [**row major**](https://en.wikipedia.org/wiki/Row-_and_column-major_order) form rather than the usual in fortran, column major.

+ The user is required delete all objects created within the library by using the provided deletor functions to avoid memory leaks.
  Furthermore, each function return an error code (see [logging and errors](./logging/) or `STDL_ERR_OK` if everything went well.
  For example,

  ```c
  #include <stdlite/utils/base_parser.h>
  
  int an_integer, error;
  stdl_lexer* lexer = NULL;
  FILE* f = open("tmp.txt");
  
  if(f != NULL) {
    if(stdl_lexer_new(&lexer, f) == STDL_ERR_OK) {
      if(stdl_parser_get_integer(lx, &an_integer) == STDL_ERR_OK) {
        printf("read: %d\n", an_integer);
      }
    
      stdl_lexer_delete(lexer); // deletor
    }
    
    fclose(f);
  }
  ```

+ Apart from the creator functions (`stdl_*_new()`) and stuffs that comes from files (FCHK or MOLDEN), the principle for outputs is "bring-your-own-space". 
  In other words, unless otherwise mentioned, outputs of functions (especially arrays) must have been allocated to their correct size.

## General workflow

The general workflow to use `libstdlite` in order to compute a given property is the following:

![General workflow](../assets/activity_diagram_API.svg)


## Specific topics