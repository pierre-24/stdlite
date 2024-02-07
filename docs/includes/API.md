<!-- to be included in API/index.md -->

## Basic principles

+ Unless otherwise mentioned, everything inputs and outputs are (to be) given in [atomic units](https://en.wikipedia.org/wiki/Atomic_units) (of energy [hartree], of length [bohr], etc).

+ Every object, enum, or function that belongs to the public library API starts with `stdl_`.

+ Every array is stored in the usual [**row major**](https://en.wikipedia.org/wiki/Row-_and_column-major_order) form rather than the usual in fortran, column major.

+ The user is required delete all objects created within the library by using the provided deletor functions to avoid memory leaks.
  Furthermore, each function return an error code (see [error codes](./errors/stdl_error_code_)) or `STDL_ERR_OK` if everything went well.
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

+ Apart from the creator functions (`stdl_*_new()`) and cases where the size of the output is not known (such as in the selection of CSFs), the principle for outputs is "bring-your-own-space". 
  In other words, unless otherwise mentioned, outputs of functions (especially arrays) must have been allocated to their correct size.
  
## Specific topics