#ifndef STDLITE_H
#define STDLITE_H

/// Free `a` if not `NULL`
#define STDL_FREE_IF_USED(a) if((a) != NULL) free(a)

/// Return if `a` is not `STDL_ERR_OK`
#define STDL_RETURN_ON_ERROR(a) if((a) != STDL_ERR_OK) return (a)

#endif //STDLITE_H
