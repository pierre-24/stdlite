//
// Created by pierre on 1/22/24.
//

#ifndef STDLITE_UTILS_H
#define STDLITE_UTILS_H

#define STDL_FREE_IF_USED(a) if((a) != NULL) free(a)
#define STDL_RETURN_ON_ERROR(a) if((a) != STDL_ERR_OK) return (a)

#endif //STDLITE_UTILS_H
