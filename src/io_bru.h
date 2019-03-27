#ifndef IO_BRU_H
#define IO_BRU_H

#include "native.h"

#ifdef __cplusplus
extern "C" {
#endif

/* translates an CGNS file to the native mesh format */
void io_bru_write(MeshStruct* mesh);

#ifdef __cplusplus
}
#endif

#endif /* IO_BRU_H */
