#ifndef IO_CGNS_H
#define IO_CGNS_H

#include <cgnslib.h>

#include "native.h"
#include "solution.h"

#ifdef __cplusplus
extern "C" {
#endif

/* translates an CGNS file to the native mesh format */
void io_cgns_read(const char* filename, MeshStruct* mesh);

/* Fetch the mesh vertices coordinates directly to the native
 * mesh data structure. Allocates mesh struct memory.
 */
void io_cgns_read_nodes(int cgfile, int baseid, int zoneid,
                        cgsize_t* size, MeshStruct* mesh);

/* Fetch mesh elements connectivity directly to the native
 * mesh data structure. Allocates mesh struct memory.
 * The the native mesh zones is also defined.
 */
void io_cgns_read_elems_connectivity(int cgfile, int baseid, int zoneid,
                                     int nsections, cgsize_t* size, MeshStruct* mesh);


/* translates the native mesh format to an CGNS file */
void io_cgns_write(const char* filename, MeshStruct* mesh, SolutionStruct* sol);

#ifdef __cplusplus
}
#endif

#endif /* IO_CGNS_H */
