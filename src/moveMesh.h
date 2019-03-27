
#ifndef MOVE_MESH_H
#define MOVE_MESH_H

#include "rbf_functions.h"

// functions to effectivelly move the mesh for a given prescribed movement:
void moveMesh( MeshStruct * mesh );
void getBoundaryNodes( MeshStruct * mesh );

#endif
