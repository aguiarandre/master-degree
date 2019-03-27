#ifndef SOLUTION_H
#define SOLUTION_H

#include "config.h"
#include "native.h"

typedef struct
{
	int idx;
	char* name;
	char** componentName;
	char* var_type;
	char* num_type;
	int dim;
	int size;
	int location; /* 0 nodes, 1 elems, 2 faces, 3 bfaces*/
	char* data;
} Variable;

/* Main solution structure. */
typedef struct
{
	int num_vars;
	Variable* vars;
	const char* mesh_filename; /* the mesh this solution belongs to. */
} SolutionStruct;


void solution_read(const char* filename, MeshStruct* mesh, SolutionStruct* sol);

void solution2nodes(MeshStruct *mesh, SolutionStruct* sol, SolutionStruct* nodalSol);

#endif // SOLUTION_H
