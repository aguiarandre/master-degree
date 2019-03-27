#ifndef CONFIG_H
#define CONFIG_H


#define MAXENTITY 20   /* maximum number of entities recognised from IGES file */
#define MAXPOINT 500
#define MAXSPLINES 10
#define EPS 1E-12
#define RBF_TYPE 1
#define CHANGES 0.001 // 0.01, cylinder. 0.001, airfoil.
#define MESH_ORDER 2
#define MAX_ITERATION 200

/* Project information. */
#define PROJECT_DESCRIPTION "High Order Curved Mesh with RBF propagation"
#define PROJECT_VERSION "0.2.0"
#define PROJECT_AUTHOR "Andre Aguiar`"
#define PROJECT_AUTHOR_EMAIL "ufabc.andre@gmail.com"
#define BUILD_DATE "04/20/2017"
#define BUILD_TIME "08:49:17 PM BRT"

/* Build properties. */
#define PLATFORM_X86_64
#define HOST_LINUX
#define COMPILER_GNU
#define BUILD_DEBUG
#define LIB_STATIC

/* Third party libraries support. */
#define HAVE_HDF5
#define HAVE_CGNS
#define HAVE_OP2

/* Type definitions. */
#define REAL_IS_DOUBLE
typedef double Real;
typedef unsigned int uint;

/* program parameters. */
#define DIM 2
#define STRING_MAX_LENGHT 256
#define MAX_ELEM_TYPES 20
#if DIM == 2
//  #define MAX_ELEM_NODES 4
  #define MAX_ELEM_NBS 4  
  #if MESH_ORDER == 1 
	#define MAX_ELEM_NODES 6
	#define MAX_FACE_NODES 2
  #elif MESH_ORDER == 2
	#define MAX_ELEM_NODES 9
	#define MAX_FACE_NODES 3
  #elif MESH_ORDER == 3
	#define MAX_ELEM_NODES 16
	#define MAX_FACE_NODES 4
  #elif MESH_ORDER == 4
  	#define MAX_ELEM_NODES 25
	#define MAX_FACE_NODES 5
  #endif
#elif DIM == 3
  #define MAX_ELEM_NODES 8
  #define MAX_ELEM_NBS   6
  #define MAX_FACE_NODES 4
#endif


/* general headers. */
#include <stdbool.h>

#endif /* CONFIG_H */
