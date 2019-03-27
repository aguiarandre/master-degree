/* This header holds the mesh structure in its native
 * format, used by the solver
 */

#ifndef NATIVE_H
#define NATIVE_H

#include <hdf5.h>
#include <hdf5_hl.h>

#include "config.h"

/* define uninitialized constant value. */
static int const UNINIT = -1;

/* --------------------- Mesh element types --------------------- */

/* Mesh element type.
 * description: TYPE_NODES, where
 * TYPE: element geometric type,
 * NODES: how many nodes form this element
 * Max size of supported types is MAX_ELEM_TYPES.
 */
typedef enum
{
  pNONE = 0,
  pNODE,
  pBAR_2,	// Bar element - 1st order
  pTRI_3,	// Triangular - 1st order
  pTRI_6,	// Triangular - 2nd order
  pQUAD_4,	// Quadrilateral - 1st order
  pQUAD_9,	// Quadrilateral - 2nd order
  pTETRA_4,	// Tetrahedral - 1st order
  pTETRA_10,	// Tetrahedral - 2nd order
  pPYRA_5,	// Pyramidal - 1st order
  pPENTA_6,	// Pentahedral (extruded triangle) - 1st order
  pHEXA_8, 	// Hexahedron - 1st order
  pMIXED	// More than 1 type of element
} ElemType;

/* For a given mixed element zone, this struct holds a
 * particular element type start index and size in the
 * ordered map of elements.
 */
typedef struct { // so if you have mixed type of elements, this is helpful, otherwise start = 1, num_elems = n_total
		
  ElemType elem_type;  // for this specific type of element
  uint start;          // define an unsigned integer - index of starting element of this element type
  uint num_elems;      // number of elements of this type
} ElemTypeList;




/* --------------------- Mesh mappings --------------------- */

/* struct to map faces to nodes and elements. */
typedef struct  
//given the number of the *face*, obtain its vertex index and its connectivity table. 
{
  uint* nodes;		//pointer to index of the vertex of the ith element of the table (p)  (the index or the x,y,z position?)
  uint  elem_left,	//declaring integers - index of left and right elements, connectivity table (t) ?
        elem_right;
} MapFaces;

/* struct to map boundary faces to nodes and elements. */
typedef struct

{
  uint* nodes;		//pointer to location of the x,y,z points of this element
  uint  elem_left;	//index of the one neighbour this element has
} MapBFaces;

/* struct to map elements to nodes connectivity */
typedef struct
// given the *element*, obtain its x,y,z position; in
{
  uint* nodes;  
  uint* nbs;	/* element neighbors */
  uint* faces;  /* face Number that form this element */
  bool* isFaceBoundary;
  
} MapElems;

/* struct to create mesh edges in 3D, mainly for 
 * visualization. */


typedef struct
{
  uint n1, n2;	
} WireMesh;

/* --------------------- Mesh data --------------------- */

/* mesh data sizes. */

// Total of mesh elements

typedef struct
{
  uint nodes;				// number of vertices
  uint elems;				// number of elements
  uint faces;				// number of internal faces (or total?)
  uint bfaces;				// number of boundary faces
  uint zones;				// number of different zones
  uint elem_type_count[MAX_ELEM_TYPES]; /* counter of element types in the mesh. */
  uint splines;
  
} MeshSize;

/* Mesh zone information. */

// Characteristics of Zone

typedef struct
{											// io_native_read.c
  uint          idx; /* zone index */							// 1,2,3 .... 
  char*         name; /* name of this zone */ 						// pointer to the name of zone
  short         dim; /* zone index physical dimension */				
  bool          is_boundary; /* determines if this is a domain boundary			// important	
                                zone (true) or an internal region (false) */
  ElemType      elem_type; /* zone elements type */					// type of element of this zone
  bool          homogeneous; /* true if zone has only one type of elements. */		// if elem_type != pMixed, then 1
  short         num_types; /* number of different element types in this zone. */	// if homoegenous != 1, num_types = sizeof(elem_type) ?
  ElemTypeList* type_list; /* stores the index of sorted element types. */		// pointer to the types of elements
  uint          start; /* zone elements counter index start */				// starting index
  uint          num_elems; /* number of elements in this zone */			// ending index - starting index
  int 		    bc_id; /* boundary condition identifier */				// type of boundary condition, e.g. 1 = wall, 2 = farfield, etc
  void*         bc_data; /* pointer to a boundary data struct */			// ??
  

} MeshZone;


typedef struct
{
	int nodeNumber;
	double x;
	double y;
	double z; 
	
} pFixedNodeStruct;

typedef struct
{
	int nodeNumber;
	double x;
	double y;
	double z; 
	
} fixedNodeStruct;

typedef struct
{
	int nodeNumber;
	double x;
	double y;
	double z; 
	
} movingNodeStruct;

typedef struct 
{
	int * isSet;
	
	pFixedNodeStruct * fixedNodesBefore;
	fixedNodeStruct * fixedNodes;
	movingNodeStruct * movingNodes;
	
	int nFixedNodes;
	int nMovingNodes;
	
	double ** C;
	double ** A;
	double ** cInverse;
	double ** H;

} rbfStruct; 
/* Main mesh structure. */

// Characteristics of Mesh

typedef struct
{
  int idx; /* mesh index for the current simulation */					// 
  int dim; /* mesh physical dimension */						// 
  const char* filename; /* mesh file in which we read the mesh. */			// ok
  ElemType elem_type;									// structure giving the type of elements of the mesh
  MeshSize sizes; /* mesh size data */							// ??
  MeshZone* zones; /* pointer to an array of zones */					// leads to struct above
  Real** vertices; /* pointer to an array of nodes coordinates */			// pointer that leads to (x,y,z) of a vertice (node)
  MapFaces*  mFaces;									// pointer to structure containing characteristics of internal faces
  MapBFaces* mBFaces;									// pointer to structure containing characteristics of boundary faces
  MapElems*  mElems;									// pointer to structure containing characteristics of elements
  rbfStruct* rbf;
  double* lookUpArray;
  double        time;
  double deltat;
  
} MeshStruct;


/* --------------------- Function prototypes --------------------- */

#ifdef __cplusplus
extern "C" {
#endif

void help(void);

void help_datproc(void);


/* Mesh initialization and destruction functions */

void mesh_init(MeshStruct* mesh);

void mesh_finish(MeshStruct* mesh);

void mesh_zone_init(MeshZone* zone);

void mesh_zone_finish(MeshZone* zone);

void mesh_report(MeshStruct* mesh);


/* Mesh derived data functions */

uint mesh_compute_faces_number(MeshStruct* mesh);

void get_bound(MeshStruct* mesh);

uint get_hash_key( short dim, int ipr, uint hashspace,
                   uint* a, uint* b, uint* c, uint* d );

void get_elem_faces( MeshStruct* mesh, uint zoneid, uint elem,
                     short* num_faces, uint** nodes);


/* Create the faces to nodes map and determine face neighbors,
 * even for boundary ones. Note that it is assumed that the
 * BFaces to nodes map is already available in the mesh structure.
 * Typical mesh generators will provide such information.
 */
void mesh_create_faces(MeshStruct* mesh);

void get_elem_nodes(uint elem, uint** nodes);

void mesh_elems_by_type(ElemType type, MeshStruct* mesh, uint* data);

void mesh_set_bcs(MeshStruct* mesh);


/* IO functions - native */

void mesh_write_meta_info(const char* filename, MeshStruct* mesh);

void mesh_write_data_native(const char* filename, MeshStruct* mesh);

void mesh_read_meta_info(const char* filename, MeshStruct* mesh);

void mesh_read_data_native(const char* filename, MeshStruct* mesh);


/* Helper hdf5 functions */

void hdf5_read_scalar(hid_t file, char const *name, char const *type, void* data);

void hdf5_read_data_dim(hid_t file, char const *name, uint* sizes);

void hdf5_read_data(hid_t file, char const *name, char const *type, void* data);

void hdf5_read_vlstring(hid_t file, const char* name, uint* sizes, void* data);

void hdf5_write_scalar(hid_t file, const char* name, const char* type,
                       const void* data);

void hdf5_write_data(hid_t file, const char* name, const char* type,
                     uint* sizes, bool attr_size_bytes, const void* data);

void hdf5_write_vlstring(hid_t file, const char* name, uint* sizes, const void* data);

/* Helper utility functions */

uint utils_next_prime(uint offset);

int compare_uint(const void* val1, const void* val2 );

uint unique (uint* first, uint* last);

uint binary_search( uint* array, uint size, uint value );


/* IO functions - GUI */

void mesh_export_obj(MeshStruct* mesh);

uint hash_3D_wireframe(WireMesh* edges, uint elems, uint size);

/* High-Order Mesh generation functions */

void igs_read(int argc, char * igesFile, MeshStruct* mesh);

void propagateDeformation( MeshStruct* mesh );



#ifdef __cplusplus
}
#endif

#endif /* NATIVE_H */
