#include "native.h"

#include <stdlib.h> // malloc
#include <string.h>


/* Read meta information about the mesh. Recover data for
 * a MeshStruct struc.
*/
void mesh_read_meta_info(const char* filename, MeshStruct* mesh)
{
		
	hid_t file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );

	// Read mesh meta data - Size information
	{
		typedef struct
		{
			uint nodes;
			uint elems;
			uint faces;
			uint bfaces;
			uint zones;
		} MetaSize;
		
		MetaSize* data = (MetaSize*)malloc( sizeof(MetaSize) );
	
		size_t dst_size = sizeof( MetaSize );
		
		size_t dst_offset[5] = { HOFFSET( MetaSize, nodes ),
								 HOFFSET( MetaSize, elems ),
								 HOFFSET( MetaSize, faces ),
								 HOFFSET( MetaSize, bfaces ),
								 HOFFSET( MetaSize, zones )};
		
		size_t dst_sizes[5] = { sizeof( data[0].nodes),
								sizeof( data[0].elems),
								sizeof( data[0].faces),
								sizeof( data[0].bfaces),
								sizeof( data[0].zones)};
		
		H5TBread_table( file, "meta_sizes", dst_size, dst_offset, dst_sizes, data );
		
		mesh->sizes.nodes  = data[0].nodes;
		mesh->sizes.elems  = data[0].elems;
		mesh->sizes.faces  = data[0].faces;
		mesh->sizes.bfaces = data[0].bfaces;
		mesh->sizes.zones  = data[0].zones;
		
		free( data );
	}
	
	
	// Read mesh meta data - Zones information
	{
		typedef struct
		{
			uint idx;
			char name[STRING_MAX_LENGHT];
			int  dim;
			int  is_boundary;
			char elem_type[STRING_MAX_LENGHT];
			int homogeneous;
			int num_types;
			int start;
			int num_elems;
			int	bc_id;
		} MetaZone;
	
		MetaZone* data = (MetaZone*)malloc( mesh->sizes.zones*sizeof(MetaZone) );
		
		size_t dst_size = sizeof( MetaZone );
		
		size_t dst_offset[10] = { HOFFSET( MetaZone, idx ),
								HOFFSET( MetaZone, name ),
								HOFFSET( MetaZone, dim ),
								HOFFSET( MetaZone, is_boundary ),
								HOFFSET( MetaZone, elem_type ),
							  HOFFSET( MetaZone, homogeneous ),
							  HOFFSET( MetaZone, num_types ),
							  HOFFSET( MetaZone, start ),
							  HOFFSET( MetaZone, num_elems ),
							  HOFFSET( MetaZone, bc_id )};

		size_t dst_sizes[10] = { sizeof( data[0].idx ),
								sizeof( data[0].name ),
								sizeof( data[0].dim ),
								sizeof( data[0].is_boundary ),
								sizeof( data[0].elem_type ),
							  sizeof( data[0].homogeneous ),
							  sizeof( data[0].num_types ),
							  sizeof( data[0].start ),
							  sizeof( data[0].num_elems ),
							  sizeof( data[0].bc_id )};
							  
		H5TBread_table( file, "meta_zones", dst_size, dst_offset, dst_sizes, data );
		
		// populate the mesh zone struct
		
		mesh->zones = (MeshZone*)malloc( mesh->sizes.zones*sizeof(MeshZone) );
		
		for ( int z = 0; z < mesh->sizes.zones; z++ )
		{
			mesh_zone_init( &mesh->zones[z] );

			mesh->zones[z].idx  = data[z].idx;

			mesh->zones[z].name = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
			strcpy(mesh->zones[z].name, data[z].name);

			mesh->zones[z].dim  = data[z].dim;

			if ( data[z].is_boundary == 0 ) mesh->zones[z].is_boundary = false;
			if ( data[z].is_boundary == 1 ) mesh->zones[z].is_boundary = true;

			if ( strcmp(data[z].elem_type,"BAR_2") == 0 ) {
				mesh->zones[z].elem_type = pBAR_2;
			} else if ( strcmp(data[z].elem_type,"TRI_3") == 0 ) {
				mesh->zones[z].elem_type = pTRI_3;
			} else if ( strcmp(data[z].elem_type,"TRI_6") == 0 ) {
				mesh->zones[z].elem_type = pTRI_6;
			} else if ( strcmp(data[z].elem_type,"QUAD_4") == 0 ) {
				mesh->zones[z].elem_type = pQUAD_4;
			} else if ( strcmp(data[z].elem_type,"QUAD_9") == 0 ) {
				mesh->zones[z].elem_type = pQUAD_9;
			} else if ( strcmp(data[z].elem_type,"TETRA_4") == 0 ) {
				mesh->zones[z].elem_type = pTETRA_4;
			} else if ( strcmp(data[z].elem_type,"TETRA_10") == 0 ) {
				mesh->zones[z].elem_type = pTETRA_10;
			} else if ( strcmp(data[z].elem_type,"PYRA_5") == 0 ) {
				mesh->zones[z].elem_type = pPYRA_5;
			} else if ( strcmp(data[z].elem_type,"PENTA_6") == 0 ) {
				mesh->zones[z].elem_type = pPENTA_6;
			} else if ( strcmp(data[z].elem_type,"HEXA_8") == 0 ) {
				mesh->zones[z].elem_type = pHEXA_8;
			} else if ( strcmp(data[z].elem_type,"MIXED") == 0 ) {
				mesh->zones[z].elem_type = pMIXED;
			} else {
				mesh->zones[z].elem_type = pNONE;
			}
			
			if ( data[z].homogeneous == 0 ) mesh->zones[z].homogeneous = false;
			if ( data[z].homogeneous == 1 ) mesh->zones[z].homogeneous = true;
			
			mesh->zones[z].num_types = data[z].num_types;
			mesh->zones[z].start 		 = data[z].start;
			mesh->zones[z].num_elems = data[z].num_elems;
			mesh->zones[z].bc_id 		 = data[z].bc_id;
		}
		free( data );
	}
	
	
	// Read mesh meta data - Element type information
	{
		typedef struct
		{
			char elem_type_name[STRING_MAX_LENGHT];
			int size;
			int nodes;
		} MetaElem;
		
		int num_valid_elems = 10;
		MetaElem* data = (MetaElem*)malloc( num_valid_elems*sizeof(MetaElem) );
		
		size_t dst_size = sizeof( MetaElem );
		
		size_t dst_offset[3] = { HOFFSET( MetaElem, elem_type_name ),
								HOFFSET( MetaElem, size ),
								HOFFSET( MetaElem, nodes )};

		size_t dst_sizes[3] = { sizeof( data[0].elem_type_name ),
								sizeof( data[0].size ),
								sizeof( data[0].nodes )};
							  
		H5TBread_table( file, "meta_elems", dst_size, dst_offset, dst_sizes, data );
		
		for ( int i = 0; i < num_valid_elems; i++ )
		{
			if ( strcmp(data[i].elem_type_name,"BAR_2") == 0 ) {
				mesh->sizes.elem_type_count[2] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"TRI_3") == 0 ) {
				mesh->sizes.elem_type_count[3] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"TRI_6") == 0 ) {
				mesh->sizes.elem_type_count[4] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"QUAD_4") == 0 ) {
				mesh->sizes.elem_type_count[5] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"QUAD_9") == 0 ) {
				mesh->sizes.elem_type_count[6] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"TETRA_4") == 0 ) {
				mesh->sizes.elem_type_count[7] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"TETRA_10") == 0 ) {
				mesh->sizes.elem_type_count[8] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"PYRA_5") == 0 ) {
				mesh->sizes.elem_type_count[9] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"PENTA_6") == 0 ) {
				mesh->sizes.elem_type_count[10] = data[i].size;
			} else if ( strcmp(data[i].elem_type_name,"HEXA_8") == 0 ) {
				mesh->sizes.elem_type_count[11] = data[i].size;
			}
		}
		
		free( data );
	}
	
	// close file
	H5Fclose( file );
}


void mesh_read_data_native(const char* filename, MeshStruct* mesh)
{
  printf("\nReading native mesh file.\n");
  fflush(stdout);

  hid_t file;

  // Open the HDF5 file in read-only mode, using default properties.
  file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* register the filename in the mesh struct. */
  mesh->filename = filename;
	
  /* read mesh dimension. */
  {
	uint dim;
	hdf5_read_data_dim(file, "dat_coords", &dim);
	mesh->dim = dim;
  }
  

  /* read mesh coordinates. */
  {
    uint sizes[2];
    sizes[0] = mesh->sizes.nodes;
    sizes[1] = mesh->dim;

	double* data = (double*)malloc(sizes[0]*sizes[1]*sizeof(double));
	
	hdf5_read_data(file, "dat_coords", "double", data);
	
	// copy the node coordinates to the mesh struct.
	
	mesh->vertices = (Real**)malloc(sizes[0]*sizeof(Real*));
	if ( mesh->vertices ) {
		for (uint i = 0; i < sizes[0]; i++){
			mesh->vertices[i] = (Real*)malloc( mesh->dim*sizeof(Real) );
		}
	} else {
		printf("\nError: could not allocate memory for mesh vertices.\n");
		exit(1);
	}
	
    if ( mesh->dim == 2 ) {
      for( uint i = 0; i < sizes[0]; i++ ) {
        mesh->vertices[i][0] = data[i*2]  ;
        mesh->vertices[i][1] = data[i*2+1];
      }
    }
    else if ( mesh->dim == 3 ) {
      for( uint i = 0; i < sizes[0]; i++ ) {
        mesh->vertices[i][0] = data[i*3  ];
        mesh->vertices[i][1] = data[i*3+1];
        mesh->vertices[i][2] = data[i*3+2];
      }
    }
	
	free(data);
  }

  /* read mesh mappings */
  {
	int num_quad;
	hdf5_read_scalar(file, "num_quad", "int", &num_quad );
	mesh->sizes.elem_type_count[pQUAD_4] = num_quad;
	
	uint sizes[2];
    sizes[0] = num_quad;
    sizes[1] = 4;
    uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));
	
	hdf5_read_data(file, "map_ElemNode_quad", "int", data);
	
	// copy the connectivity table to the mesh struct.
	
	mesh->mElems = (MapElems*)malloc( mesh->sizes.elems*sizeof(MapElems) );
	if ( mesh->mElems ) {
		for (uint i = 0; i < mesh->sizes.elems; i++){
			mesh->mElems[i].nodes = (uint*)malloc( MAX_ELEM_NODES*sizeof(uint) );
			mesh->mElems[i].nbs   = NULL;
			mesh->mElems[i].faces = NULL;
		}
	} else {
		printf("\nError: could not allocate memory for mesh elements connectivity.\n");
		exit(1);
	}
	
	uint k = 0;
	for (uint i = 0; i < mesh->sizes.elems; i++)
	{
		for(uint j = 0; j < sizes[1]; j++ )
		{
			mesh->mElems[i].nodes[j] = data[k++];
		}
	}
	
	free(data);
  }

  
  /*
  Note that the following maps are ignored:
    map_FaceNode
    map_FaceElem
    map_BFaceElem
  Such information is not required for the CGNS output.
  */
  
  /* read boundary face to nodes mapping, mbface. */
  {
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      sizes[0] = mesh->sizes.bfaces;
      sizes[1] = 2;
      uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

	  hdf5_read_data(file, "map_BFaceNode", "int", data);
	  
	mesh->mBFaces = (MapBFaces*)malloc( sizes[0]*sizeof(MapBFaces) );
	if ( mesh->mBFaces ) {
		for (uint i = 0; i < mesh->sizes.bfaces; i++){
			mesh->mBFaces[i].nodes = (uint*)malloc( 2*sizeof(uint) );
	}
	} else {
		printf("\nError: could not allocate memory for mesh boundary faces.\n");
		exit(1);
	}
	  
      for( uint i = 0; i < sizes[0]; i++ ) {
        mesh->mBFaces[i].nodes[0] = data[i*2  ];
        mesh->mBFaces[i].nodes[1] = data[i*2+1];
      }
	  
      free(data);
    }
    else if ( mesh->dim == 3 ) {
      
    }
  }
  
  H5Fclose (file);
}
