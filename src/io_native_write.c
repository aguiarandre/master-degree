#include "native.h"

#include <stdlib.h> // malloc
#include <string.h>
#include <stddef.h>


/* Write meta information about the mesh. The idea is to record
 * the data defined in MeshStruct so we can read back the mesh
 * at a latter point in time.
*/
void mesh_write_meta_info(const char* filename, MeshStruct* mesh)
{

  // Create new HDF5 file using default properties.
  hid_t file = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

  /* *************************************** */
  /* write mesh size meta information */
  {
		
		typedef struct
		{
			uint nodes;
			uint elems;
			uint faces;
			uint bfaces;
			uint zones;
		} MetaSize;
		
	/* Calculate the size and the offsets of our struct members in memory */
	size_t dst_size = sizeof( MetaSize );
	size_t dst_offset[5] = { HOFFSET( MetaSize, nodes ),
							 HOFFSET( MetaSize, elems ),
                             HOFFSET( MetaSize, faces ),
                             HOFFSET( MetaSize, bfaces ),
							 HOFFSET( MetaSize, zones )};

	MetaSize p_data[1] = {
	{mesh->sizes.nodes, mesh->sizes.elems,
	 mesh->sizes.faces, mesh->sizes.bfaces,
	 mesh->sizes.zones}
	};

	/* Define field information */
	const char *field_names[5]  =
	{ "Nodes", "Elems", "Faces", "BFaces", "Zones" };
	hid_t      field_type[5];
	hsize_t    chunk_size = 10;
	int        *fill_data = NULL;
	int        compress  = 0;

	/* Initialize field_type */
	field_type[0] = H5T_NATIVE_INT;
	field_type[1] = H5T_NATIVE_INT;
	field_type[2] = H5T_NATIVE_INT;
	field_type[3] = H5T_NATIVE_INT;
	field_type[4] = H5T_NATIVE_INT;
	
	H5TBmake_table( "meta_sizes", file, "meta_sizes", 5, 1,
                     dst_size, field_names, dst_offset, field_type,
                     chunk_size, fill_data, compress, p_data );
  }
  
  
  /* *************************************** */
  /* write mesh zones meta information */
  {
	typedef struct
	{
	  uint          idx;
	  char          name[STRING_MAX_LENGHT];
	  int           dim;
	  int           is_boundary;
	  char          elem_type[STRING_MAX_LENGHT];
	  int           homogeneous;
	  int           num_types;
	  int           start;
	  int           num_elems;
	  int			bc_id;
	} MetaZone;
	
	/* Calculate the size and the offsets of our struct members in memory */
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

	int num_zones = mesh->sizes.zones;
	MetaZone* p_data = (MetaZone*)malloc( num_zones*sizeof(MetaZone) );
	
	for( int z = 0; z < num_zones; z++ )
	{
		p_data[z].idx  = mesh->zones[z].idx;
		strcpy(p_data[z].name, mesh->zones[z].name);
		p_data[z].dim  = mesh->zones[z].dim;
		p_data[z].is_boundary = (int)mesh->zones[z].is_boundary;
		
		switch( (int)mesh->zones[z].elem_type )
		{
			case 2:
				strcpy(p_data[z].elem_type,"BAR_2");
				break;
			case 3:
				strcpy(p_data[z].elem_type,"TRI_3");
				break;
			case 4:
				strcpy(p_data[z].elem_type,"TRI_6");
				break;
			case 5:
				strcpy(p_data[z].elem_type,"QUAD_4");
				break;
			case 6:
				strcpy(p_data[z].elem_type,"QUAD_9");
				break;
			case 7:
				strcpy(p_data[z].elem_type,"TETRA_4");
				break;
			case 8:
				strcpy(p_data[z].elem_type,"TETRA_10");
				break;
			case 9:
				strcpy(p_data[z].elem_type,"PYRA_5");
				break;
			case 10:
				strcpy(p_data[z].elem_type,"PENTA_6");
				break;
			case 11:
				strcpy(p_data[z].elem_type,"HEXA_8");
				break;
			case 12:
				strcpy(p_data[z].elem_type,"MIXED");
				break;
			default:
				break;
		}
		
		p_data[z].homogeneous = (int)mesh->zones[z].homogeneous;
		p_data[z].num_types   = mesh->zones[z].num_types;
		p_data[z].start       = mesh->zones[z].start;
		p_data[z].num_elems   = mesh->zones[z].num_elems;
		p_data[z].bc_id		  = mesh->zones[z].bc_id;
	}
	
	/* Define field information */
	const char *field_names[10]  =
	{ "idx", "name", "dim", "is_boundary", "elem_type",
	  "homogeneous", "num_types", "start", "num_elems", "bc_id" };
	hid_t      field_type[10];
	hid_t      string_type;
	hsize_t    chunk_size = 10;
	int        *fill_data = NULL;
	int        compress  = 0;

	/* Initialize field_type */
	string_type = H5Tcopy( H5T_C_S1 );
	H5Tset_size( string_type, STRING_MAX_LENGHT );
	field_type[0] = H5T_NATIVE_INT;
	field_type[1] = string_type;
	field_type[2] = H5T_NATIVE_INT;
	field_type[3] = H5T_NATIVE_INT;
	field_type[4] = string_type;
	field_type[5] = H5T_NATIVE_INT;
	field_type[6] = H5T_NATIVE_INT;
	field_type[7] = H5T_NATIVE_INT;
	field_type[8] = H5T_NATIVE_INT;
	field_type[9] = H5T_NATIVE_INT;

	H5TBmake_table( "meta_zones", file, "meta_zones", 10, num_zones,
                     dst_size, field_names, dst_offset, field_type,
                     chunk_size, fill_data, compress, p_data );

	free(p_data);
	H5Tclose( string_type );
  }
  
  
  /* *************************************** */
  /* write mesh elements meta information */
	{
		typedef struct
		{
			char elem_type_name[STRING_MAX_LENGHT];
			int size;
			int nodes;
		} MetaElem;
		
	/* Calculate the size and the offsets of our struct members in memory */
	size_t dst_size = sizeof( MetaElem );
	size_t dst_offset[3] = { HOFFSET( MetaElem, elem_type_name ),
							 HOFFSET( MetaElem, size ),
               HOFFSET( MetaElem, nodes )};
		
	int num_type_elems = MAX_ELEM_TYPES;
	MetaElem* p_data = (MetaElem*)malloc( num_type_elems*sizeof(MetaElem) );
	
	int j = 0;
	for( int i = 0; i < num_type_elems; i++ )
	{
		switch(i)
		{
			case 2:
				strcpy(p_data[j].elem_type_name,"BAR_2");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 2;
				j++;
				break;
			case 3:
				strcpy(p_data[j].elem_type_name,"TRI_3");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 3;
				j++;
				break;
			case 4:
				strcpy(p_data[j].elem_type_name,"TRI_6");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 6;
				j++;
				break;
			case 5:
				strcpy(p_data[j].elem_type_name,"QUAD_4");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 4;
				j++;
				break;
			case 6:
				strcpy(p_data[j].elem_type_name,"QUAD_9");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 9;
				j++;
				break;
			case 7:
				strcpy(p_data[j].elem_type_name,"TETRA_4");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 4;
				j++;
				break;
			case 8:
				strcpy(p_data[j].elem_type_name,"TETRA_10");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 10;
				j++;
				break;
			case 9:
				strcpy(p_data[j].elem_type_name,"PYRA_5");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 5;
				j++;
				break;
			case 10:
				strcpy(p_data[j].elem_type_name,"PENTA_6");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 6;
				j++;
				break;
			case 11:
				strcpy(p_data[j].elem_type_name,"HEXA_8");
				p_data[j].size = mesh->sizes.elem_type_count[i];
				p_data[j].nodes = 8;
				j++;
				break;
			default:
				break;
		}
	}
	int num_valid_elems = j--;

	/* Define field information */
	const char *field_names[3]  =
	{ "ElemType", "Size", "Nodes" };
	hid_t      field_type[3];
	hid_t 		 string_type;
	hsize_t    chunk_size = 10;
	int        *fill_data = NULL;
	int        compress  = 0;

	/* Initialize field_type */
	string_type = H5Tcopy( H5T_C_S1 );
	H5Tset_size( string_type, STRING_MAX_LENGHT );
	field_type[0] = string_type;
	field_type[1] = H5T_NATIVE_INT;
	field_type[2] = H5T_NATIVE_INT;
	
	H5TBmake_table( "meta_elems", file, "meta_elems", 3, num_valid_elems,
                     dst_size, field_names, dst_offset, field_type,
                     chunk_size, fill_data, compress, p_data );

		free(p_data);
		H5Tclose( string_type );
	}
  
  
  /* *************************************** */
  /* write mesh meta information */
  {
		
		typedef struct
		{
			char filename[STRING_MAX_LENGHT];
			int idx;
			int dim;
		} MetaMesh;
		
	/* Calculate the size and the offsets of our struct members in memory */
	size_t dst_size = sizeof( MetaMesh );
	size_t dst_offset[3] = { HOFFSET( MetaMesh, filename ),
								HOFFSET( MetaMesh, idx ),
                HOFFSET( MetaMesh, dim )};

	MetaMesh p_data[1];
	strcpy(p_data[0].filename,mesh->filename);
	p_data[0].idx = mesh->idx;
	p_data[0].dim = mesh->dim;

	/* Define field information */
	const char *field_names[3]  =
	{ "FileName", "Index", "Dimension" };
	hid_t      field_type[3];
	hid_t  		 string_type;
	hsize_t    chunk_size = 10;
	int        *fill_data = NULL;
	int        compress  = 0;

	/* Initialize field_type */
	string_type = H5Tcopy( H5T_C_S1 );
	H5Tset_size( string_type, STRING_MAX_LENGHT );
	field_type[0] = string_type;
	field_type[1] = H5T_NATIVE_INT;
	field_type[2] = H5T_NATIVE_INT;

	
	H5TBmake_table( "meta_mesh", file, "meta_mesh", 3, 1,
                     dst_size, field_names, dst_offset, field_type,
                     chunk_size, fill_data, compress, p_data );

	H5Tclose( string_type );
  }
  
  
  H5Fclose( file );
}


/* Here we output the native mesh to disk,
 * including the MeshStruct data. The hdf5
 * format is used as the native file format.
 */
void mesh_write_data_native(const char* filename, MeshStruct* mesh)
{
  printf("\nWriting native mesh file.\n");
  fflush(stdout);

  hid_t file;

  // Open HDF5 file using default properties.
  file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );

  
  /* *************************************** */
  /* write mesh coordinates */

  {
    uint sizes[2];
    sizes[0] = mesh->sizes.nodes;
    sizes[1] = mesh->dim;

    double* data = (double*)malloc(sizes[0]*sizes[1]*sizeof(double));
    if ( mesh->dim == 2 ) {
      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i*2  ] = mesh->vertices[i][0];
        data[i*2+1] = mesh->vertices[i][1];
      }
    }
    else if ( mesh->dim == 3 ) {
      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i*3  ] = mesh->vertices[i][0];
        data[i*3+1] = mesh->vertices[i][1];
        data[i*3+2] = mesh->vertices[i][2];
      }
    }

    hdf5_write_data(file, "dat_coords", "double", sizes, 1, data);

    free(data);
  }

  /* *************************************** */
  /* write mesh mappings */
  
  /* write number of a particular mesh element type. */
  {
    uint size;

    if ( mesh->dim == 2 ) {

      // test for triangles in the mesh
      size = mesh->sizes.elem_type_count[pTRI_3];
      if ( size > 0 ) hdf5_write_scalar(file, "num_tri" , "int", &size );
      // test for quads in the mesh
      size = mesh->sizes.elem_type_count[pQUAD_4];
      if ( size > 0 ) hdf5_write_scalar(file, "num_quad" , "int", &size );

    }
    else if ( mesh->dim == 3 ) {

      // test for tetras in the mesh
      size = mesh->sizes.elem_type_count[pTETRA_4];
      if ( size > 0 ) hdf5_write_scalar(file, "num_tetra" , "int", &size );
      // test for pyra in the mesh
      size = mesh->sizes.elem_type_count[pPYRA_5];
      if ( size > 0 ) hdf5_write_scalar(file, "num_pyra" , "int", &size );
      // test for prism in the mesh
      size = mesh->sizes.elem_type_count[pPENTA_6];
      if ( size > 0 ) hdf5_write_scalar(file, "num_prism" , "int", &size );
      // test for hexa in the mesh
      size = mesh->sizes.elem_type_count[pHEXA_8];
      if ( size > 0 ) hdf5_write_scalar(file, "num_hexa" , "int", &size );

    }
  }

  /* mesh element to nodes mapping, melem,
   * stored by element type. */
  {
    uint elem_count;
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      // write mesh triangles
      elem_count = mesh->sizes.elem_type_count[pTRI_3];
      if ( elem_count > 0 ) {

        /* Get a list of element indexes in the mesh. */
        sizes[0] = elem_count;
        sizes[1] = 3;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pTRI_3, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_tri", "int", sizes, 0, data);

        free(indices);
        free(data);
      }
      // write mesh quads
      elem_count = mesh->sizes.elem_type_count[pQUAD_4];
      if ( elem_count > 0 ) {
        sizes[0] = elem_count;
        sizes[1] = 4;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pQUAD_4, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_quad", "int", sizes, 0, data);

        free(indices);
        free(data);
      }

    }
    else if ( mesh->dim == 3 ) {

      // write mesh tetra
      elem_count = mesh->sizes.elem_type_count[pTETRA_4];
      if ( elem_count > 0 ) {
        sizes[0] = elem_count;
        sizes[1] = 4;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pTETRA_4, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_tetra", "int", sizes, 0, data);

        free(indices);
        free(data);
      }
      // write mesh pyras
      elem_count = mesh->sizes.elem_type_count[pPYRA_5];
      if ( elem_count > 0 ) {
        sizes[0] = elem_count;
        sizes[1] = 5;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pPYRA_5, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_pyra", "int", sizes, 0, data);

        free(indices);
        free(data);
      }
      // write mesh prisms
      elem_count = mesh->sizes.elem_type_count[pPENTA_6];
      if ( elem_count > 0 ) {
        sizes[0] = elem_count;
        sizes[1] = 6;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pPENTA_6, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_prism", "int", sizes, 0, data);

        free(indices);
        free(data);
      }
      // write mesh hexas
      elem_count = mesh->sizes.elem_type_count[pHEXA_8];
      if ( elem_count > 0 ) {
        sizes[0] = elem_count;
        sizes[1] = 8;
        uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

        uint* indices = (uint*)malloc(sizes[0]*sizeof(uint));
        mesh_elems_by_type(pHEXA_8, mesh, indices);

        for( uint i = 0; i < sizes[0]; i++ ) {
          for ( short j = 0; j < sizes[1]; j++ ) {
            uint idx = i*sizes[1]+j;
            data[idx] = mesh->mElems[ indices[i] ].nodes[j];
          }
        }
        hdf5_write_data(file, "map_ElemNode_hexa", "int", sizes, 0, data);

        free(indices);
        free(data);
      }

    }
  }


  /* mesh face to nodes mapping, mface. */
  {
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      sizes[0] = mesh->sizes.faces;
      sizes[1] = 2;
      uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i*2  ] = mesh->mFaces[i].nodes[0];
        data[i*2+1] = mesh->mFaces[i].nodes[1];
      }
      hdf5_write_data(file, "map_FaceNode", "int", sizes, 0, data);

      free(data);

    }
    else if ( mesh->dim == 3 ) {

      /* internal faces in 3D are either a tri or a quad. */
      
      
      

    }

  }


  /* mesh face to cells mapping, mfelem. */
  {
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      sizes[0] = mesh->sizes.faces;
      sizes[1] = 2;
      uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i*2  ] = mesh->mFaces[i].elem_left;
        data[i*2+1] = mesh->mFaces[i].elem_right;
      }
      hdf5_write_data(file, "map_FaceElem", "int", sizes, 0, data);

      free(data);

    }
    else if ( mesh->dim == 3 ) {
      
      
      
    }

  }


  /* mesh boundary face to nodes mapping, mbface. */
  {
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      sizes[0] = mesh->sizes.bfaces;
      sizes[1] = 2;
      uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i*2  ] = mesh->mBFaces[i].nodes[0];
        data[i*2+1] = mesh->mBFaces[i].nodes[1];
      }
      hdf5_write_data(file, "map_BFaceNode", "int", sizes, 0, data);

      free(data);

    }
    else if ( mesh->dim == 3 ) {
      
    }
  }


  /* mesh boundary face to elems mapping, mbfelem. */
  {
    uint sizes[2];

    if ( mesh->dim == 2 ) {

      sizes[0] = mesh->sizes.bfaces;
      sizes[1] = 1;
      uint* data = (uint*)malloc(sizes[0]*sizes[1]*sizeof(uint));

      for( uint i = 0; i < sizes[0]; i++ ) {
        data[i] = mesh->mBFaces[i].elem_left;
      }
      hdf5_write_data(file, "map_BFaceElem", "int", sizes, 0, data);

      free(data);

    }
    else if ( mesh->dim == 3 ) {

    }
  }


  /* *************************************** */
  /* finish mesh file writting. */

  H5Fclose( file );

  printf("\nDone.\n");
  fflush(stdout);
}


/* Find the index of all mesh elements of the
 * requested type, by searching the mesh zones,
 * and return them in data.
 */
void mesh_elems_by_type(ElemType type, MeshStruct* mesh, uint* data)
{
  if ( mesh->sizes.elem_type_count[type] == 0 ) return;

  uint idx = 0;

  for ( uint z = 0; z < mesh->sizes.zones; z++ ) {

    if ( mesh->zones[z].homogeneous ) {
      /* we are in an homogeneous mesh zone. */
      if ( mesh->zones[z].elem_type == type ) {
        uint beg = mesh->zones[z].start;
        uint end = beg + mesh->zones[z].num_elems;
        for ( uint i = beg; i < end; i++ )
            data[idx++] = i;
      }
    }
    else {
      /* we are in a mixed mesh zone. */
      for( short n = 0; n < mesh->zones[z].num_types; n++ ) {
        if ( mesh->zones[z].type_list[n].elem_type == type ) {
          uint beg = mesh->zones[z].type_list[n].start;
          uint end = beg + mesh->zones[z].type_list[n].num_elems;
          for ( uint i = beg; i < end; i++ )
              data[idx++] = i;
        }
      }
    }

  }
}
