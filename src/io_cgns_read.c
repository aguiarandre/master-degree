#include <stdio.h>
#include <string.h> /* memcpy() */
#include <stdlib.h>
#include <math.h>

#include "io_cgns.h"

void io_cgns_read(const char* filename, MeshStruct* mesh)
{
  int cgfile, nbases, nzones, nsections;
  cgsize_t size[3];
  char zone_name[STRING_MAX_LENGHT];

  if ( cg_open( filename, CG_MODE_READ, &cgfile ) ) cg_error_exit();

  cg_nbases( cgfile, &nbases );

  if ( nbases > 1 ) {
    printf("CGNS file nbases greater than one is not supported.\n");
    exit(1);
  }

  for ( int baseid = 1; baseid <= nbases; baseid++ )
  {

    cg_nzones( cgfile, baseid, &nzones );

    if ( nzones > 1 ) {
      printf("CGNS file nzones greater than one is not supported.\n");
      exit(1);
    }

    for ( int zoneid = 1; zoneid <= nzones; zoneid++ )
    {
      cg_zone_read( cgfile, baseid, zoneid, zone_name, size );

      cg_nsections( cgfile, baseid, zoneid, &nsections );

      printf("\nCGNS zone name     = %s"     , zone_name);
      printf("\nNumber of nodes    = %d", (uint) size[0]);
      printf("\nNumber of elements = %d", (uint) size[1]);
      printf("\nNumber of parts    = %d", nsections);
      fflush(stdout);

	  mesh->filename = filename;
      mesh->sizes.nodes = size[0];
      mesh->sizes.elems = size[1];
      mesh->sizes.zones = nsections;

      io_cgns_read_nodes(cgfile, baseid, zoneid, size, mesh);

      io_cgns_read_elems_connectivity(cgfile, baseid, zoneid,
                                      nsections, size, mesh);

    }

  }

  cg_close (cgfile);

  printf("\nDone.\n");
  fflush(stdout);
}


void io_cgns_read_nodes(int cgfile, int baseid, int zoneid,
                        cgsize_t* size, MeshStruct* mesh)
{
  printf("\nReading mesh nodes...");
  fflush(stdout);
  
  



  mesh->vertices = (Real**)malloc(size[0]*sizeof(Real*));
  if ( mesh->vertices ) {
    for (uint i = 0; i < size[0]; i++){
      mesh->vertices[i] = (Real*)malloc( 3*sizeof(Real) ); // x y z
    }
    

    
    
  } else {
    printf("\nError: could not allocate memory for mesh vertices.\n");
    exit(1);
  }

  /* temporary array for transfer */
  Real* coord = (Real*)malloc( size[0]*sizeof(Real) );
  if ( coord == NULL ) {
    printf("\nError: could not allocate memory for mesh vertices.\n");
    exit(1);
  }

  cgsize_t start = 1;

  /* Read nodes coordinates */
  cg_coord_read( cgfile, baseid, zoneid, "CoordinateX",
                 RealDouble, &start, &size[0], &coord[0]);

  for (uint i = 0; i < size[0]; i++)
    mesh->vertices[i][0] = coord[i];

  cg_coord_read( cgfile, baseid, zoneid, "CoordinateY",
                 RealDouble, &start, &size[0], &coord[0] );

  for (uint i = 0; i < size[0]; i++)
    mesh->vertices[i][1] = coord[i];

  /* Apparently, mesh generators always write the CoordinateZ
   * field, regardless of the mesh dimension. Check here if
   * we are dealing with a 3D or 2D planar mesh.
   */
  cg_coord_read( cgfile, baseid, zoneid, "CoordinateZ",
                 RealDouble, &start, &size[0], &coord[0] );

  mesh->dim = 2;
  Real eps = 1.0e-8;
  for (uint i = 1; i < size[0]; i++) {
    if ( fabs(coord[i] - coord[i-1]) > eps ) mesh->dim = 3;
  }

  if ( mesh->dim == 3 ) {
    for (uint i = 0; i < size[0]; i++)
      mesh->vertices[i][2] = coord[i];
  }

  free(coord);
}

void io_cgns_read_elems_connectivity(int cgfile, int baseid, int zoneid,
                                     int nsections, cgsize_t* size, MeshStruct* mesh)
{
  printf("\nReading elements connectivity...");
  fflush(stdout);

  /* go through the sections once to obtain sizes and
   * element types to allocate the appropriate arrays.
   * We also define the order in which the zones should
   * be created. Note that CGNS begins arrays indexing
   * with 1, so we offset that to start at 0.
   */

  int sec_order[nsections],
      sec_start[nsections];

  bool is_boundary[nsections];

  int num_boundary_zones = 0,
      num_internal_zones = 0;

  uint num_boundary_faces = 0;

  for ( int secid = 0; secid < nsections; secid++ ) {

    char section_name[STRING_MAX_LENGHT];
    ElementType_t etype;
    cgsize_t start;
    cgsize_t end;
    int nbndry, parent_flag;
    int cgsec = secid+1;

    cg_section_read( cgfile, baseid, zoneid, cgsec, section_name,
                     &etype, &start, &end, &nbndry, &parent_flag );

    sec_order[secid] = secid;
    sec_start[secid] = start;
    uint sec_size = end - start +1;

    /* Determine if this is a boundary zone. Note that
     * all zones are initialized with is_boundary = false.
     */
    is_boundary[secid] = false;

    if ( mesh->dim == 2 && etype == BAR_2 ) {
      is_boundary[secid] = true;
      num_boundary_zones += 1;
      num_boundary_faces += sec_size;
    }

    if ( mesh->dim == 3 ){

      if ( etype == TRI_3 || etype == QUAD_4 ) {

        /* homogeneous element type. */
        is_boundary[secid] = true;
        num_boundary_zones += 1;
        num_boundary_faces += sec_size;

      }
      /* Note that MIXED is a CGNS enum. */
      else if ( etype == MIXED ) {

        /* temporary array to hold element connectivity. */
        cgsize_t ElementDataSize;
        cgsize_t ParentData ;

        cg_ElementDataSize( cgfile, baseid, zoneid, cgsec, &ElementDataSize );
        cgsize_t* Elements = (cgsize_t*)malloc(ElementDataSize*sizeof(cgsize_t));
        cg_elements_read( cgfile, baseid, zoneid, cgsec, &Elements[0], &ParentData );

        /* check for mixed element types. */
        int node_count;
        uint connIndex = 0;

        for ( uint i = 0; i < sec_size; i++ ) {
          /* get current element type. */
          cg_npe( Elements[connIndex], &node_count );
          /* check if we have shell elements.]
           * Note that TETRA_4 is a CGNS enum defined in cgnslib.h */
          if ( Elements[connIndex] < TETRA_4 ) {
            is_boundary[secid] = true;
          } else {
            is_boundary[secid] = false;
            break;
          }
          connIndex += node_count +1;
        }
        if ( is_boundary[secid] ) {
          num_boundary_zones += 1;
          num_boundary_faces += sec_size;
        }
        /* deallocate temporary array. */
        free( Elements );
      }
    } /* mesh->dim == 3 if */

  } /* secid for loop */

  num_internal_zones = nsections - num_boundary_zones;

  /* sort sections by elements index. Brute force sort should
   * be OK here. We wont have a huge number of sections.
   */
  /* sort sections by index. */
  for ( int i = 0; i < nsections; i++ ){
    for ( int j= 0; j < nsections; j++ ){
      if ( i != j ) {
        int tmp, aux;
        if ( sec_start[i] < sec_start[j] ) {
          tmp = sec_start[j];
          aux = sec_order[j];
          sec_start[j] = sec_start[i];
          sec_start[i] = tmp;
          sec_order[j] = sec_order[i];
          sec_order[i] = aux;
        }
      }
    }
  }

  /* index the zones by type. */
  int zone_internal[num_internal_zones],
      zone_boundary[num_boundary_zones];

  int j = 0, k = 0;
  for ( int secid = 0; secid < nsections; secid++ ){
    int index = sec_order[secid];
    if ( is_boundary[index] ) {
      zone_boundary[j++] = index;
    } else {
      zone_internal[k++] = index;
    }
  }

  /* --------------------- allocate structs --------------------- */

  /* allocate the mesh zones structs. Note that CGNS zones and
   * native mesh zones are not the same thing.
   */
  mesh->zones = (MeshZone*)malloc(nsections*sizeof(MeshZone));
  if ( mesh->zones == NULL ) {
    printf("\nError: could not allocate memory for mesh zones.\n");
    exit(1);
  }
  mesh->sizes.bfaces = num_boundary_faces;

  /* allocate mesh mappings */

  mesh->mBFaces = (MapBFaces*)malloc( mesh->sizes.bfaces*sizeof(MapBFaces) );
  if ( mesh->mBFaces ) {
    for (uint i = 0; i < mesh->sizes.bfaces; i++){
      mesh->mBFaces[i].nodes = (uint*)malloc( MAX_FACE_NODES*sizeof(uint) );
      for ( short j = 0; j < MAX_FACE_NODES; j++ ) {
        mesh->mBFaces[i].nodes[j] = -1;
      }
    }
  } else {
    printf("\nError: could not allocate memory for mesh boundary faces.\n");
    exit(1);
  }

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

  /* --------------------- read internal sections --------------------- */

  /* second pass, actual reading of the connectivity of sections.
   * Note that we follow an specific order for the zone creations
   * in the native structure.
   */

  /* counter for mixed mesh zones. */
  uint zone_type_start = 0;

  for ( int secid = 0; secid < num_internal_zones; secid++ ) {

    char section_name[STRING_MAX_LENGHT];
    ElementType_t etype;
    cgsize_t start, end;
    uint size;
    int nbndry, parent_flag;
    int cgsec = zone_internal[secid]+1;

    cg_section_read( cgfile, baseid, zoneid, cgsec, section_name,
                     &etype, &start, &end, &nbndry, &parent_flag );

    mesh_zone_init( &mesh->zones[secid] );
    mesh->zones[secid].idx = secid;
    mesh->zones[secid].name = (char*)malloc(STRING_MAX_LENGHT*sizeof(char));
    sprintf(mesh->zones[secid].name, "%s",section_name);
    size = (uint)(end - start +1);
    mesh->zones[secid].num_elems = size;
    mesh->zones[secid].start = start-1;
    mesh->zones[secid].is_boundary = false;

    int node_count;

    /* figure out which is the type of element of this section */
    switch ( etype ) {
      case TRI_3:
        mesh->zones[secid].elem_type = pTRI_3;
        mesh->sizes.elem_type_count[pTRI_3] += size;
        node_count = 3;
        break;
      case QUAD_4:
        mesh->zones[secid].elem_type = pQUAD_4;
        mesh->sizes.elem_type_count[pQUAD_4] += size;
        node_count = 4;
        break;
      case TETRA_4:
        mesh->zones[secid].elem_type = pTETRA_4;
        mesh->sizes.elem_type_count[pTETRA_4] += size;
        node_count = 4;
        break;
      case PYRA_5:
        mesh->zones[secid].elem_type = pPYRA_5;
        mesh->sizes.elem_type_count[pPYRA_5] += size;
        node_count = 5;
        break;
      case PENTA_6:
        mesh->zones[secid].elem_type = pPENTA_6;
        mesh->sizes.elem_type_count[pPENTA_6] += size;
        node_count = 6;
        break;
      case HEXA_8:
        mesh->zones[secid].elem_type = pHEXA_8;
        mesh->sizes.elem_type_count[pHEXA_8] += size;
        node_count = 8;
        break;
      case MIXED:
        mesh->zones[secid].elem_type = pMIXED;
        node_count = 0;
        break;
      default:
        mesh->zones[secid].elem_type = pNONE;
        node_count = 0;
        break;
    }

    /* Determine zone dimension. */
    if ( mesh->dim == 2 )
      mesh->zones[secid].dim = 2;
    if ( mesh->dim == 3 )
      mesh->zones[secid].dim = 3;

    /* Some mesh generators don't write the correct CGNS type for
     * some elements, possibly due to an old cgns lib version.
     * If the zone type was not found above, search for it again here.
     * We also check if all the elements in this zone have the same type.
     */
    bool homogeneous;

    if ( etype == TETRA_4 || etype == PYRA_5 || etype == PENTA_6 || etype == HEXA_8 ||
         etype == TRI_3   || etype == QUAD_4 ) {
      mesh->zones[secid].homogeneous = true;
      homogeneous = true;
    } else {
      mesh->zones[secid].homogeneous = false;
      homogeneous = false;
    }

    /* Read the zone elements connectivity. For non-mixed families,
      * the nodes are listed all together, without a type prefix.
      */
    cgsize_t ElementDataSize;
    cgsize_t ParentData ;

    cg_ElementDataSize( cgfile, baseid, zoneid, cgsec, &ElementDataSize );

    /* temporary array to hold element connectivity. */
    cgsize_t* Elements = (cgsize_t*)malloc(ElementDataSize*sizeof(cgsize_t));

    cg_elements_read( cgfile, baseid, zoneid, cgsec, &Elements[0], &ParentData );

    /* fetch current zone connectivity data to the
      * native mesh data structure. */

    if ( homogeneous ) {

      uint connIndex = 0;

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        for (uint j = 0; j < node_count; j++) {
          /* we need to offset the index, since cgns starts at 1. */
          mesh->mElems[i].nodes[j] = Elements[connIndex+j] -1;
        }
        connIndex += node_count;
      }

      zone_type_start += mesh->zones[secid].num_elems;

    } else {

      /* Zones with mixed element types will be populated differently.
       * We shall order the elements in the map by its type, from lower
       * to greater element type. This will allow us to loop over elements
       * of the same type once at a time, avoiding a bunch of if checks
       * by the solver. This also will save us memory and disk space since
       * we will have only to save the start and end indexes of such common
       * regions, instead of an array[num_elems] with element types.
       * 
       * CGNS typically have elements already ordered in its mixed sections.
       * However, I don't think thats a requirement of the standard and so we
       * will check it here. The mixed sections connectivity contents are:
       *
       * ElemType
       * node1
       * node2
       * node3
       * ElemType
       * node1
       * node2
       * node3
       * node4
       * ...
       * So we must first read "ElemType" to know how many vertices to read.
       * This is acomplished by the cg_npe() function.
       */

      /* get current element node count/type. Note that we
      * are in an internal mesh zone.
      */
      uint elem_type_count[MAX_ELEM_TYPES] = {0};
      uint connIndex = 0;

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        cg_npe( Elements[connIndex], &node_count );
        if ( mesh->dim == 2 ) {
          if ( node_count == 3 ) elem_type_count[ pTRI_3  ] += 1;
          if ( node_count == 4 ) elem_type_count[ pQUAD_4 ] += 1;
        }
        else if ( mesh->dim == 3 ) {
          if ( node_count == 4 ) elem_type_count[ pTETRA_4 ] += 1;
          if ( node_count == 5 ) elem_type_count[ pPYRA_5  ] += 1;
          if ( node_count == 6 ) elem_type_count[ pPENTA_6 ] += 1;
          if ( node_count == 8 ) elem_type_count[ pHEXA_8  ] += 1;
        }
        connIndex += node_count +1;
      }

      /* create temporary array for sorted indexes. */

      short zone_num_types = 0;

      uint *sorted_tri, *sorted_quad;
      sorted_tri = sorted_quad = NULL;

      uint *sorted_tetra, *sorted_pyra, *sorted_penta, *sorted_hexa;
      sorted_tetra = sorted_pyra = sorted_penta = sorted_hexa = NULL;

      if ( mesh->dim == 2 ) {

        if ( elem_type_count[pTRI_3] ) { 
          sorted_tri = (uint*)malloc( elem_type_count[pTRI_3]*sizeof(uint) );
          zone_num_types += 1;
        }
        if ( elem_type_count[pQUAD_4] ) {
          sorted_quad = (uint*)malloc( elem_type_count[pQUAD_4]*sizeof(uint) );
          zone_num_types += 1;
        }
      }
      else if ( mesh->dim == 3 ) {

        if ( elem_type_count[pTETRA_4] ) {
          sorted_tetra = (uint*)malloc( elem_type_count[pTETRA_4]*sizeof(uint) );
          zone_num_types += 1;
        }
        if ( elem_type_count[pPYRA_5] ) {
          sorted_pyra = (uint*)malloc( elem_type_count[pPYRA_5]*sizeof(uint) );
          zone_num_types += 1;
        }
        if ( elem_type_count[pPENTA_6] ) {
          sorted_penta = (uint*)malloc( elem_type_count[pPENTA_6]*sizeof(uint) );
          zone_num_types += 1;
        }
        if ( elem_type_count[pHEXA_8] ) {
          sorted_hexa = (uint*)malloc( elem_type_count[pHEXA_8]*sizeof(uint) );
          zone_num_types += 1;
        }
      }

      /* setup the zones mixed struct data. */
      mesh->zones[secid].num_types = zone_num_types;
      mesh->zones[secid].type_list = (ElemTypeList*)malloc( zone_num_types*sizeof(ElemTypeList) );

      /* initialize element type indexes. */

      uint idx_tri, idx_quad;
      idx_tri = idx_quad = 0;

      uint idx_tetra, idx_pyra, idx_penta, idx_hexa;
      idx_tetra = idx_pyra = idx_penta = idx_hexa = 0;

      /* sort out the zone element indexes by type. */
      connIndex = 0;
      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        cg_npe( Elements[connIndex], &node_count );
        /* Note that we use the idx and then increment it. */
        if ( mesh->dim == 2 ) {
          if ( node_count == 3 ) sorted_tri [ idx_tri++  ] = connIndex;
          if ( node_count == 4 ) sorted_quad[ idx_quad++ ] = connIndex;
        }
        else if ( mesh->dim == 3 ) {
          if ( node_count == 4 ) sorted_tetra[ idx_tetra++ ] = connIndex;
          if ( node_count == 5 ) sorted_pyra [ idx_pyra++  ] = connIndex;
          if ( node_count == 6 ) sorted_penta[ idx_penta++ ] = connIndex;
          if ( node_count == 8 ) sorted_hexa [ idx_hexa++  ] = connIndex;
        }
        connIndex += node_count +1;
      }

      /* Now we combine the various sorted indexes array
       * into a single one. We also setup the type_list struct
       * for this zone.
       */
      uint* sorted_index = (uint*)malloc( mesh->zones[secid].num_elems*sizeof(uint) );
      if ( sorted_index == NULL ) {
        printf("\nError: could not allocate memory for sorted_index array.\n");
        exit(1);
      }
      uint idx_sorted = 0;
      short typeid = 0;
      if ( mesh->dim == 2 ) {
        if ( elem_type_count[pTRI_3] ) {
          mesh->sizes.elem_type_count[pTRI_3] += elem_type_count[pTRI_3];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pTRI_3;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pTRI_3];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_tri[0], elem_type_count[pTRI_3]*sizeof(uint) );
          idx_sorted += elem_type_count[pTRI_3];
          zone_type_start += elem_type_count[pTRI_3];
          free(sorted_tri);
        }
        if ( elem_type_count[pQUAD_4] ) {
          mesh->sizes.elem_type_count[pQUAD_4] += elem_type_count[pQUAD_4];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pQUAD_4;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pQUAD_4];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_quad[0], elem_type_count[pQUAD_4]*sizeof(uint) );
          idx_sorted += elem_type_count[pQUAD_4];
          zone_type_start += elem_type_count[pQUAD_4];
          free(sorted_quad);
        }
      }
      else if ( mesh->dim == 3 ) {
        if ( elem_type_count[pTETRA_4] ) {
          mesh->sizes.elem_type_count[pTETRA_4] += elem_type_count[pTETRA_4];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pTETRA_4;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pTETRA_4];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_tetra[0], elem_type_count[pTETRA_4]*sizeof(uint) );
          idx_sorted += elem_type_count[pTETRA_4];
          zone_type_start += elem_type_count[pTETRA_4];
          free(sorted_tetra);
        }
        if ( elem_type_count[pPYRA_5] ) {
          mesh->sizes.elem_type_count[pPYRA_5] += elem_type_count[pPYRA_5];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pPYRA_5;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pPYRA_5];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_pyra[0], elem_type_count[pPYRA_5]*sizeof(uint) );
          idx_sorted += elem_type_count[pPYRA_5];
          zone_type_start += elem_type_count[pPYRA_5];
          free(sorted_pyra);
        }
        if ( elem_type_count[pPENTA_6] ) {
          mesh->sizes.elem_type_count[pPENTA_6] += elem_type_count[pPENTA_6];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pPENTA_6;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pPENTA_6];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_penta[0], elem_type_count[pPENTA_6]*sizeof(uint) );
          idx_sorted += elem_type_count[pPENTA_6];
          zone_type_start += elem_type_count[pPENTA_6];
          free(sorted_penta);
        }
        if ( elem_type_count[pHEXA_8] ) {
          mesh->sizes.elem_type_count[pHEXA_8] += elem_type_count[pHEXA_8];
          /* setup type_list struct */
          mesh->zones[secid].type_list[typeid].elem_type = pHEXA_8;
          mesh->zones[secid].type_list[typeid].start     = zone_type_start;
          mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pHEXA_8];
          typeid++;
          /* copy memory block to sorted index */
          memcpy( &sorted_index[idx_sorted], &sorted_hexa[0], elem_type_count[pHEXA_8]*sizeof(uint) );
          idx_sorted += elem_type_count[pHEXA_8];
          zone_type_start += elem_type_count[pHEXA_8];
          free(sorted_hexa);
        }
      }

      /* Finally, lets read the elements connectivity sorted by type. */

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {

        /* get current element node count/type. */
        cg_npe( Elements[ sorted_index[i] ], &node_count );

        /* we must offset the node_count value to account for the
         * element type flag. */
        for ( int j = 1; j <= node_count; j++ ) {
          /* we need to offset the index, since cgns starts at 1. */
          mesh->mElems[i].nodes[j-1] = Elements[ sorted_index[i]+j ] -1;
        }
      }

      /* free temporary arrays. */
      free(sorted_index);

    } /* homogeneous */

    free( Elements );

  } /* num_internal_zones for loop */


  /* --------------------- read bounday sections --------------------- */
  int zid = 0;
  uint bstart = 0;
  zone_type_start = 0; /* reset variable. */

  /* we start the section loop right after the internal zones index. */
  for ( int secid = num_internal_zones; secid < mesh->sizes.zones; secid++ ) {

    char section_name[STRING_MAX_LENGHT];
    ElementType_t etype;
    cgsize_t start, end;
    uint size;
    int nbndry, parent_flag;
    int cgsec = zone_boundary[zid]+1;

    cg_section_read( cgfile, baseid, zoneid, cgsec, section_name,
                     &etype, &start, &end, &nbndry, &parent_flag );

    mesh_zone_init( &mesh->zones[secid] );
    mesh->zones[secid].idx = secid;
    mesh->zones[secid].name = (char*)malloc(STRING_MAX_LENGHT*sizeof(char));
    sprintf(mesh->zones[secid].name,"%s", section_name);
    size = (uint)(end - start +1);
    mesh->zones[secid].num_elems = size;
    mesh->zones[secid].start = bstart; /* boundary faces have their own array. */
    bstart += size; /* setup next zone start index. */
    mesh->zones[secid].is_boundary = true;

    int node_count;

    /* figure out which is the type of element of this section */
    switch ( etype ) {
      case BAR_2:
        mesh->zones[secid].elem_type = pBAR_2;
        mesh->sizes.elem_type_count[pBAR_2] += size;
        node_count = 2;
        break;
      case TRI_3:
        mesh->zones[secid].elem_type = pTRI_3;
        mesh->sizes.elem_type_count[pTRI_3] += size;
        node_count = 3;
        break;
      case QUAD_4:
        mesh->zones[secid].elem_type = pQUAD_4;
        mesh->sizes.elem_type_count[pQUAD_4] += size;
        node_count = 4;
        break;
      case MIXED:
        mesh->zones[secid].elem_type = pMIXED;
        node_count = 0;
        break;
      default:
        mesh->zones[secid].elem_type = pNONE;
        node_count = 0;
        break;
    }

    /* Determine zone dimension. */
    if ( mesh->dim == 2 )
      mesh->zones[secid].dim = 1;
    if ( mesh->dim == 3 )
      mesh->zones[secid].dim = 2;

    /* Some mesh generators don't write the correct CGNS type for
     * some elements, possibly due to an old cgns lib version.
     * If the zone type was not found above, search for it again here.
     * We also check if all the elements in this zone have the same type.
     */
    bool homogeneous;

    if ( etype == TRI_3 || etype == QUAD_4 || ( ( etype == BAR_2 ) && mesh->dim == 2 ) ) {
      mesh->zones[secid].homogeneous = true;
      homogeneous = true;
    } else {
      mesh->zones[secid].homogeneous = false;
      homogeneous = false;
    }

    /* Read the zone elements connectivity. For non-mixed families,
      * the nodes are listed all together, without a type prefix.
      */
    cgsize_t ElementDataSize;
    cgsize_t ParentData ;

    cg_ElementDataSize( cgfile, baseid, zoneid, cgsec, &ElementDataSize );

    /* temporary array to hold element connectivity. */
    cgsize_t* Elements = (cgsize_t*)malloc(ElementDataSize*sizeof(cgsize_t));

    cg_elements_read( cgfile, baseid, zoneid, cgsec, &Elements[0], &ParentData );

    /* fetch current zone connectivity data to the
      * native mesh data structure. */

    if ( homogeneous ) {

      uint connIndex = 0;
      uint start = mesh->zones[secid].start;

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        for (uint j = 0; j < node_count; j++) {
          /* we need to offset the index, since cgns starts at 1. */
          mesh->mBFaces[start+i].nodes[j] = Elements[connIndex+j] -1;
        }
        connIndex += node_count;
      }

      zone_type_start += mesh->zones[secid].num_elems;

    } else {

      /* We are in a boundary zone with mixed elements. The same
       * remarks for mixed internal mesh zones is valid here.
       */

      if ( mesh->dim == 2 ) {
        /* Note that mixed boundary zones can only happen with mesh->dim == 3. */
        printf("\nError: found a mixed element boundary section for DIM=2.\n");
        exit(1);
      }

      /* get current element node count/type. Note that we
      * are in a boundary mesh zone.
      */
      uint elem_type_count[MAX_ELEM_TYPES] = {0};
      uint connIndex = 0;

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        cg_npe( Elements[connIndex], &node_count );

        if ( node_count == 3 ) elem_type_count[ pTRI_3  ] += 1;
        if ( node_count == 4 ) elem_type_count[ pQUAD_4 ] += 1;

        connIndex += node_count +1;
      }

      /* create temporary array for sorted indexes. */

      short zone_num_types = 0;
      uint *sorted_tri, *sorted_quad;
      sorted_tri = sorted_quad = NULL;

      if ( elem_type_count[pTRI_3] ) {
        sorted_tri = (uint*)malloc( elem_type_count[pTRI_3]*sizeof(uint) );
        zone_num_types += 1;
      }
      if ( elem_type_count[pQUAD_4] ) {
        sorted_quad = (uint*)malloc( elem_type_count[pQUAD_4]*sizeof(uint) );
        zone_num_types += 1;
      }

      /* setup the zones mixed struct data. */
      mesh->zones[secid].num_types = zone_num_types;
      mesh->zones[secid].type_list = (ElemTypeList*)malloc( zone_num_types*sizeof(ElemTypeList) );

      /* initialize element type indexes. */
      uint idx_tri, idx_quad;
      idx_tri = idx_quad = 0;

      /* sort out the zone element indexes by type. */
      connIndex = 0;
      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {
        cg_npe( Elements[connIndex], &node_count );
        /* Note that we use the idx and then increment it. */
        if ( node_count == 3 ) sorted_tri [ idx_tri++  ] = connIndex;
        if ( node_count == 4 ) sorted_quad[ idx_quad++ ] = connIndex;
        connIndex += node_count +1;
      }

      /* Now we combine the various sorted indexes array
       * into a single one. We also setup the type_list struct
       * for this zone.
       */
      uint* sorted_index = (uint*)malloc( mesh->zones[secid].num_elems*sizeof(uint) );
      if ( sorted_index == NULL ) {
        printf("\nError: could not allocate memory for sorted_index array.\n");
        exit(1);
      }
      uint idx_sorted = 0;
      short typeid = 0;

      if ( elem_type_count[pTRI_3] ) {
        mesh->sizes.elem_type_count[pTRI_3] += elem_type_count[pTRI_3];
        /* setup type_list struct */
        mesh->zones[secid].type_list[typeid].elem_type = pTRI_3;
        mesh->zones[secid].type_list[typeid].start     = zone_type_start;
        mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pTRI_3];
        typeid++;
        /* copy memory block to sorted index */
        memcpy( &sorted_index[idx_sorted], &sorted_tri[0], elem_type_count[pTRI_3]*sizeof(uint) );
        idx_sorted += elem_type_count[pTRI_3];
        zone_type_start += elem_type_count[pTRI_3];
        free(sorted_tri);
      }
      if ( elem_type_count[pQUAD_4] ) {
        mesh->sizes.elem_type_count[pQUAD_4] += elem_type_count[pQUAD_4];
        /* setup type_list struct */
        mesh->zones[secid].type_list[typeid].elem_type = pQUAD_4;
        mesh->zones[secid].type_list[typeid].start     = zone_type_start;
        mesh->zones[secid].type_list[typeid].num_elems = elem_type_count[pQUAD_4];
        typeid++;
        /* copy memory block to sorted index */
        memcpy( &sorted_index[idx_sorted], &sorted_quad[0], elem_type_count[pQUAD_4]*sizeof(uint) );
        idx_sorted += elem_type_count[pQUAD_4];
        zone_type_start += elem_type_count[pQUAD_4];
        free(sorted_quad);
      }

      /* Finally, lets read the elements connectivity sorted by type. */

      uint start = mesh->zones[secid].start;

      for ( uint i = 0; i < mesh->zones[secid].num_elems; i++ ) {

        /* get current element node count/type. */
        cg_npe( Elements[ sorted_index[i] ], &node_count );

        for ( int j = 1; j <= node_count; j++ ) {
          /* we need to offset the index, since cgns starts at 1. */
          mesh->mBFaces[start+i].nodes[j-1] = Elements[ sorted_index[i]+j ] -1;
        }

      }

      /* free temporary arrays. */
      free(sorted_index);

    } /* homogeneous */

    free( Elements );

    zid++;

  } /* num_boundary_zones for loop */

}
