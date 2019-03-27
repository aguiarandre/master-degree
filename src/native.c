#include "native.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void help_datproc(void) {
  printf("\nHELP: \n");
  printf(" Description:\n");
  printf(" Data processing program.\n");
  printf("          argv[0] = \"executable program file name\"\n");
  printf("          argv[1] = \"operation mode\"\n");
  printf("\n");
  printf(" Program operation mode. Valid values are:\n");
  printf("          0 - Convert CGNS mesh to native format, save to disk (hdf5) \n");
  printf("              argv[2] input data file name. \n");
  printf("             	argv[3] output data file name. \n");
  printf("          1 - Convert native (hdf5) mesh and solution to CGNS format \n");
  printf("            	argv[2] input mesh file name.\n");
  printf("             	argv[3] solution data file name.\n");
  printf("             	argv[4] output data file name.\n");
  printf("          2 - Convert CGNS to BRU3D format\n");
  printf("             	argv[2] input data file name. \n");
  printf("             	argv[3] output data file name.\n");
  printf("          3 - Convert linear mesh to curved mesh.\n");
  printf("             	argv[2] input mesh file name (.cgns/hdf5) \n");
  printf("             	argv[3] input CAD file name (.iges)\n");
  printf("\n");
  printf("Usage:\n");
  printf("Convert CGNS to HDF5  : ./datproc 0 ../cgns/cgnsFileName outputHDF5Name\n");
  printf("Convert HDF5 to CGNS  : ./datproc 1 inputHDF5Name ../cgns/cgnsFileName\n");
  printf("Convert CGNS to BRU3D : ./datproc 2 inputHDF5Name ../cgns/cgnsFileName\n");
  printf("Convert to Curved Mesh: ./datproc 3 ../cgns/cgnsFileName ../iges/igesFileName\n");
  printf("\n");

}

void mesh_init(MeshStruct* mesh)

{
  mesh->idx = UNINIT; 						   // Null integers (some start at 0, hence -1)
  mesh->dim = UNINIT;
  
  mesh->sizes.nodes  = UNINIT;
  mesh->sizes.elems  = UNINIT;
  mesh->sizes.faces  = UNINIT;
  mesh->sizes.bfaces = UNINIT;
  mesh->sizes.zones  = UNINIT;
  
  for ( short i = 0; i < MAX_ELEM_TYPES; i++ ) // MAX_ELEM_TYPES = 20 (config.h)
    mesh->sizes.elem_type_count[i] = 0; 	   // counter of element types in the mesh.	

  mesh->elem_type = pNONE;					   // No type defined.
  mesh->zones     = NULL;					   
  mesh->vertices  = NULL;	// Initialize pointer to vertice
  mesh->mFaces    = NULL;
  mesh->mBFaces   = NULL;
  mesh->mElems    = NULL;
  
  mesh->lookUpArray = malloc( (MESH_ORDER-1) * sizeof(double) );
  switch( MESH_ORDER )
  {
	  case 2:
	  mesh->lookUpArray[0] = 0.5;
	  break;
	  
	  case 3:
	  mesh->lookUpArray[0] = 1.0/3.0;
	  mesh->lookUpArray[1] = 2.0/3.0;
	  break;
	  
	  case 4:
	  mesh->lookUpArray[0] = 1.0/4.0;
	  mesh->lookUpArray[1] = 2.0/4.0;
	  mesh->lookUpArray[2] = 3.0/4.0;
	  break;
	  
	  default:
	  printf("Not done for mesh order greater than 4\n");
	  exit(1);
	  break;
  }
}


/* free memory allocated with malloc */
void mesh_finish(MeshStruct* mesh)
{ 
  /* free allocated mesh zones */
  for ( int i = 0; i < mesh->sizes.zones; i++ ) 				// Loop through zones
    mesh_zone_finish(&mesh->zones[i]);

  if ( mesh->zones ) free( mesh->zones );

  /* free allocated mesh vertices */
 
  if ( mesh->vertices ) {
    for (uint i = 0; i < mesh->sizes.nodes; i++){ 				// Loop through nodes
      free( mesh->vertices[i] );
    }
    free( mesh->vertices );
  }

  /* free allocated mesh element mappings */
  if ( mesh->mElems ) {
    for (uint i = 0; i < mesh->sizes.elems; i++){				// Loop through elements
      if( mesh->mElems[i].nodes ) free( mesh->mElems[i].nodes );
      if( mesh->mElems[i].nbs   ) free( mesh->mElems[i].nbs   );
      if( mesh->mElems[i].faces ) free( mesh->mElems[i].faces );
    }
    free( mesh->mElems );
  }

  /* free allocated mesh face mappings */
  if ( mesh->mFaces ) {
    for (uint i = 0; i < mesh->sizes.faces; i++){				// Loop through faces
      if( mesh->mFaces[i].nodes ) free( mesh->mFaces[i].nodes );
    }
    free( mesh->mFaces );
  }

  /* free allocated mesh boundary face mappings */
  if ( mesh->mBFaces ) {
    for (uint i = 0; i < mesh->sizes.bfaces; i++){				// Loop through boundary faces
      if( mesh->mBFaces[i].nodes ) free( mesh->mBFaces[i].nodes );
    }
    free( mesh->mBFaces );
  }
  
  //free( mesh->lookUpArray );

}




void mesh_zone_init(MeshZone* zone)
{
  zone->idx = UNINIT;
  zone->name = NULL;
  zone->dim = UNINIT;
  zone->is_boundary = false;
  zone->num_elems = UNINIT;
  zone->elem_type = pNONE;
  zone->type_list = NULL;
  zone->homogeneous = false;
  zone->start = UNINIT;
  zone->bc_data = NULL;
}


void mesh_zone_finish(MeshZone* zone)
{
  if( zone->name      ) free( zone->name      );
  if( zone->type_list ) free( zone->type_list );
  if( zone->bc_data   ) free( zone->bc_data   );
}

/* prints data about mesh. */
void mesh_report(MeshStruct* mesh)
{
  /* print mesh size information. */
  printf("\nMesh information:");
  printf("\n dimension = %d", mesh->dim );
  printf("\n elements  = %d", mesh->sizes.elems );
  printf("\n nodes     = %d", mesh->sizes.nodes );
  printf("\n faces     = %d", mesh->sizes.faces );
  printf("\n bfaces    = %d", mesh->sizes.bfaces);
  printf("\n zones     = %d", mesh->sizes.zones );

  printf("\nZone information:");
  for (uint i = 0; i < mesh->sizes.zones; i++) {	// Loop through number of zones.
    /* get zone element type. */
    char elem_type[10];
    if ( mesh->zones[i].homogeneous ) {
      if ( mesh->zones[i].elem_type == pBAR_2   ) sprintf(elem_type, "BAR_2"  );
      if ( mesh->zones[i].elem_type == pTRI_3   ) sprintf(elem_type, "TRI_3"  );
    //  if ( mesh->zones[i].elem_type == pTRI_6   ) sprintf(elem_type, "TRI_6"  );
      if ( mesh->zones[i].elem_type == pQUAD_4  ) sprintf(elem_type, "QUAD_4" );
      if ( mesh->zones[i].elem_type == pTETRA_4 ) sprintf(elem_type, "TETRA_4");
      if ( mesh->zones[i].elem_type == pPYRA_5  ) sprintf(elem_type, "PYRA_5" );
      if ( mesh->zones[i].elem_type == pPENTA_6 ) sprintf(elem_type, "PENTA_6");
      if ( mesh->zones[i].elem_type == pHEXA_8  ) sprintf(elem_type, "HEXA_8" );
    } else {
      sprintf(elem_type, "MIXED");
    }

    /* report data. */
    printf("\n type = %8s, size = %8d, elem_type = %7s, name = %s",
            mesh->zones[i].is_boundary ? "boundary" : "interior",
            mesh->zones[i].num_elems, elem_type, mesh->zones[i].name);
  }
  printf("\nDone.\n");
  fflush(stdout);
}


uint mesh_compute_faces_number(MeshStruct* mesh)
{
  uint nfaces  = 0,
       nbfaces = 0;

  nbfaces = mesh->sizes.bfaces;

  nfaces = mesh->sizes.elem_type_count[ pTRI_3   ]*3
         + mesh->sizes.elem_type_count[ pQUAD_4  ]*4
         + mesh->sizes.elem_type_count[ pTETRA_4 ]*4
         + mesh->sizes.elem_type_count[ pPYRA_5  ]*5
         + mesh->sizes.elem_type_count[ pPENTA_6 ]*5
         + mesh->sizes.elem_type_count[ pHEXA_8  ]*6;

  return nfaces = (nfaces - nbfaces)/2 + nbfaces;
}

void mesh_create_faces(MeshStruct* mesh)
{
  /* by now, the mesh struct has the necessary data to
   * compute faces connectivity and neighbors.
   */
  printf("\nComputing mesh face data...");
  fflush(stdout);

  /* get number of faces. */
  uint nfaces = mesh_compute_faces_number(mesh);

  /* estimate the hash table size. */
  uint hashspace = utils_next_prime( mesh->sizes.elems );

  /* update the estimate to handle hash size. */
  nfaces *= 2;

  /* --------------------- allocate structs --------------------- */

  /* Temporary storage for faces creation. 
   * Note that this array begins at face counter 1, instead of
   * zero, as in the original Fortran algorithm.
   */
  MapFaces* facestemp = (MapFaces*)malloc( nfaces*sizeof(MapFaces) );
  if ( facestemp ) {
    for( uint i = 0; i < nfaces; i++ ) {
      facestemp[i].nodes = (uint*)malloc( MAX_FACE_NODES*sizeof(uint) );
      for( short j = 0; j < MAX_FACE_NODES; j++ ) {
        facestemp[i].nodes[j] = UNINIT;
      }
      facestemp[i].elem_left  = UNINIT;
      facestemp[i].elem_right = UNINIT;
    }
  }
  
   /* Allocate mesh internal faces map. */
  //mesh->mElems = (MapElems*)malloc( mesh->sizes.elems*sizeof(MapElems) );
  if ( mesh->mElems ) {
    for( uint i = 0; i < mesh->sizes.elems; i++ )
      mesh->mElems[i].faces = (uint*)malloc( 8*sizeof(uint) );
  }
  else if ( mesh->mElems == NULL ) {
    printf("\nError: could not allocate memory for mesh faces struct.\n");
    exit(1);  
  }
  
  
  
  else if ( facestemp == NULL ) {
    printf("\nError: could not allocate memory for facestemp struct.\n");
    exit(1);
  }

  /* allocate and initialize arrays. */
  uint* hash = (uint*)malloc( nfaces*sizeof(uint) );
  if ( hash == NULL ) {
    printf("\nError: could not allocate memory for hash.\n");
    exit(1);
  }

  uint* iprox = (uint*)malloc( nfaces*sizeof(uint) );
  if ( iprox == NULL ) {
    printf("\nError: could not allocate memory for iprox.\n");
    exit(1);
  }

  for ( uint i = 0; i < nfaces; i++ ) {
    hash[i] = UNINIT;
    iprox[i] = UNINIT;
  }

  /* --------------------- create mesh faces --------------------- */

  /* compute the number of digitis in mesh->sizes.nodes,
   * used latter by the hash key function. */
  short ipr = floor(log10 ( mesh->sizes.nodes) ) + 1;

  uint icol; /* collisions counter. */
  float pct = 0.0f;
  uint itl;
  uint ipont;
  uint cont, aux;
  icol = itl = 0;
  uint key;

  /* temporary array to hold an element faces connectivity. */
  uint** facends;
  facends = (uint**)malloc(MAX_FACE_NODES*sizeof(uint*));
  for(short i = 0; i < MAX_FACE_NODES; i++) {
    facends[i] = (uint*)malloc(MAX_ELEM_NBS*sizeof(uint));
    for ( short j = 0; j < MAX_ELEM_NBS; j++ ) {
      facends[i][j] = UNINIT;
    }
  }
	
	nfaces = 0;
  /* reset the faces counter. */
	
	
  /* create and add internal mesh face keys into the hash table. */
  for( uint id = 0; id < mesh->sizes.zones; id++ ) {

    if ( ! mesh->zones[id].is_boundary ) {

      uint zone_start    = mesh->zones[id].start;
      uint zone_end      = zone_start + mesh->zones[id].num_elems;

      /* loop over this zone elements. */
      for( uint elem = zone_start; elem < zone_end; elem++ ) {
		  
	
		
        short num_faces = 0;
       
        get_elem_faces( mesh, id, elem, &num_faces, facends );

        /* loop over this element faces. */
        for( short faceid = 0; faceid < num_faces; faceid++ ) {

          /* create the hash key for this face. */
          uint an, bn, cn, dn;
          an = bn = cn = dn = UNINIT;

          an = facends[0][faceid];
          bn = facends[1][faceid];
          
          if ( mesh->dim == 3 ) {
            cn = facends[2][faceid];
            dn = facends[3][faceid];
          }

          key = get_hash_key( mesh->dim, ipr, hashspace, &an, &bn, &cn, &dn );

          if ( key > hashspace ) {
            printf("\nError: increase hashspace value.\n");
            exit(1);
          }

          if ( hash[key] == UNINIT ) {
            facestemp[nfaces].nodes[0] = facends[0][faceid];
            facestemp[nfaces].nodes[1] = facends[1][faceid];
            if ( mesh->dim == 3 ) {
              facestemp[nfaces].nodes[2] = facends[2][faceid];
              facestemp[nfaces].nodes[3] = facends[3][faceid];
            }
            /* list face neighbors. The right cell is defined latter. */
            facestemp[nfaces].elem_left  = elem;   
            
            
            /* set up hash variables. */
            iprox[nfaces] = UNINIT;
            hash[key] = nfaces;
            
            
           // hash_temp[key] = nfaces;
          
         //   mesh->mElems[elem].faces[faceid] = hash[key];  // setting up mesh->mElems[elem].faces
			nfaces++;
           
          }
          else {

		//	mesh->mElems[elem].faces[faceid] = 0;
            /* face might already exist. */

            uint a, b, c, d;
            a = b = c = d = UNINIT;

            a = facestemp[ hash[key] ].nodes[0];
            b = facestemp[ hash[key] ].nodes[1];
            if ( mesh->dim == 3 ) {
              c = facestemp[ hash[key] ].nodes[2];
              d = facestemp[ hash[key] ].nodes[3];
            }
            key = get_hash_key( mesh->dim, ipr, hashspace, &a, &b, &c, &d );
			
			
            if ( a == an && b == bn && c == cn && d == dn ) {
              facestemp[ hash[key] ].elem_right = elem;
            }
            else {
				
			
              /* we got a collision. */

              cont = 1;
              icol += 1;
              ipont = iprox[ hash[key] ];
              aux = hash[key];

label_10:
              if ( ipont == UNINIT ) {
                facestemp[nfaces].nodes[0] = facends[0][faceid];
                facestemp[nfaces].nodes[1] = facends[1][faceid];
                if ( mesh->dim == 3 ) {
                  facestemp[nfaces].nodes[2] = facends[2][faceid];
                  facestemp[nfaces].nodes[3] = facends[3][faceid];
                }

                facestemp[nfaces].elem_left  = elem;
                facestemp[nfaces].elem_right = UNINIT;
				
				
				
                /* set up hash variables. */
                iprox[nfaces] = UNINIT;

                /* update pointer to the before-last position
                 * on the collision table. */
                iprox[aux] = nfaces;
				
                nfaces++;

                /* hash efficiency variables. */
                pct += ( (float)cont / ((float)aux+1.f) )*100;
                itl++;
              }
              else {
                /* include neighbor information to existing face. */
                uint a, b, c, d;
                a = b = c = d = UNINIT;

                a = facestemp[ ipont ].nodes[0];
                b = facestemp[ ipont ].nodes[1];
                if ( mesh->dim == 3 ) {
                  c = facestemp[ ipont ].nodes[2];
                  d = facestemp[ ipont ].nodes[3];
                }

                key = get_hash_key( mesh->dim, ipr, hashspace, &a, &b, &c, &d );

                if ( a == an && b == bn && c == cn && d == dn ) {
                  facestemp[ ipont ].elem_right = elem;
                  icol--;
                  pct += ( (float)cont / ((float)aux+1.f) )*100;
                  itl++;
                  
                }
                else {
                  /* search for the corresponding side on the collision table. */
                  icol++;
                  aux = ipont;
                  ipont = iprox[ipont];
                  cont++;
                  goto label_10;
                }
              }
            }
          }

        } /* faceid for loop */

      } /* elem for loop */

    } /* is_boundary if */

  } /* zone for loop */


  mesh->sizes.faces = nfaces -1;

  float efc;
  if ( itl == 0 ) efc = 100.0;
  if ( itl != 0 ) efc = 100.0 - pct/(float)itl;

  /* output hash information. */
  printf("\nnumber of faces      = %d", mesh->sizes.faces);
  printf("\nnumber of collisions = %d", icol);
  printf("\nhash efficiency      = %6.2f%%", efc);
  fflush(stdout);

  /* ---------------- setup boundary faces neighbors ---------------- */

  for ( uint bfaceid = 0; bfaceid < mesh->sizes.bfaces; bfaceid++ ) {

    /* create the hash key for this face. */

    uint an, bn, cn, dn;
    an = bn = cn = dn = UNINIT;

    an = mesh->mBFaces[bfaceid].nodes[0];
    bn = mesh->mBFaces[bfaceid].nodes[1];
    if ( mesh->dim == 3 ) {
      cn = mesh->mBFaces[bfaceid].nodes[2];
      dn = mesh->mBFaces[bfaceid].nodes[3];
    }
    key = get_hash_key( mesh->dim, ipr, hashspace, &an, &bn, &cn, &dn );

    if ( hash[key] == UNINIT ) {
      /* we have a problem. At this stage, all internal faces should exist. */
      printf("\nError: bface creation 1.\n");
      exit(1);
    }
    else {
      uint a, b, c, d;
      a = b = c = d = UNINIT;

      a = facestemp[ hash[key] ].nodes[0];
      b = facestemp[ hash[key] ].nodes[1];
      if ( mesh->dim == 3 ) {
        c = facestemp[ hash[key] ].nodes[2];
        d = facestemp[ hash[key] ].nodes[3];
      }

      key = get_hash_key( mesh->dim, ipr, hashspace, &a, &b, &c, &d );

      if ( a == an && b == bn && c == cn && d == dn ) {
        /* Note that we store the face nodes in a way that the
         * counter-clokwise ordering creates an outward normal to
         * the internal neighbor.
         */
        mesh->mBFaces[bfaceid].elem_left = facestemp[ hash[key] ].elem_left;
      }
      else {
        /* we have a collision. */
        ipont = iprox[ hash[key] ];
label_20:
        if ( ipont == UNINIT ) {
          /* Problem, the routine did not find a face. */
          printf("\nError: bface creation 2.\n");
          exit(1);
        }
        else {
          /* include neighbor information to existing face. */
          uint a, b, c, d;
          a = b = c = d = UNINIT;

          a = facestemp[ipont].nodes[0];
          b = facestemp[ipont].nodes[1];
          if ( mesh->dim == 3 ) {
            c = facestemp[ipont].nodes[2];
            d = facestemp[ipont].nodes[3];
          }

          key = get_hash_key( mesh->dim, ipr, hashspace, &a, &b, &c, &d );

          if ( a == an && b == bn && c == cn && d == dn ) {
            mesh->mBFaces[bfaceid].elem_left = facestemp[ipont].elem_left;
          }
          else {
            /* search for the corresponding side on the collision table. */
            ipont = iprox[ ipont ];
            goto label_20;
          }
        }
      }
    } /* hash[key] == 0 if */

  } /* bfaceid for loop */

  /* --------------------- set mesh face struct --------------------- */

  /* free temporary data. */
 // free( hash );
  free( iprox );
 // for(short i = 0; i < MAX_FACE_NODES; i++)
 //   free( facends[i] );
 // free(facends);

  /* remove boundary faces from the internal face set. */
  uint int_faces = 0;
  for( uint i = 0; i < mesh->sizes.faces; i++ ) {
    if ( facestemp[i].elem_right != UNINIT ) int_faces++;
  }

  /* Allocate mesh internal faces map. */
  mesh->mFaces = (MapFaces*)malloc( int_faces*sizeof(MapFaces) );
  if ( mesh->mFaces ) {
    for( uint i = 0; i < int_faces; i++ )
      mesh->mFaces[i].nodes = (uint*)malloc( MAX_FACE_NODES*sizeof(uint) );
  }
  else if ( mesh->mFaces == NULL ) {
    printf("\nError: could not allocate memory for mesh faces struct.\n");
    exit(1);
  }
 
  
   
  
  uint k = 0;
 
  for(uint i = 0; i < mesh->sizes.faces; i++ ) {
    if ( facestemp[i].elem_right != UNINIT ) {
      for (short j = 0; j < MAX_FACE_NODES; j++ ) {
        mesh->mFaces[k].nodes[j] = facestemp[i].nodes[j];
      }
      
     
	  
      mesh->mFaces[k].elem_left  = facestemp[i].elem_left;
      mesh->mFaces[k].elem_right = facestemp[i].elem_right;
      
      k++;
     
    }
    
  }

    
	uint face_count=0, elem;

	
  /* adjust mesh face size. */
  mesh->sizes.faces = int_faces;


//  ipr = floor(log10 (abs (mesh->sizes.nodes))) + 1;







  //--------------------------------------------------------------------------------------------------------------//

// Run through internal elements.

uint correct_faceid;

for (uint i = 0; i < mesh->sizes.faces; i++) {


			// I am at face 'i'. 
	
		elem = mesh->mFaces[i].elem_left; // Now I'm looking to the left element.
	for (short zoneid = 0; zoneid < mesh->sizes.zones; zoneid++) {	

		uint a, b, c, d;
	
	    a = b = c = d = UNINIT;

		a = mesh->mFaces[i].nodes[0];
		b = mesh->mFaces[i].nodes[1];
		if ( mesh->dim == 3 ) {
	      c = mesh->mFaces[i].nodes[2];
	      d = mesh->mFaces[i].nodes[3];
	    }  // Face nodes
	
		short num_faces = 0;
		get_elem_faces( mesh, zoneid, elem, &num_faces, facends ); // Need to know all of its nodes.
			
			for (short faceid = 0; faceid < num_faces; faceid++) {	
			uint an, bn, cn, dn;
			
			    an = bn = cn = dn = UNINIT;
			
			    an = facends[0][faceid];
			    bn = facends[1][faceid];
			    if ( mesh->dim == 3 ) {
			      cn = facends[2][faceid];
			      dn = facends[3][faceid];
			    } // Now I have this element nodes.
			
			//Now I need to compare the nodes I have (a,b,c,d) with the nodes I got in this loop (an, bn, cn, dn). And then I can find the faceid I need to use.
			
				if ( (an == a && bn == b && cn == c && dn == d) || (an == b && bn == a && cn == d && dn == c) ) { // Then we are at the correct faceid
					correct_faceid = faceid;
					face_count++;
				}
			}	
			
			mesh->mElems[elem].faces[correct_faceid] = i;
	}
// ----------------------------------- same thing for the right element, remember it is inverted. -----------------------------------
	
		elem = mesh->mFaces[i].elem_right;

for (short zoneid = 0; zoneid < mesh->sizes.zones; zoneid++) {	
		uint a, b, c, d;
	    a = b = c = d = UNINIT;

		a = mesh->mFaces[i].nodes[0];
		b = mesh->mFaces[i].nodes[1];
		if ( mesh->dim == 3 ) {
	      c = mesh->mFaces[i].nodes[2];
	      d = mesh->mFaces[i].nodes[3];
	    } 
	    
		short num_faces = 0;
		get_elem_faces( mesh, zoneid, elem, &num_faces, facends ); // Need to know all of its nodes.
			
			for (short faceid = 0; faceid < num_faces; faceid++) {	
				uint an, bn, cn, dn;
			
			    an = bn = cn = dn = UNINIT;
			
			    an = facends[0][faceid];
			    bn = facends[1][faceid];
			    if ( mesh->dim == 3 ) {
			      cn = facends[2][faceid];
			      dn = facends[3][faceid];
			    } // Now I have this element nodes.
			
			//Now I need to compare the nodes I have with the nodes I got in this loop. And then I can find the faceid I need to use.
			
				if ( (an == b && bn == a && cn == d && dn == c) || (an == a && bn == b && cn == c && dn == d) ) { // Then we are at the correct faceid
					correct_faceid = faceid;
			
					face_count++;
				}
			}	
			
			mesh->mElems[elem].faces[correct_faceid] = i;
	}

}

for (uint i = 0; i < mesh->sizes.bfaces; i++) {

		
	// ----------------------------------- same thing for the boundary element -----------------------------------//
	
		elem = mesh->mBFaces[i].elem_left;
		
	for (short zoneid = 0; zoneid < mesh->sizes.zones; zoneid++) {	
		uint a,b,c,d;
		
	    a = b = c = d = UNINIT;

		a = mesh->mBFaces[i].nodes[0];
		b = mesh->mBFaces[i].nodes[1];
		if ( mesh->dim == 3 ) {
	      c = mesh->mBFaces[i].nodes[2];
	      d = mesh->mBFaces[i].nodes[3];
	    } 
	    
		short num_faces = 0;
		get_elem_faces( mesh, zoneid, elem, &num_faces, facends ); // Need to know all of its nodes.
			
			for (short faceid = 0; faceid < num_faces; faceid++) {	
				uint an, bn, cn, dn;
			
			    an = bn = cn = dn = UNINIT;
			
			    an = facends[0][faceid];
			    bn = facends[1][faceid];
			    if ( mesh->dim == 3 ) {
			      cn = facends[2][faceid];
			      dn = facends[3][faceid];
			    } // Now I have this element nodes.
			
			//Now I need to compare the nodes I have with the nodes I got in this loop. And then I can find the faceid I need to use.
			
				if ( (an == a && bn == b && cn == c && dn == d) || (an == b && bn == a && cn == d && dn == c) ) { // Then we are at the correct faceid
					correct_faceid = faceid;
					face_count++;
				}
			}	
			
			mesh->mElems[elem].faces[correct_faceid] = i;
	}
}
printf("\nface count = %d\n", face_count);

    //--------------------------------------------------------------------------------------------------------------//
  
  
  

  

  
  


 


  /* finish deallocation. */
  for( short i = 0; i < MAX_FACE_NODES; i++ )
    free( facestemp[i].nodes );
  free( facestemp );

  printf("\nDone.\n");
  fflush(stdout);
}

void get_elem_faces( MeshStruct* mesh, uint zoneid, uint elem,
                     short* num_faces, uint** nodes)
{
  /* Define the element faces based on its nodes. */
  bool homogeneous = mesh->zones[zoneid].homogeneous;

  /* figure out this element type. */
  ElemType type;
  if ( homogeneous ) {
    type = mesh->zones[zoneid].elem_type;
  } else {
    for( short i = 0; i < mesh->zones[zoneid].num_types; i++ ) {
      uint start = mesh->zones[zoneid].type_list[i].start;
      uint end   = start + mesh->zones[zoneid].type_list[i].num_elems -1;
      if ( elem >= start && elem <= end ) {
        type = mesh->zones[zoneid].type_list[i].elem_type;
      }
    }
  }

  if ( type == pTRI_3 ) {
    *num_faces = 3;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[0];
    nodes[1][0] = mesh->mElems[elem].nodes[1];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[1];
    nodes[1][1] = mesh->mElems[elem].nodes[2];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[2];
    nodes[1][2] = mesh->mElems[elem].nodes[0];
  }
  
    if ( type == pTRI_6 ) {
    *num_faces = 3;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[0];	// nó 1 da face 1 
    nodes[1][0] = mesh->mElems[elem].nodes[1];	// nó 2 da face 1 
    
    nodes[1][0] = mesh->mElems[elem].nodes[1];	// nó 2 da face 1
    nodes[2][0] = mesh->mElems[elem].nodes[2];	// nó 3 da face 1
  
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[1];	//
    nodes[1][1] = mesh->mElems[elem].nodes[2];
    
    nodes[0][1] = mesh->mElems[elem].nodes[1];
    nodes[1][1] = mesh->mElems[elem].nodes[2];
    
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[2];
    nodes[1][2] = mesh->mElems[elem].nodes[0];
  
	nodes[0][2] = mesh->mElems[elem].nodes[2];
    nodes[1][2] = mesh->mElems[elem].nodes[0];
  }
  
  
  if ( type == pQUAD_4) {
    *num_faces = 4;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[0];
    nodes[1][0] = mesh->mElems[elem].nodes[1];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[1];
    nodes[1][1] = mesh->mElems[elem].nodes[2];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[2];
    nodes[1][2] = mesh->mElems[elem].nodes[3];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[3];
    nodes[1][3] = mesh->mElems[elem].nodes[0];
    
    //nodes[# of points in a face][# of faces]
  }

  if ( type == pQUAD_9) {
    *num_faces = 8;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[0];
    nodes[1][0] = mesh->mElems[elem].nodes[4];
    
    /*second face */
    nodes[0][1] = mesh->mElems[elem].nodes[4];
    nodes[1][1] = mesh->mElems[elem].nodes[1];
    
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[1];
    nodes[1][2] = mesh->mElems[elem].nodes[5];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[5];
    nodes[1][3] = mesh->mElems[elem].nodes[2];
    /* fifth face */
    nodes[0][4] = mesh->mElems[elem].nodes[2];
    nodes[1][4] = mesh->mElems[elem].nodes[6];
 
    /* sixth face */
    nodes[0][5] = mesh->mElems[elem].nodes[6];
    nodes[1][5] = mesh->mElems[elem].nodes[3];

    /* seventh face */
    nodes[0][6] = mesh->mElems[elem].nodes[3];
    nodes[1][6] = mesh->mElems[elem].nodes[7];

    /* eighth face */
    nodes[0][7] = mesh->mElems[elem].nodes[7];
    nodes[1][7] = mesh->mElems[elem].nodes[0];

  }

  if ( type == pTETRA_4 ) {
    *num_faces = 4;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[3];
    nodes[1][0] = mesh->mElems[elem].nodes[1];
    nodes[2][0] = mesh->mElems[elem].nodes[2];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[3];
    nodes[1][1] = mesh->mElems[elem].nodes[2];
    nodes[2][1] = mesh->mElems[elem].nodes[0];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[3];
    nodes[1][2] = mesh->mElems[elem].nodes[0];
    nodes[2][2] = mesh->mElems[elem].nodes[1];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[1];
    nodes[1][3] = mesh->mElems[elem].nodes[0];
    nodes[2][3] = mesh->mElems[elem].nodes[2];
  }
  if ( type == pPYRA_5 ) {
    *num_faces = 5;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[1];
    nodes[1][0] = mesh->mElems[elem].nodes[2];
    nodes[2][0] = mesh->mElems[elem].nodes[4];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[3];
    nodes[1][1] = mesh->mElems[elem].nodes[4];
    nodes[2][1] = mesh->mElems[elem].nodes[2];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[3];
    nodes[1][2] = mesh->mElems[elem].nodes[0];
    nodes[2][2] = mesh->mElems[elem].nodes[4];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[0];
    nodes[1][3] = mesh->mElems[elem].nodes[1];
    nodes[2][3] = mesh->mElems[elem].nodes[4];
    /* fifth face */
    nodes[0][4] = mesh->mElems[elem].nodes[0];
    nodes[1][4] = mesh->mElems[elem].nodes[3];
    nodes[2][4] = mesh->mElems[elem].nodes[2];
    nodes[3][4] = mesh->mElems[elem].nodes[1];
  }
  if ( type == pPENTA_6 ) {
    *num_faces = 5;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[1];
    nodes[1][0] = mesh->mElems[elem].nodes[2];
    nodes[2][0] = mesh->mElems[elem].nodes[5];
    nodes[3][0] = mesh->mElems[elem].nodes[4];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[0];
    nodes[1][1] = mesh->mElems[elem].nodes[3];
    nodes[2][1] = mesh->mElems[elem].nodes[5];
    nodes[3][1] = mesh->mElems[elem].nodes[2];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[0];
    nodes[1][2] = mesh->mElems[elem].nodes[1];
    nodes[2][2] = mesh->mElems[elem].nodes[4];
    nodes[3][2] = mesh->mElems[elem].nodes[3];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[0];
    nodes[1][3] = mesh->mElems[elem].nodes[2];
    nodes[2][3] = mesh->mElems[elem].nodes[1];
    /* fifth face */
    nodes[0][4] = mesh->mElems[elem].nodes[3];
    nodes[1][4] = mesh->mElems[elem].nodes[4];
    nodes[2][4] = mesh->mElems[elem].nodes[5];
  }
  if ( type == pHEXA_8 ) {
    *num_faces = 8;
    /* first face */
    nodes[0][0] = mesh->mElems[elem].nodes[1];
    nodes[1][0] = mesh->mElems[elem].nodes[2];
    nodes[2][0] = mesh->mElems[elem].nodes[6];
    nodes[3][0] = mesh->mElems[elem].nodes[5];
    /* second face */
    nodes[0][1] = mesh->mElems[elem].nodes[3];
    nodes[1][1] = mesh->mElems[elem].nodes[7];
    nodes[2][1] = mesh->mElems[elem].nodes[6];
    nodes[3][1] = mesh->mElems[elem].nodes[2];
    /* third face */
    nodes[0][2] = mesh->mElems[elem].nodes[0];
    nodes[1][2] = mesh->mElems[elem].nodes[4];
    nodes[2][2] = mesh->mElems[elem].nodes[7];
    nodes[3][2] = mesh->mElems[elem].nodes[3];
    /* fourth face */
    nodes[0][3] = mesh->mElems[elem].nodes[0];
    nodes[1][3] = mesh->mElems[elem].nodes[1];
    nodes[2][3] = mesh->mElems[elem].nodes[5];
    nodes[3][3] = mesh->mElems[elem].nodes[4];
    /* fifth face */
    nodes[0][4] = mesh->mElems[elem].nodes[0];
    nodes[1][4] = mesh->mElems[elem].nodes[3];
    nodes[2][4] = mesh->mElems[elem].nodes[2];
    nodes[3][4] = mesh->mElems[elem].nodes[1];
    /* sixth face */
    nodes[0][5] = mesh->mElems[elem].nodes[4];
    nodes[1][5] = mesh->mElems[elem].nodes[5];
    nodes[2][5] = mesh->mElems[elem].nodes[6];
    nodes[3][5] = mesh->mElems[elem].nodes[7];
  }
}

/* sort the nodes with ascending values in
 * preparation for face hash key.
 */
uint get_hash_key( short dim, int ipr, uint hashspace,
                   uint *a, uint *b, uint *c, uint *d )
{
  uint nmax;
  uint nmin;

  if ( dim == 2 )
  {
    int node[2];
    int aux;

    node[0] = *a;
    node[1] = *b;

    if ( node[0] > node[1] ) {
      aux = node[1];
      node[1] = node[0];
      node[0] = aux;
    }

    /* store sorted nodes. */
    *a = node[0];
    *b = node[1];

    nmax = *a + *b;
    nmin = 0;
  }

  if ( dim == 3 )
  {
    int node[4];
    int aux;

    node[0] = *a;
    node[1] = *b;
    node[2] = *c;
    node[3] = *d;

    short node_count = 4;
    if ( *d == UNINIT ) {
      /* we have a tri face. */
      node_count = 3;
    }

    for ( int j = 0; j < node_count; j++ ) {
      for ( int i = 0; i < node_count-1; i++ ) {
        if ( node[i] > node[i+1] ) {
          aux = node[i+1];
          node[i+1] = node[i];
          node[i] = aux;
        }
      }
    }

    /* store sorted nodes. */
    *a = node[0];
    *b = node[1];
    *c = node[2];
    *d = node[3];

    nmax = *a + *b;
    nmin = *c + *d;
  }

  /* return key function. */
  uint key;
  key = nmax*( pow(2,ipr) ) + nmin;
  key = ( key % hashspace ) + 1;

  return key;
}

void mesh_set_bcs(MeshStruct* mesh)
{
	/* Query user about boundary types. */
	
	{
		printf("\n");
		printf("Select the appropriate zone type\n");
		printf("0 - fluid\n");
		printf("1 - wall\n");
		printf("2 - farfield\n");
		printf("\n");
		fflush(stdout);

		char buffer[256];

		for( uint z = 0; z < mesh->sizes.zones; z++ )
		{
			printf("zone: %s, type: ",mesh->zones[z].name);
			fflush(stdout);
			fgets(buffer,256,stdin);
			mesh->zones[z].bc_id = atoi(buffer);
		}
	}

}

