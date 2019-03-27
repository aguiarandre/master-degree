#include "native.h"
#include "native.h"
#include "math.h"

#include <stdlib.h> // malloc

void mesh_export_obj(MeshStruct* mesh)
{
  char* filename = "gui.obj";
  
  FILE* fp = fopen( filename, "w" );
  
  
  if ( mesh->dim == 2 ) {
    
    /* write all mesh vertices */
    double x, y, z;
    for ( uint i = 0; i < mesh->sizes.nodes; i++ ) {
      x = mesh->vertices[i][0];
      y = mesh->vertices[i][1];
      z = 0.0f;
      fprintf( fp, "v %.16g %.16g %.16g\n", x, y, z);
    }
    
    
    /* write mesh wireframe before faces, so there is
     * no z-trashing. */
    fprintf( fp, "o #WIRE_wireframe\n" );
    for ( uint i = 0; i < mesh->sizes.faces; i++ ) {
      uint na = mesh->mFaces[i].nodes[0]+1;
      uint nb = mesh->mFaces[i].nodes[1]+1;
      fprintf( fp, "w %d %d\n", na, nb );
    }
    
    
    /* write elements connectivity */
    for ( uint z = 0; z < mesh->sizes.zones; z++ ) {
      
      /* write this zone name. */
      
      if ( ! mesh->zones[z].is_boundary ) {
        
        char int_zone[STRING_MAX_LENGHT];
        sprintf( int_zone, "#IN2D_%s", mesh->zones[z].name );
        fprintf( fp, "o %s\n", int_zone );
        
        if ( mesh->zones[z].homogeneous ) {
          
          if ( mesh->zones[z].elem_type == pTRI_3 ) {
            
            /* write connectivity direct to the obj file. */
            uint na, nb, nc;
            for ( uint i = 0; i < mesh->sizes.elems; i++ ) {
              na = mesh->mElems[i].nodes[0]+1;
              nb = mesh->mElems[i].nodes[1]+1;
              nc = mesh->mElems[i].nodes[2]+1;
              fprintf( fp, "f %d %d %d\n", na, nb, nc );
            }
          
          }
          else if ( mesh->zones[z].elem_type == pQUAD_4 ) {
            
            uint size = mesh->zones[z].num_elems;
            uint beg  = mesh->zones[z].start;
            uint end  = beg + size;
            
            uint na, nb, nc;
            
            for ( uint i = beg; i < end; i++ ) {
              
              /* create two tris from one quad. */

              na = mesh->mElems[i].nodes[0]+1;
              nb = mesh->mElems[i].nodes[1]+1;
              nc = mesh->mElems[i].nodes[2]+1;
              fprintf( fp, "f %d %d %d\n", na, nb, nc );
              
              na = mesh->mElems[i].nodes[2]+1;
              nb = mesh->mElems[i].nodes[3]+1;
              nc = mesh->mElems[i].nodes[0]+1;
              fprintf( fp, "f %d %d %d\n", na, nb, nc );
              
            }
            
          }
        
        }
        
        else {
          
          /* we are in a mixed zone. */
          
          for ( uint k = 0; k < mesh->zones[z].num_types; k++ ) {
            
            ElemType type = mesh->zones[z].type_list[k].elem_type;
            
            uint size = mesh->zones[z].type_list[k].num_elems;
            uint beg  = mesh->zones[z].type_list[k].start;
            uint end  = beg + size;
            
            uint na, nb, nc;
            
            if ( type == pTRI_3 ) {
              
              for ( uint i = beg; i < end; i++ ) {
                na = mesh->mElems[i].nodes[0]+1;
                nb = mesh->mElems[i].nodes[1]+1;
                nc = mesh->mElems[i].nodes[2]+1;
                fprintf( fp, "f %d %d %d\n", na, nb, nc );
              }
              
            }
            else if ( type == pQUAD_4 ) {
            
              for ( uint i = beg; i < end; i++ ) {
                /* create two tris from one quad. */

                na = mesh->mElems[i].nodes[0]+1;
                nb = mesh->mElems[i].nodes[1]+1;
                nc = mesh->mElems[i].nodes[2]+1;
                fprintf( fp, "f %d %d %d\n", na, nb, nc );
                
                na = mesh->mElems[i].nodes[2]+1;
                nb = mesh->mElems[i].nodes[3]+1;
                nc = mesh->mElems[i].nodes[0]+1;
                fprintf( fp, "f %d %d %d\n", na, nb, nc );
              }
            
            }
          
          }
          
        }
        
      }
      else {
        
        /* we are in a boundary zone. Note that in 2D
           there are only bar elements. */
        
        char bc_zone[STRING_MAX_LENGHT];
        sprintf( bc_zone, "#BC2D_%s", mesh->zones[z].name );
        fprintf( fp, "o %s\n", bc_zone );
        
        uint size = mesh->zones[z].num_elems;
        uint beg  = mesh->zones[z].start;
        uint end  = beg + size;
        
        uint na, nb;
        
        for ( uint i = beg; i < end; i++ ) {
          
          /* create two tris from one quad. */

          na = mesh->mBFaces[i].nodes[0]+1;
          nb = mesh->mBFaces[i].nodes[1]+1;
          fprintf( fp, "l %d %d\n", na, nb );
        }
        
        
      }
      
    }
    
    
  }

  else if ( mesh->dim == 3 ) {
    
    /* We write only the boundary faces for 3D meshes.
     * Note that for such set we dont need to write all
     * mesh vertices.
     */
    
    /* Figure out how many vertices we need to write. */
    uint bc_vertices_counter = 0;
    for ( uint i = 0; i < mesh->sizes.bfaces; i++ ) {
      
      uint nd4 = mesh->mBFaces[i].nodes[3];
      if ( nd4 == UNINIT ) {
        bc_vertices_counter += 3;
      }
      else {
        bc_vertices_counter += 4;
      }
      
    }
    
    /* get list of boundary faces vertices */

    uint* bcvertices_idx = (uint*)malloc(bc_vertices_counter*sizeof(uint));
    
    uint k = 0;
    for ( uint z = 0; z < mesh->sizes.zones; z++ ) {

      if ( mesh->zones[z].is_boundary ) {
        
        if ( mesh->zones[z].homogeneous ) {
          
          ElemType type = mesh->zones[z].elem_type;
          
          uint size = mesh->zones[z].num_elems;
          uint beg  = mesh->zones[z].start;
          uint end  = beg + size;
          short nodes;
          
          if ( type == pTRI_3  ) nodes = 3;
          if ( type == pQUAD_4 ) nodes = 4;
            
          for ( uint i = beg; i < end; i++ ) {
            for ( short j = 0; j < nodes; j++ )
              bcvertices_idx[k++] = mesh->mBFaces[i].nodes[j];
          }

        }
        else {
          
          /* we are in a mixed zone. */
          
          for ( uint m = 0; m < mesh->zones[z].num_types; m++ ) {
            
            ElemType type = mesh->zones[z].type_list[m].elem_type;
            
            uint size = mesh->zones[z].type_list[m].num_elems;
            uint beg  = mesh->zones[z].type_list[m].start;
            uint end  = beg + size;
            short nodes;
            
            if ( type == pTRI_3  ) nodes = 3;
            if ( type == pQUAD_4 ) nodes = 4;
              
            for ( uint i = beg; i < end; i++ ) {
              for ( short j = 0; j < nodes; j++ )
                bcvertices_idx[k++] = mesh->mBFaces[i].nodes[j];
            }
            
          }
          
        }
        
      }
      
    }
    

    /* remove the duplicate vertices entries.
     * note that compare_uint(), unique() and
     * binary_search() are defined in utils.h */
    qsort( bcvertices_idx, bc_vertices_counter, sizeof(uint), compare_uint );
    bc_vertices_counter = unique( &bcvertices_idx[0], &bcvertices_idx[bc_vertices_counter] );
    
    /* write unique bc vertices. */
    double x, y, z;
    for ( uint i = 0; i < bc_vertices_counter; i++ ) {
      uint idx = bcvertices_idx[i];
      x = mesh->vertices[idx][0];
      y = mesh->vertices[idx][1];
      z = mesh->vertices[idx][2];
      //fprintf( fp, "v %.16g %.16g %.16g\n", x, y, z);
      fprintf( fp, "v %.5f %.5f %.5f\n", x, y, z);
    }
    
    
    /* write mesh wireframe before faces, so there is
     * no z-trashing. */
    
    fprintf( fp, "o #WIRE_wireframe\n" );
    for ( uint z = 0; z < mesh->sizes.zones; z++ ) {

      if ( mesh->zones[z].is_boundary ) {
        
        if ( mesh->zones[z].homogeneous ) {
          
          ElemType type = mesh->zones[z].elem_type;
          
          uint size = mesh->zones[z].num_elems;
          uint beg  = mesh->zones[z].start;
          uint end  = beg + size;
          
          if ( type == pTRI_3 ) {
            
            WireMesh* edges = (WireMesh*)malloc( 3*size*sizeof(WireMesh) );
            
            k = 0;
            uint na, nb, nc;
            for ( uint i = beg; i < end; i++ ) {

              na = mesh->mBFaces[i].nodes[0];
              nb = mesh->mBFaces[i].nodes[1];
              nc = mesh->mBFaces[i].nodes[2];

              edges[k].n1  = na;
              edges[k].n2  = nb;
              k++;
              edges[k].n1  = nb;
              edges[k].n2  = nc;
              k++;
              edges[k].n1  = nc;
              edges[k].n2  = na;
              k++;
              
            }
            
            uint num_edges = k;

            num_edges = hash_3D_wireframe(edges, size, num_edges);

            for ( uint i = 0; i < num_edges; i++ ) {
              na = edges[i].n1;
              nb = edges[i].n2;
              /* find the corresponding renumbered vertices. Since the array 
               * is sorted, we do a bindary search for efficiency. */
              na = binary_search( bcvertices_idx, bc_vertices_counter, na ) +1;
              nb = binary_search( bcvertices_idx, bc_vertices_counter, nb ) +1;
              fprintf( fp, "w %d %d\n", na, nb );
            }

            free( edges );
            
          }
          else if ( type == pQUAD_4 ) {
            
            WireMesh* edges = (WireMesh*)malloc( 4*size*sizeof(WireMesh) );
            
            k = 0;
            uint na, nb, nc, nd;
            for ( uint i = beg; i < end; i++ ) {

              na = mesh->mBFaces[i].nodes[0];
              nb = mesh->mBFaces[i].nodes[1];
              nc = mesh->mBFaces[i].nodes[2];
              nd = mesh->mBFaces[i].nodes[3];

              edges[k].n1  = na;
              edges[k].n2  = nb;
              k++;
              edges[k].n1  = nb;
              edges[k].n2  = nc;
              k++;
              edges[k].n1  = nc;
              edges[k].n2  = nd;
              k++;
              edges[k].n1  = nd;
              edges[k].n2  = na;
              k++;
              
            }
            
            uint num_edges = k;
            
            num_edges = hash_3D_wireframe(edges, size, num_edges);

            for ( uint i = 0; i < num_edges; i++ ) {
              na = edges[i].n1;
              nb = edges[i].n2;
              /* find the corresponding renumbered vertices. Since the array 
               * is sorted, we do a bindary search for efficiency. */
              na = binary_search( bcvertices_idx, bc_vertices_counter, na ) +1;
              nb = binary_search( bcvertices_idx, bc_vertices_counter, nb ) +1;
              fprintf( fp, "w %d %d\n", na, nb );
            }
            
            free( edges );
            
          }
        
        }
        else {
          
          /* we are in a mixed zone. */
          
          for ( uint m = 0; m < mesh->zones[z].num_types; m++ ) {

            ElemType type = mesh->zones[z].type_list[m].elem_type;

            uint size = mesh->zones[z].type_list[m].num_elems;
            uint beg  = mesh->zones[z].type_list[m].start;
            uint end  = beg + size;

            if ( type == pTRI_3 ) {
              
              WireMesh* edges = (WireMesh*)malloc( 3*size*sizeof(WireMesh) );
              
              k = 0;
              uint na, nb, nc;
              for ( uint i = beg; i < end; i++ ) {

                na = mesh->mBFaces[i].nodes[0];
                nb = mesh->mBFaces[i].nodes[1];
                nc = mesh->mBFaces[i].nodes[2];

                edges[k].n1  = na;
                edges[k].n2  = nb;
                k++;
                edges[k].n1  = nb;
                edges[k].n2  = nc;
                k++;
                edges[k].n1  = nc;
                edges[k].n2  = na;
                k++;
                
              }
              
              uint num_edges = k;

              num_edges = hash_3D_wireframe(edges, size, num_edges);

              for ( uint i = 0; i < num_edges; i++ ) {
                na = edges[i].n1;
                nb = edges[i].n2;
                /* find the corresponding renumbered vertices. Since the array 
                * is sorted, we do a bindary search for efficiency. */
                na = binary_search( bcvertices_idx, bc_vertices_counter, na ) +1;
                nb = binary_search( bcvertices_idx, bc_vertices_counter, nb ) +1;
                fprintf( fp, "w %d %d\n", na, nb );
              }

              free( edges );
              
            }
            else if ( type == pQUAD_4 ) {
              
              WireMesh* edges = (WireMesh*)malloc( 4*size*sizeof(WireMesh) );
              
              k = 0;
              uint na, nb, nc, nd;
              for ( uint i = beg; i < end; i++ ) {

                na = mesh->mBFaces[i].nodes[0];
                nb = mesh->mBFaces[i].nodes[1];
                nc = mesh->mBFaces[i].nodes[2];
                nd = mesh->mBFaces[i].nodes[3];

                edges[k].n1  = na;
                edges[k].n2  = nb;
                k++;
                edges[k].n1  = nb;
                edges[k].n2  = nc;
                k++;
                edges[k].n1  = nc;
                edges[k].n2  = nd;
                k++;
                edges[k].n1  = nd;
                edges[k].n2  = na;
                k++;
                
              }
              
              uint num_edges = k;

              num_edges = hash_3D_wireframe(edges, size, num_edges);

              for ( uint i = 0; i < num_edges; i++ ) {
                na = edges[i].n1;
                nb = edges[i].n2;
                /* find the corresponding renumbered vertices. Since the array 
                 * is sorted, we do a bindary search for efficiency. */
                na = binary_search( bcvertices_idx, bc_vertices_counter, na ) +1;
                nb = binary_search( bcvertices_idx, bc_vertices_counter, nb ) +1;
                fprintf( fp, "w %d %d\n", na, nb );
              }
              
              free( edges );
              
            }
            
          }
          
        }
        
        
      }
      
    }
    
    
    /* write bc faces connectivity. */
    
    for ( uint z = 0; z < mesh->sizes.zones; z++ ) {
      
      if ( mesh->zones[z].is_boundary ) {
        
        /* write zone name. */
        char bc_zone[STRING_MAX_LENGHT];
        sprintf( bc_zone, "#BC3D_%s", mesh->zones[z].name );
        fprintf( fp, "o %s\n", bc_zone );
        
        if ( mesh->zones[z].homogeneous ) {
          
          ElemType type = mesh->zones[z].elem_type;
          
          uint size = mesh->zones[z].num_elems;
          uint beg  = mesh->zones[z].start;
          uint end  = beg + size;
          
          if ( type == pTRI_3 ) {
            
            uint na, nb, nc;
            for ( uint i = beg; i < end; i++ ) {
              na = mesh->mBFaces[i].nodes[0];
              nb = mesh->mBFaces[i].nodes[1];
              nc = mesh->mBFaces[i].nodes[2];
              na = binary_search( bcvertices_idx, bc_vertices_counter, na ) +1;
              nb = binary_search( bcvertices_idx, bc_vertices_counter, nb ) +1;
              nc = binary_search( bcvertices_idx, bc_vertices_counter, nc ) +1;
              fprintf( fp, "f %d %d %d\n", na, nb, nc );
            }
            
          }
        
        }
        
      }
      
    }
    
    
    /* free local buffers */
    free( bcvertices_idx );
    
  }
  
  fclose(fp);
  
}

/* This function creates unique edges from BC zones
 * elements connectivity, in order to produce a 2D
 * wireframe for GUI processing 3D meshes.
 * The data is sorted in place in the edges[] variable
 * and the function returns the number of unique edges for
 * the current BC zone. Note that the caller is responsible
 * for memory management of edges[].
 */
uint hash_3D_wireframe(WireMesh* edges, uint elems, uint size)
{
  
  short ipr = floor(log10 (size) ) + 1;
  uint hashspace = utils_next_prime( elems );
  
  uint* iprox = (uint*)malloc( size*2*sizeof(uint) );
  uint* hash  = (uint*)malloc( size*2*sizeof(uint) );
  for ( uint i = 0; i < size; i++ ) {
     hash[i] = UNINIT;
    iprox[i] = UNINIT;
  }

  uint aux;
  uint ipont;
  uint key;
  
  uint nedges = 0;
  
  for( uint i = 0; i < size; i++ ) {

    uint an = edges[i].n1;
    uint bn = edges[i].n2;

    key = get_hash_key( 2, ipr, hashspace, &an, &bn, 0, 0 );

    if ( key > hashspace ) {
      printf("\nError: increase hashspace value.\n");
      exit(1);
    }

    if ( hash[key] == UNINIT )
    {
      edges[ nedges ].n1 = an;
      edges[ nedges ].n2 = bn;
      
      /* set up hash variables. */
      iprox[nedges] = UNINIT;
      hash[key] = nedges;
      nedges++;
      
    }
    else {

      /* face might already exist. */

      uint a = edges[ hash[key] ].n1;
      uint b = edges[ hash[key] ].n2;

      key = get_hash_key( 2, ipr, hashspace, &a, &b, 0, 0 );

      if ( ! ( a == an && b == bn ) ) {
        /* we got a collision. */
        ipont = iprox[ hash[key] ];
        aux = hash[key];

label_30:
        if ( ipont == UNINIT ) {

          edges[ nedges ].n1 = an;
          edges[ nedges ].n2 = bn;

          /* set up hash variables. */
          iprox[nedges] = UNINIT;

          /* update pointer to the before-last position
           * on the collision table. */
          iprox[aux] = nedges;

          nedges++;
        }
        else {
          uint a = edges[ ipont ].n1;
          uint b = edges[ ipont ].n2;

          /* just sort the nodes */
          key = get_hash_key( 2, ipr, hashspace, &a, &b, 0, 0 );

          if ( ! ( a == an && b == bn ) ) {

            /* search for the corresponding side on the collision table. */
            aux = ipont;
            ipont = iprox[ipont];

            goto label_30;
          }
        }
        
      }
    }

  }
  
  free( hash );
  free( iprox );
  
  return nedges;
  
}
