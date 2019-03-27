#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "io_cgns.h"

/* This routine */

void io_cgns_write(const char* filename, MeshStruct* mesh, SolutionStruct* sol)
{
  printf("\nWriting CGNS file.\n");
  fflush(stdout);
  
  int cgfile;
  int baseid, zoneid, sectionid;
  char basename[STRING_MAX_LENGHT];
  char zonename[STRING_MAX_LENGHT];
  
  cg_open(filename, CG_MODE_WRITE, &cgfile);
  
  strcpy(basename, "Base");
  int icelldim=mesh->dim;
  int iphysdim=mesh->dim;
  
  cg_base_write(cgfile, basename, icelldim, iphysdim, &baseid);
  
  
	/* create zone */
	
	cgsize_t sizes[3];
	sizes[0] = mesh->sizes.nodes;
	sizes[1] = mesh->sizes.elems;
	sizes[2] = mesh->sizes.bfaces;
	
	strcpy(zonename, "Zone");
	
    cg_zone_write(cgfile, baseid, zonename, sizes, Unstructured, &zoneid);
	
	
	/* write mesh coordinates */
	{
		int coordid;
		double* data;	
		data = (double*)malloc( mesh->sizes.nodes*sizeof(double) );
		
		for( uint i = 0; i < mesh->sizes.nodes; i++ )
		{
			data[i] = mesh->vertices[i][0];
		}
		cg_coord_write(cgfile, baseid, zoneid, RealDouble, "CoordinateX", data, &coordid);
		
		for( uint i = 0; i < mesh->sizes.nodes; i++ )
		{
			data[i] = mesh->vertices[i][1];
		}
		cg_coord_write(cgfile, baseid, zoneid, RealDouble, "CoordinateY", data, &coordid);

		if ( mesh->dim == 3 )
		{
			for( uint i = 0; i < mesh->sizes.nodes; i++ )
			{
				data[i] = mesh->vertices[i][2];
			}
			cg_coord_write(cgfile, baseid, zoneid, RealDouble, "CoordinateZ", data, &coordid);
		}
		free(data);
	}
	
	
	/* write element connectivity. Note that CGNS
	   connectivity indices have base 1.
	*/
	{
		uint elem_count;
		ElementType_t etype;
		uint sizes[2];
		
		if ( mesh->dim == 2 ) {
		  
		  char zonename[STRING_MAX_LENGHT];
		  
		  
		  for ( int z = 0; z < mesh->sizes.zones; z++ )
		  {

			if ( mesh->zones[z].is_boundary == false )
			{
				// write volumetric connectivity.
				
				strcpy(zonename, mesh->zones[z].name);
				
			  // write mesh quads
			  etype = QUAD_4;
			  elem_count = mesh->sizes.elem_type_count[pQUAD_4];
			  int k = 0;

			  if ( elem_count > 0 ) {
				sizes[0] = elem_count;
				sizes[1] = 4;
				cgsize_t* data = (cgsize_t*)malloc(sizes[0]*sizes[1]*sizeof(cgsize_t));

				for( uint i = 0; i < sizes[0]; i++ ) {
				  for ( short j = 0; j < sizes[1]; j++ ) {
					data[k++] = (cgsize_t)(mesh->mElems[i].nodes[j] + 1);
				  }
				}
				cg_section_write(cgfile, baseid, zoneid, zonename, etype, 0,
					elem_count-1, 0, data, &sectionid);

				free(data);
			  }
			} else {
				// write boundary connectivity.

				strcpy(zonename, mesh->zones[z].name);
				
				etype = BAR_2;
				elem_count = mesh->zones[z].num_elems;
				int k = 0;

				sizes[0] = elem_count;
				sizes[1] = 2;
				cgsize_t* data = (cgsize_t*)malloc( sizes[0]*sizes[1]*sizeof(cgsize_t) );
				
				uint beg = mesh->zones[z].start;
				uint end = beg + mesh->zones[z].num_elems;

				for( uint i = beg; i < end; i++ ) {
				  for ( short j = 0; j < sizes[1]; j++ ) {
					data[k++] = (cgsize_t)(mesh->mBFaces[i].nodes[j] + 1);
				  }
				}
				cg_section_write(cgfile, baseid, zoneid, zonename, etype, 0,
					mesh->zones[z].num_elems-1, 0, data, &sectionid);

				free(data);
			}
		  }
		  
		  
		}
	}

	/* interpolate the solution to mesh nodes, if necessary. */
	SolutionStruct nodalSol;
	solution2nodes(mesh,sol,&nodalSol);


	// write solution
	{
		int flowid, fieldid;
		
		// define flow solution node name
		char solname[33];
		strcpy(solname,"FlowSolution");

		// create flow solution node
		cg_sol_write(cgfile, baseid, zoneid, solname, Vertex, &flowid);
		
		for( int i = 0; i < nodalSol.num_vars; i++ )
		{
			if ( nodalSol.vars[i].dim > 1 )
			{
				// we have to write a component at a time.
				for ( int j = 0; j < nodalSol.vars[i].dim; j++ )
				{
					double* data = (double*)malloc( nodalSol.vars[i].size*sizeof(double) );
					
					for( int k = 0; k < nodalSol.vars[i].size; k++ )
					{
						int idx = (k*nodalSol.vars[i].dim + j)*sizeof(double);
						data[k] = *(double*)(nodalSol.vars[i].data+idx);
					}

					cg_field_write(cgfile, baseid, zoneid, flowid, RealDouble,
						nodalSol.vars[i].componentName[j], data, &fieldid);

					free(data);
				}
			}
		}
	}

	
	/* close file. */
    cg_close(cgfile);
}
