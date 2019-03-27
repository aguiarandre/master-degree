#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "solution.h"

void solution_read(const char* filename, MeshStruct* mesh, SolutionStruct* sol)
{

  printf("\nReading solution file.\n");
  fflush(stdout);

  // Open the HDF5 file in read-only mode, using default properties.
  hid_t file = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	

	/* Check if we are reading the mesh we think we are... */
	
	/* Get solution mesh filename. */
	{
		uint sizes[2];
		sizes[0] = 1;
		sizes[1] = 1;
		
		char** data = (char**)malloc( sizes[0]*sizeof(char*) );

		hdf5_read_vlstring(file, "mesh_filename", sizes, data);

		sol->mesh_filename = (char*)malloc(STRING_MAX_LENGHT*sizeof(char));
		sol->mesh_filename = data[0];

		free( data );

		if ( (strcmp(sol->mesh_filename,mesh->filename) != 0) )
		{
			printf("\nError: Current solution was not created for the \"%s\" mesh file.", mesh->filename);
			exit(1);
		}
	}


	/* Get mesh sizes for cross-referencing/check. */
	{
		typedef struct
		{
			uint nodes;
			uint elems;
			uint faces;
			uint bfaces;
			uint zones;
			uint variables;
		} MetaSize;
		
		MetaSize* data = (MetaSize*)malloc( sizeof(MetaSize) );
	
		size_t dst_size = sizeof( MetaSize );
		
		size_t dst_offset[6] = { HOFFSET( MetaSize, nodes ),
								 HOFFSET( MetaSize, elems ),
								 HOFFSET( MetaSize, faces ),
								 HOFFSET( MetaSize, bfaces ),
								 HOFFSET( MetaSize, zones ),
								 HOFFSET( MetaSize, variables )};
		
		size_t dst_sizes[6] = { sizeof( data[0].nodes),
								sizeof( data[0].elems),
								sizeof( data[0].faces),
								sizeof( data[0].bfaces),
								sizeof( data[0].zones),
								sizeof( data[0].variables)};
		
		H5TBread_table( file, "meta_sizes", dst_size, dst_offset, dst_sizes, data );
		
		int error = 0;
		if ( mesh->sizes.nodes != data[0].nodes ) {
			error = 1;
		} else if ( mesh->sizes.elems != data[0].elems ) {
			error = 1;
		} else if ( mesh->sizes.faces != data[0].faces ) {
			error = 1;
		} else if ( mesh->sizes.bfaces != data[0].bfaces ) {
			error = 1;
		} else if ( mesh->sizes.zones  != data[0].zones ) {
			error = 1;
		}
		
		if ( error > 0 ) {
			printf("\nError: Solution data size does not match the one in \"%s\" mesh file.", mesh->filename);
			exit(1);
		}
		
		sol->num_vars = data[0].variables;
		
		free( data );
	}
	
	
	/* Get variables meta information. */
	{
		typedef struct
		{
			int idx;
			char name[STRING_MAX_LENGHT];
			char var_type[STRING_MAX_LENGHT]; // scalar,vector,etc
			char num_type[STRING_MAX_LENGHT]; // int,double,etc
			int dim; // number of components
			int size;
			char location[STRING_MAX_LENGHT];
		} MetaVariable;
		
		MetaVariable* data = (MetaVariable*)malloc( sol->num_vars*sizeof(MetaVariable) );
	
		size_t dst_size = sizeof( MetaVariable );
		
		size_t dst_offset[7] = { HOFFSET( MetaVariable, idx ),
								 HOFFSET( MetaVariable, name ),
								 HOFFSET( MetaVariable, var_type ),
								 HOFFSET( MetaVariable, num_type ),
								 HOFFSET( MetaVariable, dim ),
								 HOFFSET( MetaVariable, size ),
								 HOFFSET( MetaVariable, location )};
		
		size_t dst_sizes[7] = { sizeof( data[0].idx),
								sizeof( data[0].name),
								sizeof( data[0].var_type),
								sizeof( data[0].num_type),
								sizeof( data[0].dim),
								sizeof( data[0].size),
								sizeof( data[0].location)};
		
		H5TBread_table( file, "meta_var", dst_size, dst_offset, dst_sizes, data );
		
		
		/* allocate and populate the variables in the solution struct. */
		
		sol->vars = (Variable*)malloc( sol->num_vars*sizeof(Variable) );
		
		for ( int i = 0; i < sol->num_vars; i++ )
		{
			sol->vars[i].idx = data[i].idx;
			
			sol->vars[i].name = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
			strcpy(sol->vars[i].name,data[i].name);
			
			sol->vars[i].var_type = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
			strcpy(sol->vars[i].var_type,data[i].var_type);
			
			sol->vars[i].num_type = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
			strcpy(sol->vars[i].num_type,data[i].num_type);
			
			sol->vars[i].dim = data[i].dim;
			
			sol->vars[i].size = data[i].size;
			
			if ( strcmp(data[i].location,"nodes") == 0 ) {
				sol->vars[i].location = 0;
			} else if ( strcmp(data[i].location,"elems") == 0 ) {
				sol->vars[i].location = 1;
			} else if ( strcmp(data[i].location,"faces") == 0 ) {
				sol->vars[i].location = 2;
			} else if ( strcmp(data[i].location,"bfaces") == 0 ) {
				sol->vars[i].location = 3;
			}
		}
		free( data );
	}
	
	
	/* Get total number of variables (scalar fields) */
	int num_total_vars = 0;
	for( int i = 0; i < sol->num_vars; i++ )
	{
		for ( int n = 0; n < sol->vars[i].dim; n++ ) num_total_vars++;
	}
	
	/* Get variable component's name */
	
	{
		typedef struct
		{
			int idx;
			int component;
			char name[STRING_MAX_LENGHT];
		} MetaVarNames;
		
		MetaVarNames* data = (MetaVarNames*)malloc( num_total_vars*sizeof(MetaVarNames) );
	
		size_t dst_size = sizeof( MetaVarNames );
		
		size_t dst_offset[3] = { HOFFSET( MetaVarNames, idx ),
								 HOFFSET( MetaVarNames, component ),
								 HOFFSET( MetaVarNames, name )};
		
		size_t dst_sizes[3] = { sizeof( data[0].idx),
								sizeof( data[0].component),
								sizeof( data[0].name)};
		
		H5TBread_table( file, "meta_var_names", dst_size, dst_offset, dst_sizes, data );
		
		for( int i = 0; i < sol->num_vars; i++ )
		{
			int flag = 0;
			for ( int j = 0; j < num_total_vars; j++)
			{
				if ( sol->vars[i].idx == data[j].idx )
				{
					if ( flag == 0 )
					{
						sol->vars[i].componentName = (char**)malloc( sol->vars[i].dim*sizeof(char*) );
						flag = 1;
					}
					
					for ( int n = 0; n < sol->vars[i].dim; n++ )
					{
						if ( n == data[j].component )
						{
							sol->vars[i].componentName[n] = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
							strcpy(sol->vars[i].componentName[n],data[j].name);
						}
					}
				}
			}
		}
		
		free(data);
	}
	
	/* Get actual data for the variable. */

	{
		for( int i = 0; i < sol->num_vars; i++ )
		{
			if ( (strcmp(sol->vars[i].num_type,"double") == 0) )
			{
				int size = sol->vars[i].dim*sol->vars[i].size;
				double* data = (double*)malloc( size*sizeof(double) );
				
				char var_name_id[STRING_MAX_LENGHT];
				strcpy(var_name_id,"data_");
				strcat(var_name_id,sol->vars[i].name);

				hdf5_read_data(file, var_name_id, "double", data);

				sol->vars[i].data = (char*)malloc( size*sizeof(double) );
				memcpy( sol->vars[i].data, data, size*sizeof(double) );
				free(data);
			}
		}
	}
	
	H5Fclose (file);
}


void solution2nodes(MeshStruct *mesh, SolutionStruct* sol, SolutionStruct* nodalSol)
{
	/* The variables that are not located at mesh nodes are
	   interpolated to them. The nodal averaging method presented
	   by "Neal Tilson Frink, Three-Dimesional Upwind Scheme for
	    Solving the Euler Equations on Unstructured Tetrahedral
		Grids, sep. 1994 (PhD Thesis)" is used here.
	*/
	
	int num_total_vars = 0; // total number of variables (scalars + vector components)
	
	for( int i = 0; i < sol->num_vars; i++ )
	{
		for ( int n = 0; n < sol->vars[i].dim; n++ ) num_total_vars++;
	}
	
	
	/* Get list of unique variables that are not located at mesh nodes. */
	
	for( int i = 0; i < sol->num_vars; i++ )
	{
		if ( sol->vars[i].location != 0 )
		{
			for( int j = i+1; j < sol->num_vars; j++ )
			{
				if ( j != i )
				{
					if ( sol->vars[i].dim == sol->vars[j].dim )
					{
						for ( int n = 0; n < sol->vars[i].dim; n++ )
						{
							if ( strcmp(sol->vars[i].componentName[n],sol->vars[j].componentName[n]) == 0 )
							{
								num_total_vars--;
							}
						}
					}
				}
			}
		}
	}
	
	nodalSol->num_vars = sol->num_vars;
	
	nodalSol->vars = (Variable*)malloc( nodalSol->num_vars*sizeof(Variable) );
	
	
	/* compute mesh elements centroid */
	
	double* xc = (double*)calloc( mesh->sizes.elems, sizeof(double) );
	double* yc = (double*)calloc( mesh->sizes.elems, sizeof(double) );
	
	for( int i = 0; i < mesh->sizes.elems; i++ )
	{
		int num_nodes = 0;
		for( int n = 0; n < 4; n++ )
		{
			if ( mesh->mElems[i].nodes[n] == UNINIT ) break;
			xc[i] += mesh->vertices[ mesh->mElems[i].nodes[n] ][0];
			yc[i] += mesh->vertices[ mesh->mElems[i].nodes[n] ][1];
			num_nodes++;
		}
		xc[i] /= (double)num_nodes;
		yc[i] /= (double)num_nodes;
	}
	
	
	/* compute boundary faces centroid */
	double* xc_bf = (double*)calloc( mesh->sizes.bfaces, sizeof(double) );
	double* yc_bf = (double*)calloc( mesh->sizes.bfaces, sizeof(double) );
	
	for( int i = 0; i < mesh->sizes.bfaces; i++ )
	{
		for( int n = 0; n < 2; n++ )
		{
			xc_bf[i] += mesh->vertices[ mesh->mBFaces[i].nodes[n] ][0];
			yc_bf[i] += mesh->vertices[ mesh->mBFaces[i].nodes[n] ][1];
		}
		xc_bf[i] /= 2.0f;
		yc_bf[i] /= 2.0f;
	}
	

	int repeated_var = 0;
	
	for( int i = 0; i < sol->num_vars; i++ )
	{
		/* setup the nodalSol variable struct. */
		
		nodalSol->vars[i].idx = sol->vars[i].idx;
		nodalSol->vars[i].dim = sol->vars[i].dim;
		nodalSol->vars[i].size = mesh->sizes.nodes;
		nodalSol->vars[i].location = 0;
		
		nodalSol->vars[i].name = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
		strcpy(nodalSol->vars[i].name,sol->vars[i].name);
		
		nodalSol->vars[i].componentName = (char**)malloc( nodalSol->vars[i].dim*sizeof(char*) );
		for( int n = 0; n < nodalSol->vars[i].dim; n++ )
		{
			nodalSol->vars[i].componentName[n] = (char*)malloc( STRING_MAX_LENGHT*sizeof(char) );
			strcpy(nodalSol->vars[i].componentName[n],sol->vars[i].componentName[n]);
		}
	

		if ( sol->vars[i].location == 1 )
		{
			/* variable is stored at element */

			double* num  = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );
			double* den  = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );		
			double* data = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );
		
			for( int n = 0; n < sol->vars[i].dim; n++ )
			{
				for( int e = 0; e < mesh->sizes.elems; e++ )
				{
					int num_nodes = 4;
					int vid = (e*sol->vars[i].dim + n)*sizeof(double);
					double* q = (double*)(sol->vars[i].data+vid);
					for ( int j = 0; j < num_nodes; j++ )
					{
						/* compute solution at nodes */
						int idx = mesh->mElems[e].nodes[j];
						double xnd = mesh->vertices[ idx ][0];
						double ynd = mesh->vertices[ idx ][1];
						
						double ri = sqrt( (xc[e]-xnd)*(xc[e]-xnd) + (yc[e]-ynd)*(yc[e]-ynd) );
						num[idx*sol->vars[i].dim+n] += (*q)/ri;
						den[idx*sol->vars[i].dim+n] += 1.f/ri;
					}
				}
				
				for( int j = 0; j < mesh->sizes.nodes; j++ )
				{
					data[j*sol->vars[i].dim+n]  = num[j*sol->vars[i].dim+n]/den[j*sol->vars[i].dim+n];
				}
			}
			
			int size = mesh->sizes.nodes*sol->vars[i].dim*sizeof(double);
			nodalSol->vars[i].data = (char*)malloc( size );
			memcpy( nodalSol->vars[i].data, data, size );

			free(num);
			free(den);
			free(data);
		}
		else if ( sol->vars[i].location == 2 )
		{
			/* variable is stored at faces. */
			
		}
		else if ( sol->vars[i].location == 3 )
		{
			/* variable is stored at bfaces. Note that bface variables
			   take precedence over other interpolations. The idea is
			   that bface values are computed by boundary conditions and
			   should not be interpolated.
			 */
			
			double* num  = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );
			double* den  = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );		
			double* data = (double*)calloc( mesh->sizes.nodes*sol->vars[i].dim, sizeof(double) );
			
			for( int n = 0; n < sol->vars[i].dim; n++ )
			{
				for( int f = 0; f < mesh->sizes.bfaces; f++ )
				{
					int num_nodes = 2;
					int vid = (f*sol->vars[i].dim + n)*sizeof(double);
					double* q = (double*)(sol->vars[i].data+vid);
					for ( int j = 0; j < num_nodes; j++ )
					{
						/* compute solution at nodes */
						int idx = mesh->mBFaces[f].nodes[j];
						double xnd = mesh->vertices[ idx ][0];
						double ynd = mesh->vertices[ idx ][1];
						
						double ri = sqrt( (xc_bf[f]-xnd)*(xc_bf[f]-xnd) + (yc_bf[f]-ynd)*(yc_bf[f]-ynd) );
						num[idx*sol->vars[i].dim+n] += (*q)/ri;
						den[idx*sol->vars[i].dim+n] += 1.f/ri;
					}
				}

				/* check if this variable already exists */
				for ( int j = 0; j < i; j++ )
				{
					if ( strcmp(sol->vars[i].componentName[n],sol->vars[j].componentName[n]) == 0 )
					{
						data = (double*)(nodalSol->vars[j].data);
						repeated_var++;
						break;
					}
				}
				
				for( int f = 0; f < mesh->sizes.bfaces; f++ )
				{
					int num_nodes = 2;
					for ( int j = 0; j < num_nodes; j++ )
					{
						int idx = mesh->mBFaces[f].nodes[j];
						data[idx*sol->vars[i].dim+n] = num[idx*sol->vars[i].dim+n]/den[idx*sol->vars[i].dim+n];
					}
				}
			}
			
			int size = mesh->sizes.nodes*sol->vars[i].dim*sizeof(double);
			nodalSol->vars[i].data = (char*)malloc( size );
			memcpy( nodalSol->vars[i].data, data, size );

			free(num);
			free(den);
			free(data);
		}
	}
	
	if ( repeated_var == sol->vars[0].dim ) nodalSol->num_vars--;
	
	free(xc); free(yc);
	free(xc_bf); free(yc_bf);
}
