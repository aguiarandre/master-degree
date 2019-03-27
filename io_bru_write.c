#include <stdio.h>
#include <string.h> /* memcpy() */
#include <stdlib.h>
#include <math.h>

#include "io_bru.h"

void io_bru_write(MeshStruct* mesh)
{
	FILE* fp;
	
	printf("Writing fort.2 ...\n");
	fflush(stdout);
	
	// Fort.2

	fp = fopen( "fort.2", "w" );
	
	fprintf(fp,"%d\n",mesh->sizes.nodes);
	
	for(uint i = 0; i < mesh->sizes.nodes; i++ )
	{
		double x = mesh->vertices[i][0];
		double y = mesh->vertices[i][1];
		double z = mesh->vertices[i][2];
		fprintf(fp,"%f %f %f\n",x,y,z);
	}
	
	fclose(fp);
	
	
	// Fort.3
	
	printf("Writing fort.3 ...\n");
	fflush(stdout);
	
	fp = fopen( "fort.3", "w" );
	
	fprintf(fp, "%d\n",mesh->sizes.elems);

		
	for(uint z = 0; z < mesh->sizes.zones; z++ )
	{
		if ( ! mesh->zones[z].is_boundary )
		{
			if ( mesh->zones[z].homogeneous )
			{
				if ( mesh->zones[z].elem_type == pTETRA_4 )
				{
					uint beg = mesh->zones[z].start;
					uint end = beg + mesh->zones[z].num_elems;
					for ( uint i = beg; i < end; i++ )
					{
						fprintf(fp, "4");
						for ( short j = 0; j < 4; j++ ) {
							uint node = mesh->mElems[i].nodes[j]+1;
							fprintf(fp, " %d", node);
						}
						fprintf(fp,"\n");
					}
				}
				else if ( mesh->zones[z].elem_type == pHEXA_8 )
				{
					uint beg = mesh->zones[z].start;
					uint end = beg + mesh->zones[z].num_elems;
					for ( uint i = beg; i < end; i++ )
					{
						fprintf(fp, "8");
						for ( short j = 0; j < 8; j++ ) {
							uint node = mesh->mElems[i].nodes[j]+1;
							fprintf(fp, " %d", node);
						}
						fprintf(fp,"\n");
					}
				}
				else
				{
					printf("ERROR: element type not supported.\n");
					exit(1);
				}
			}
			else
			{
				for ( short k = 0; k < mesh->zones[z].num_types; k++ )
				{
					if ( mesh->zones[z].type_list[k].elem_type == pTETRA_4 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
							fprintf(fp, "4");
							for ( short j = 0; j < 4; j++ ) {
								uint node = mesh->mElems[i].nodes[j]+1;
								fprintf(fp, " %d", node);
							}
							fprintf(fp,"\n");
						}
					}
					else if ( mesh->zones[z].type_list[k].elem_type == pPYRA_5 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
							fprintf(fp, "5");
							for ( short j = 0; j < 5; j++ ) {
								uint node = mesh->mElems[i].nodes[j]+1;
								fprintf(fp, " %d", node);
							}
							fprintf(fp,"\n");
						}
					}
					else if ( mesh->zones[z].type_list[k].elem_type == pPENTA_6 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
							fprintf(fp, "6");
							for ( short j = 0; j < 6; j++ ) {
								uint node = mesh->mElems[i].nodes[j]+1;
								fprintf(fp, " %d", node);
							}
							fprintf(fp,"\n");
						}
					}
					else if ( mesh->zones[z].type_list[k].elem_type == pHEXA_8 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
							fprintf(fp, "8");
							for ( short j = 0; j < 8; j++ ) {
								uint node = mesh->mElems[i].nodes[j]+1;
								fprintf(fp, " %d", node);
							}
							fprintf(fp,"\n");
						}
					}
				else
				{
					printf("ERROR: element type not supported.\n");
					exit(1);
				}
			}
			}

		}
	}
	
	fclose(fp);
	
	
	// Fort.5
	
	printf("Writing fort.5 ...\n");
	
	printf("\n");
	printf("Select the appropriate zone type\n");
	printf("1 - wall\n");
	printf("2 - farfield\n");
	printf("3 - wall\n");
	printf("4 - farfield\n");
	printf("5 - wall\n");
	printf("6 - farfield\n");
	printf("7 - wall\n");
	printf("8 - farfield\n");
	printf("8 - farfield\n");
	printf("8 - farfield\n");
	printf("8 - farfield\n");
	printf("8 - farfield\n");
	printf("8 - farfield\n");
	printf("\n");
	fflush(stdout);
	
	fp = fopen( "fort.5", "w" );
	
	fprintf(fp, "%d\n",mesh->sizes.bfaces);
	
	for(uint z = 0; z < mesh->sizes.zones; z++)
	{
		if ( mesh->zones[z].is_boundary )
		{
			char* bc_name = mesh->zones[z].name;

			// Query user about boundary types.

		char buffer[256];

			printf("zone: %s, type: ",bc_name);
			fflush(stdout);
			fgets(buffer,256,stdin);
			mesh->zones[z].bc_id = atoi(buffer);
			
			
			if ( mesh->zones[z].homogeneous )
			{
				if ( mesh->zones[z].elem_type == pTRI_3 )
				{
					uint beg = mesh->zones[z].start;
					uint end = beg + mesh->zones[z].num_elems;
					for ( uint i = beg; i < end; i++ )
					{
						fprintf(fp, "%d",mesh->zones[z].bc_id);
						for ( short j = 0; j < 3; j++ ) {
							uint node = mesh->mBFaces[i].nodes[j]+1;
							fprintf(fp," %d",node);
						}
						fprintf(fp,"\n");
					}
				}
				else if ( mesh->zones[z].elem_type == pQUAD_4 )
				{
					uint beg = mesh->zones[z].start;
					uint end = beg + mesh->zones[z].num_elems;
					for ( uint i = beg; i < end; i++ )
					{
						fprintf(fp, "%d",mesh->zones[z].bc_id);
						for ( short j = 0; j < 4; j++ ) {
							uint node = mesh->mBFaces[i].nodes[j]+1;
							fprintf(fp," %d",node);
						}
						fprintf(fp,"\n");
					}
				}
				else
				{
					printf("ERROR: element type not supported.\n");
					exit(1);
				}
			}
			else
			{
				for ( short k = 0; k < mesh->zones[z].num_types; k++ )
				{
					if ( mesh->zones[z].type_list[k].elem_type == pTRI_3 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
						fprintf(fp, "%d",mesh->zones[z].bc_id);
						for ( short j = 0; j < 3; j++ ) {
							uint node = mesh->mBFaces[i].nodes[j]+1;
							fprintf(fp," %d",node);
						}
						fprintf(fp,"\n");
						}
					}
					else if ( mesh->zones[z].type_list[k].elem_type == pQUAD_4 )
					{
						uint beg = mesh->zones[z].type_list[k].start;
						uint end = beg + mesh->zones[z].type_list[k].num_elems;
						for ( uint i = beg; i < end; i++ )
						{
						fprintf(fp, "%d",mesh->zones[z].bc_id);
						for ( short j = 0; j < 4; j++ ) {
							uint node = mesh->mBFaces[i].nodes[j]+1;
							fprintf(fp," %d",node);
						}
						fprintf(fp,"\n");
						}
					}
				else
				{
					printf("ERROR: element type not supported.\n");
					exit(1);
				}
			}
			}

		}
	}
	
	fclose(fp);
	
	printf("Finished ...\n\n");
	fflush(stdout);
	
}

