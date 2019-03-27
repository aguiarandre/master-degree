#include "native.h"
#include <stdlib.h> // malloc
#include <string.h>
#include <stddef.h>
#include <math.h>
#include "rbf_functions.h"
#include "moveMesh.h"

         
#define ind2d(i,j) i*tam+j
void writeMeshMovement ( MeshStruct * mesh , int nMesh );

void moveMesh( MeshStruct * mesh )
{ 
	printf("\n");
	printf("Using RBF to move the mesh ...\n");
    double T = 0.0;
    double deltat = 0.2;
    
    int nMesh = 0;
    for ( double time = 0.0 ; time <= T ; time = time + deltat, nMesh++ )
    {
        mesh->time = time;
        mesh->deltat = deltat;
        
        /** ALLOCATE MEMORY */
        initializeData( mesh );
	
        /* CREATING NECESSARY DATA */
        getBoundaryNodes( mesh );
        
        
        getMovingNodes( mesh );
        
        
        /* MATRICES FROM THESE DATA */
        populateMatrices( mesh );
        
        /* CALCULATE DISPLACEMENTS */
        applyDisplacement( mesh );
        
        writeMeshMovement( mesh, nMesh );
        
    
	
        /* FREE UP DATA */
        destroyData( mesh );
    
    }        
	
	return;
}



/* MESH_ORDER needs to be defined. Can be hard coded. It is the 
polynomial order of the faces of your unstructured mesh. Default = 1 */
 

void getBoundaryNodes( MeshStruct * mesh )
{
    int startingElem, endingElem;
	int nodeNumber = 0, nodeCounter = 0;
    double theta = 80.0 * 3.14 / 180.0;
    double x, y, x_rot, y_rot;
    double t = mesh->time; 
    double deltat = mesh->deltat;
    double b1, b2, b3, b1_b, b2_b, b3_b;
    double c_rot = 3.0;
    
    b1 = t * t * ( t * t - 4.0 * t + 4.0);
    b2 = t * t * ( 3.0 - t ) / 4.0;
    b3 = t * t * t * ( - 8.0 * t * t * t + 51.0 * t * t - 111.0 * t + 84.0 ) / 16.0;
    
    b1_b = 0.0; b2_b = 0.0; b3_b = 0.0;
    
    if (t > deltat)
    {
        t = t - deltat;
        b1_b = t * t * ( t * t - 4.0 * t + 4.0);
        b2_b = t * t * ( 3.0 - t ) / 4.0;
        b3_b = t * t * t * ( - 8.0 * t * t * t + 51.0 * t * t - 111.0 * t + 84.0 ) / 16.0;
        
        
    }
    
    double dy = (b2 - b2_b);
    
    printf(" getting boundary nodes.. t = %lf, b1 = %.2lf , dy = %.2lf, b3 = %.2lf \n", t, b1 , dy, b3);
    
    //b1 = 1.0; b1_b = 0.0;
    
    //b2 = 1.0;
    //b2_b = 0.0;
    
    // Run through boundary faces and get nodeNumbers.
	for (int zoneId = 0; zoneId < mesh->sizes.zones; zoneId++ )
	{
	
    	startingElem = mesh->zones[ zoneId ].start;
		endingElem 	 = mesh->zones[ zoneId ].start + mesh->zones[ zoneId ].num_elems ;
		
		// if it is a indeed a boundary ...
		if ( mesh->zones[ zoneId ].is_boundary == true )
		{
			// run through its boundary faces,
			for (int bFaceNumber = startingElem; bFaceNumber < endingElem ; bFaceNumber++)
			{
				// collect all nodes it has.
				
				for (int nodeId = 0; nodeId < 2 + (MESH_ORDER-1); nodeId++) 
				{
					// store the number of each node.
					nodeNumber = mesh->mBFaces[ bFaceNumber ].nodes[ nodeId ] ;
					
					// Now, only enter here if you have not previously entered!
					if ( !mesh->rbf->isSet[ nodeNumber ] )
					{
						mesh->rbf->isSet[ nodeNumber ] = 1; // leave a flag saying you have been here.
						
						// Store data:
						mesh->rbf->fixedNodesBefore[ nodeCounter ].nodeNumber = nodeNumber;
						mesh->rbf->fixedNodesBefore[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0];
						mesh->rbf->fixedNodesBefore[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1];
						mesh->rbf->fixedNodesBefore[ nodeCounter ].z = 0.0;
						
                        // Make the movement at these nodes.
                        // if it is the airfoil:
                        if ( mesh->zones[ zoneId ].bc_id == 1 )
                        {
                            mesh->rbf->fixedNodes[ nodeCounter ].nodeNumber = nodeNumber;
                            x = mesh->vertices[ nodeNumber ][0];
                            y = mesh->vertices[ nodeNumber ][1];
                            
                            // apply rotation to x & y:

                            x_rot = ( x - 1.0 / 3.0 ) * cos( theta * (b1 - b1_b) ) + ( y - b2 ) * sin( theta * (b1 - b1_b) )  ;
                            y_rot = ( y - b2 ) * cos( theta * (b1 - b1_b) ) - (x - 1.0 / 3.0 ) * sin( theta * (b1 - b1_b) )  ;
                            
                    
                            mesh->rbf->fixedNodes[ nodeCounter ].x = (1.0/3.0 + x_rot) ;
                            mesh->rbf->fixedNodes[ nodeCounter ].y = b2 + y_rot + dy ;
                            
                            mesh->vertices[ nodeNumber ][0] = mesh->rbf->fixedNodes[ nodeCounter ].x ;
                            mesh->vertices[ nodeNumber ][1] = mesh->rbf->fixedNodes[ nodeCounter ].y ;
                            
                            mesh->rbf->fixedNodes[ nodeCounter ].z = 0.0;
                        }
                        else // remain fixed
                        {
                            mesh->rbf->fixedNodes[ nodeCounter ].nodeNumber = nodeNumber;
                            mesh->rbf->fixedNodes[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0] ;
                            mesh->rbf->fixedNodes[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1] ;
                            mesh->rbf->fixedNodes[ nodeCounter ].z = 0.0;
                        }
                        
                        
                        
						// Go to next node.
						nodeCounter++;
						
					} // End if node isSet 
					
				} // End nodeId loop
						
			} // End for boundary elements 
			
		} // End if boundary
		
	} // End zones
	
	// Store number of fixed nodes.
    
	mesh->rbf->nFixedNodes = nodeCounter;
	
	return;
}


void writeMeshMovement ( MeshStruct * mesh , int nMesh )
{
	int nElems = mesh->sizes.elems + mesh->sizes.bfaces;	// as per gmsh way of writing elements
	double x,y,z;
	int nodeNumber;
    
    int nNodesFace = (mesh->sizes.faces + mesh->sizes.bfaces) * (MESH_ORDER-1);
    
    int nInternalNodes = mesh->sizes.elems * (MESH_ORDER-1) * (MESH_ORDER-1);
    
    // The total # of nodes is:
    int nNodes = mesh->sizes.nodes +    // normal nodes
                 nNodesFace +           // HO nodes 
                 nInternalNodes;        // Cell nodes	
	
    
    char str[80];
    sprintf( str, "mesh_test_%d.gmsh", nMesh );
	FILE * fw = fopen( str ,"w");
   
    
    // List of nodes:
	fprintf(fw, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n", nNodes);
	
	for (int i = 0; i < nNodes; i++ )
	{
		nodeNumber = i;
		x = mesh->vertices[nodeNumber][0];
		y = mesh->vertices[nodeNumber][1];
		z = 0.0;
		
		fprintf(fw, "%d %.15lf %.15lf %.15lf\n", nodeNumber + 1, x, y, z);
    }
	
    
    // Conectivity Table:
	int elemNumber, elemType, elemBType, nTags, bcId , elementaryTag, nNodesPerCell, beg, end;
	
	fprintf(fw, "$EndNodes\n$Elements\n%d\n", nElems);
	
    
	if ( MESH_ORDER == 2 )
	{
		elemType = 10; // quad8 = 16; tri6 = 9; quad9 = 10; 3-node boundary face  = 8
		elemBType = 8;
	}
	else if ( MESH_ORDER == 3 )
	{	
		elemType = 36;
		elemBType = 26;
	}
	else if ( MESH_ORDER == 4 )
	{
		elemType = 37;
		elemBType = 27;
	}
	nTags = 2;
	elementaryTag = 0;	
	
	// 				# of nodes + # of HO nodes per face * nFaces + # of cell nodes
	nNodesPerCell = 4 + 4 * (MESH_ORDER-1) + (MESH_ORDER-1)*(MESH_ORDER-1);
	
	
	for (int zoneId = 0; zoneId < mesh->sizes.zones; zoneId++ )
	{
		if( !mesh->zones[ zoneId ].is_boundary )
		{
			bcId = mesh->zones[ zoneId ].bc_id;
			
			beg = mesh->zones[ zoneId ].start;
			end = mesh->zones[ zoneId ].start + mesh->zones[ zoneId ].num_elems; 
			
			
			// Then, loop through internal cells.
			for (elemNumber =  beg; elemNumber < end ; elemNumber++)
			{
				fprintf(fw, "%d %d %d %d %d", elemNumber+1, elemType, nTags, bcId , elementaryTag );


				// Now, loop throughout the faces of this cell
				for (int nodeId = 0; nodeId < nNodesPerCell; nodeId++ )
                 	fprintf(fw, " %d", mesh->mElems[ elemNumber ].nodes[ nodeId ] + 1 );
				// skip line
				fprintf(fw, "\n");
				
				// and that's it for internal elements (boundaries are included here, though)
				
			}
			
			// Now, loop through boundary faces:
		}
		else // if is_boundary
		{
			bcId = mesh->zones[ zoneId ].bc_id;

			beg = mesh->zones[ zoneId ].start ;
			end = mesh->zones[ zoneId ].start + mesh->zones[ zoneId ].num_elems;

			for (int bFaceNumber = beg, faceNumber; bFaceNumber < end; bFaceNumber++ )
			{
				faceNumber = elemNumber + bFaceNumber + 1;
				
				fprintf(fw, "%d %d %d %d %d", faceNumber , elemBType, nTags, bcId , elementaryTag );
			    
				// For each of these b-faces, run throughout their (3) nodes.
				for (int nodeId = 0; nodeId < 2 + (MESH_ORDER-1); nodeId++)
			        fprintf(fw, " %d", mesh->mBFaces[ bFaceNumber ].nodes[ nodeId ] + 1 );
            	// Skip line:
				fprintf(fw, "\n");
            }
		
		}
		
	}
	
	fprintf(fw,"$EndElements");
	fprintf(fw, "\n$PhysicalNames\n");
	fprintf(fw, "%d", mesh->sizes.zones);
	
	for (short i = 0; i < mesh->sizes.zones; i++) 
		fprintf(fw, "\n%d %d \"%s\"", mesh->zones[i].dim, mesh->zones[i].idx, mesh->zones[i].name); 
	
	fprintf(fw, "\n$EndPhysicalNames");
	
	fclose( fw );
    
	
	return;
}

