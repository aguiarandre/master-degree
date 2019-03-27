#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "native.h"
#include "BSpline.h"

double bsplineBasis( double u, int paramsId, int splineId, BSpline * spline, int m )
{
    double N, temp1, temp2;
    
    
    if ( m == 0 ) 
    {
        if (  ( u <= spline[ splineId ].knots[ paramsId + 1 ] ) && ( u >= spline[ splineId ].knots[ paramsId ] ) )
            N = 1.0;
        else 
            N = 0.0;
    }
    else
    {
        if ( fabs( spline[ splineId ].knots[ m + paramsId ] - spline [ splineId ].knots[ paramsId ] ) < EPS )
            temp1 = 0.0;
        else
            temp1 = ( u - spline[ splineId ].knots[ paramsId ] ) / ( spline[ splineId ].knots[ m + paramsId ] - spline[ splineId ].knots[ paramsId ] ) ;
        
        if ( fabs( spline[ splineId ].knots[ paramsId + m + 1 ] - spline[ splineId ].knots[ paramsId + 1 ] ) < EPS ) 
            temp2 = 0.0;
        else
            temp2 = ( spline[ splineId ].knots[ paramsId + m + 1] - u ) / ( spline[ splineId ].knots[ paramsId + m + 1 ] - spline[ splineId ].knots[ paramsId + 1 ] );
            
        N = temp1 * bsplineBasis( u , paramsId, splineId, spline, m-1 ) + temp2 * bsplineBasis( u, paramsId + 1, splineId, spline, m-1 ) ;
        
    }
    
    return N;
    
}

void reconstructSpline( BSpline * spline, int splineId, double parameter )
{
    double foo;
    int k, m;

    k = spline[ splineId ].k;
    m = spline[ splineId ].m;

   
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0, sum_den = 0.0;
        
        foo = spline[ splineId ].knots[0] + ( spline[splineId].knots[ 1 + k + m ] - spline[ splineId ].knots[0] ) * ( parameter );
        
        for (int paramsId = 0; paramsId <= k ; paramsId++ )
        {            
			
            // somatoria em k:
            double basis = bsplineBasis( foo, paramsId, splineId, spline, m ); 
            double weight = spline[ splineId ].weights[paramsId];
            
            sum_x = sum_x + spline[ splineId ].x_ctrl[paramsId] * weight * basis;
            sum_y = sum_y + spline[ splineId ].y_ctrl[paramsId] * weight * basis;
            sum_z = sum_z + spline[ splineId ].z_ctrl[paramsId] * weight * basis;
            
            sum_den = sum_den + weight * basis;
            
        }
        if ( fabs(sum_den) < EPS ) 
            printf("ERROR, sum_den = 0.0 \n");
        
        spline[ splineId ].x = sum_x / sum_den;
        spline[ splineId ].y = sum_y / sum_den;
        spline[ splineId ].z = sum_z / sum_den;
        
    
    return;
}

bool getParameter( double x, double y, BSpline * spline, int splineId, double * t  )
{
    double x1min, x1max, x2min, x2max;
    double t0, t1, tmed;
    double epsx = 1e-15, epsy = 0.001; // 0.02 for cylinder and NACA CHANGES
    t0 = 0.0; tmed = 0.5; t1 = 1.0;
    double diffx;
    int nIter = 0;

    
    
    while ( fabs(t0 - t1) > epsx && nIter < MAX_ITERATION )
    {
        if ( spline[ splineId ].k < 2 ) 
            return false; //TODO
        
        reconstructSpline( spline, splineId, t0 );
        x1min = spline[ splineId ].x;
        
        reconstructSpline( spline, splineId, tmed );
        x2min = spline[ splineId ].x;
        x1max = x2min;
        
        reconstructSpline( spline, splineId, t1 );
        x2max = spline[ splineId ].x;
        
        if ( (x >= x1min && x <= x1max) || (x >= x1max && x < x1min) )
        {
            //t0 = t0;
            t1 = tmed;
            tmed = (t0 + t1) / 2.0;
            diffx = fabs( x - x1max );
        }
        else if ( (x > x2min && x <= x2max) || (x >= x2max && x <= x2min) )
        {
            t0 = tmed;
            //t1 = t1;
            tmed = (t0 + t1) / 2.0;
            diffx = fabs( x - x1max );
        } 
        else // try near x1min and x2max
        {
            x1min = x1min-0.1;
            x2max = x2max+0.1;
            if ( (x >= x1min && x <= x1max) || (x >= x1max && x < x1min) )
            {
                //t0 = t0;
                t1 = tmed;
                tmed = (t0 + t1) / 2.0;
                diffx = fabs( x - x1max );
            }
            else if ( (x > x2min && x <= x2max) || (x >= x2max && x <= x2min) )
            {
                t0 = tmed;
                //t1 = t1;
                tmed = (t0 + t1) / 2.0;
                diffx = fabs( x - x1max );
            } 
            else
            {
                // false spline
                return false;
            }
        }
        
        nIter++;
    }

            
    reconstructSpline( spline, splineId, tmed );
    double yTest = spline[ splineId ].y;
    
    if ( (yTest > 0.0 && y < 0.0) || (yTest < 0.0 && y > 0.0) ) return false;
    
    if ( fabs( yTest - y ) < epsy )
    {
        *t = tmed;
        return true;
    }   
    
    return false;
}

void splineInit( BSpline * spline )
{
  for ( int i = 0; i < MAXSPLINES ; i++ )
  { 
      spline[i].x_ctrl = malloc( MAXPOINT * sizeof(double) );
      spline[i].y_ctrl = malloc( MAXPOINT * sizeof(double) );
      spline[i].z_ctrl = malloc( MAXPOINT * sizeof(double) );

      spline[i].knots = malloc( MAXPOINT * sizeof(double) );
      spline[i].weights = malloc( MAXPOINT * sizeof(double) );

  }
  
  return;
}

void initializeCurvedVars( MeshStruct * mesh )
{
    
    int nNodesFace = (mesh->sizes.faces + mesh->sizes.bfaces) * (MESH_ORDER-1);
    int nNodesCell = 4 + 4 * (MESH_ORDER-1) + (MESH_ORDER-1) * (MESH_ORDER-1);
    int nInternalNodes = mesh->sizes.elems * (MESH_ORDER-1) * (MESH_ORDER-1);
    
    // The total # of nodes is:
    int nNodes = mesh->sizes.nodes +    // normal nodes
                 nNodesFace +           // HO nodes 
                 nInternalNodes;        // Cell nodes
    
    mesh->vertices = realloc( mesh->vertices, nNodes * sizeof(mesh->vertices) );
    
    for ( int i = mesh->sizes.nodes; i < nNodes; i++)
        mesh->vertices[i] = malloc( 2 * sizeof(double) );  

    for ( int i = 0; i < mesh->sizes.faces; i++ )
        mesh->mFaces[i].nodes = realloc( mesh->mFaces[i].nodes, 2 + (MESH_ORDER-1) * sizeof(double) ); 
    
    for ( int i = 0; i < mesh->sizes.bfaces ; i++)
        mesh->mBFaces[i].nodes = realloc( mesh->mBFaces[i].nodes, 2 + (MESH_ORDER-1) * sizeof(double) );
    
    
    for ( int i = 0; i < mesh->sizes.elems; i++ )
        mesh->mElems[i].nodes = realloc ( mesh->mElems[i].nodes, nNodesCell * sizeof(double) );
    
    return;
}

void curveMesh( MeshStruct * mesh, BSpline * spline )
{
    // Based on the MESH_ORDER, we can realloc memory for some vars.
    initializeCurvedVars( mesh );

    reOrderBoundaries( mesh );
    
    /* 
     * Usual 1st order nodes are already generated. 
     * They finish their count at (mesh->sizes.nodes - 1)
     * So, the nodeNumber starts at: */
    
    int nodeNumber = mesh->sizes.nodes;
    
    nodeNumber = getWallNodes( spline, mesh, nodeNumber );
    nodeNumber = getInnerNodes( mesh, nodeNumber );
    nodeNumber = getBoundNodes( mesh, nodeNumber );
    nodeNumber = getCellNodes ( mesh, nodeNumber );
    
    getConnectivity( mesh );
    
    printf("Curved Mesh Generated ... \n");
    
    return;
}

void reOrderBoundaries( MeshStruct * mesh )
{
    // Only works for quadrangular elements.
    
    // run through the boundary faces.
    for (int bfaceId = 0; bfaceId < mesh->sizes.bfaces; bfaceId++)
    {
        // get its node in order.
        int node1, node2;
        node1 = mesh->mBFaces[ bfaceId ].nodes[0];
        node2 = mesh->mBFaces[ bfaceId ].nodes[1];
        int found = 0;
        
        // check whether you can find this bface in any element.
        for (int elemId = 0; elemId < mesh->sizes.elems; elemId++)
        {
            int a, b, c, d;
            a = mesh->mElems[ elemId ].nodes[0];
            b = mesh->mElems[ elemId ].nodes[1];
            c = mesh->mElems[ elemId ].nodes[2];
            d = mesh->mElems[ elemId ].nodes[3];
            
            if ( node1 == a )
                if ( node2 == b )
                {
                    found = 1; break;
                }
                    
            if ( node1 == b )
                if ( node2 == c )
                {
                    found = 1; break;
                }
                    
            if ( node1 == c )
                if ( node2 == d )
                {
                    found = 1; break;
                }
                    
            if ( node1 == d )
                if ( node2 == a )
                {
                    found = 1; break;
                }   
        }
        
        // Swap nodes:
        if ( !found )
        {
            mesh->mBFaces[ bfaceId ].nodes[0] = mesh->mBFaces[ bfaceId ].nodes[0]^mesh->mBFaces[ bfaceId ].nodes[1];
            mesh->mBFaces[ bfaceId ].nodes[1] = mesh->mBFaces[ bfaceId ].nodes[0]^mesh->mBFaces[ bfaceId ].nodes[1];
            mesh->mBFaces[ bfaceId ].nodes[0] = mesh->mBFaces[ bfaceId ].nodes[1]^mesh->mBFaces[ bfaceId ].nodes[0];
        }
    }
    
    return;
}

int getWallNodes( BSpline * spline, MeshStruct * mesh, int nodeNumber )
{
    // Generate HO nodes at wall boundaries:
    for( int zoneId = 0; zoneId < mesh->sizes.zones; zoneId++ )
    {
        // Find curved wall
        if ( mesh->zones[ zoneId ].bc_id == 1 ) 
        {
            int beg = mesh->zones[ zoneId ].start;
            int end = mesh->zones[ zoneId ].start + mesh->zones[ zoneId ].num_elems;
            
            // Run through its bfaces.
            for (int bfaceId = beg, node; bfaceId < end; bfaceId++) 
            {
                                
                double x1,x2,y1,y2;
                double t1, t2;
                int paramFound[2];
                paramFound[0] = -1; paramFound[1] = -1;
                
                // x,y of 1st node.
                node = mesh->mBFaces[ bfaceId ].nodes[0];
                x1 = mesh->vertices[ node ][0];
                y1 = mesh->vertices[ node ][1];
                
                
                // x,y of 2nd node.
                node = mesh->mBFaces[ bfaceId ].nodes[1];
                x2 = mesh->vertices[ node ][0];
                y2 = mesh->vertices[ node ][1];

                    
                // look up for these 2 nodes at all splines:
                for ( int splineId = 0; splineId < mesh->sizes.splines; splineId++ )
                {
                    if ( paramFound[0] == -1 )
                    {
                        if ( getParameter( x1, y1, spline, splineId, &t1 ) )
                            paramFound[0] = splineId;
                        else
                            paramFound[0] = -1;
                    }
                    if ( paramFound[1] == -1 )
                    {
                        if ( getParameter( x2, y2, spline, splineId, &t2 ) )
                            paramFound[1] = splineId;
                        else
                            paramFound[1] = -1;
                    }
                } // end lookup through splines
                
                double x, y;
                
                if ( paramFound[0] == -1 || paramFound[1] == -1 )
                {
                    printf("Parameter wasn't found in any spline. TODO\n");
                }
                else // if both parameters were found.
                {
                    
                    for (int nodeId = 0, splineId; nodeId < (MESH_ORDER-1); nodeId++) 
                    {   
                        splineId = projectNodes( t1, t2, paramFound, spline, nodeId , mesh);
                        
                        mesh->vertices[ nodeNumber ][0] = spline[ splineId ].x ;
                        mesh->vertices[ nodeNumber ][1] = spline[ splineId ].y ;
                        
                        mesh->mBFaces[ bfaceId ].nodes[ 2 + nodeId ] = nodeNumber;
                                                
                        nodeNumber++;
                    }
                    
                } // end if found 
                        
            } //end bfaceId				
            
        }
                
    }
    
    return nodeNumber;
    
}

int projectNodes( double t1, double t2, int * paramFound, BSpline * spline, int nodeId, MeshStruct * mesh )
{
    // paramFound[0] = splineId of parameter t1
    // paramFound[1] = splineId of parameter t2
    // nodeId = used as index of lookUpArray.
    int splineId1 = paramFound[0];
    int splineId2 = paramFound[1];
    
    int splineId; 
    double t;
    
    if ( MESH_ORDER < 1 )
    { 
        printf("Working curved mesh orders: 2, 3, 4\n"); exit(1);
    }
    
   /* double * lookUpArray = (double*) malloc( (MESH_ORDER - 1) * sizeof(double) );
    
    switch( MESH_ORDER )
    {
        case 2:
        lookUpArray[0] = 0.5; 
        break;
        
        case 3:
        lookUpArray[0] = 0.276393205;
        lookUpArray[1] = 0.723606797;
        break;
        
        case 4:
        lookUpArray[0] = 0.172673164;
        lookUpArray[1] = 0.5;
        lookUpArray[2] = 0.827326836;
        break;
        
        default:
        printf("Working curved mesh orders: 2, 3, 4\n");
        exit(1);
    }    */
    
    if ( splineId1 == splineId2 )
    {
                 
        t = t1 + mesh->lookUpArray[ nodeId ] * (t2 - t1);
            
        splineId = splineId1; // does not matter. 
        if( fabs( t - 0.83402519391394792 ) < EPS ) 
            printf("AQUI");
        
        reconstructSpline( spline, splineId, t );
                                              
    }
    else 
    {
        // There are 4 cases:
        if ( t1 > 0.7 && t2 < 0.3 ) // they are in correct usual order
        {
            double dist = ( (1.0 + t2) - t1 );
            t = t1 + mesh->lookUpArray[ nodeId ] * dist;
            if ( t > 1.0 )
            {
                t = t - 1.0; splineId = splineId2; 
            }
            else 
                splineId = splineId1;

        }
        else if ( t1 > 0.7 && t2 > 0.7 )
        {
            double dist = ( 1.0 + (1.0 - t2) - t1 );
            t = t1 + mesh->lookUpArray[ nodeId ] * dist;
            if ( t > 1.0 )
            {
                t = t - 1.0;
                t = 1.0 - t; 
                splineId = splineId2; 
            }
            else 
                splineId = splineId1;
        }
        else if ( t1 < 0.3 && t2 < 0.3 )
        {
            double dist = t1 + t2;
            t = t1 - mesh->lookUpArray[ nodeId ] * dist;
            if ( t < 0.0 )
            {
                t = -t;
                splineId = splineId2;
            }
            else
                splineId = splineId1;
        }
        else if (t1 < 0.3 && t2 > 0.7 )
        {
            double dist = ( ( 1.0 - t2 ) + t1 );
            t = t1 - mesh->lookUpArray[ nodeId ] * dist;
            if ( t < 0.0 )
            {
                t = 1.0 + t;
                splineId = splineId2;
            }
            else 
                splineId = splineId1;
        }
        else
        {
            printf("There's probably an error. You're finding too separated parameters within different splines\n");
            splineId = splineId2; 
        }
        
       // if( fabs( t - 0.83402519391394792 ) < EPS ) 
         //   printf("AQUI");
        
        reconstructSpline( spline, splineId, t );

    }
    
    //free( lookUpArray );
    
    
// A solucao (x,y,z) esta dentro de spline[ splineId ].x,y,z

return splineId;
}

int getInnerNodes( MeshStruct * mesh, int nodeNumber )
{
    


    for (int faceId = 0; faceId <  mesh->sizes.faces; faceId++) 
    { 

        double x1,x2,y1,y2, avg_x, avg_y;

        for (int nodeId = 0; nodeId < (MESH_ORDER-1); nodeId++)
        {
            x1 = mesh->vertices[ mesh->mFaces[faceId].nodes[0] ][0]; 	// vertices[node#][x or y]
            y1 = mesh->vertices[ mesh->mFaces[faceId].nodes[0] ][1]; 	// mesh->mFaces[ facenumber ].nodes[0,1]
            
            x2 = mesh->vertices[ mesh->mFaces[faceId].nodes[1] ][0];
            y2 = mesh->vertices[ mesh->mFaces[faceId].nodes[1] ][1];
            
            avg_x = x1 + (x2 - x1) * mesh->lookUpArray[ nodeId ];
            avg_y = y1 + (y2 - y1) * mesh->lookUpArray[ nodeId ];
            
            mesh->mFaces[ faceId ].nodes[ 2 + nodeId ] = nodeNumber; //até int_faces são os pontos médios normais.
                    
            mesh->vertices[ nodeNumber ][0] = avg_x;
            mesh->vertices[ nodeNumber ][1] = avg_y;
            
            nodeNumber++;
            
        }

    } 
 //   free( lookUpArray );
    return nodeNumber;
}

int getBoundNodes( MeshStruct * mesh, int nodeNumber ) 
{
    /*double * lookUpArray = (double*) malloc( (MESH_ORDER - 1) * sizeof(double) );
    
    switch( MESH_ORDER )
    {
        case 2:
        lookUpArray[0] = 0.5; 
        break;
        
        case 3:
        lookUpArray[0] = 0.276393205;
        lookUpArray[1] = 0.723606797;
        break;
        
        case 4:
        lookUpArray[0] = 0.172673164;
        lookUpArray[1] = 0.5;
        lookUpArray[2] = 0.827326836;
        break;
        
        default:
        printf("Working curved mesh orders: 2, 3, 4\n");
        exit(1);
    }    */

    
    for ( int zoneId = 0; zoneId < mesh->sizes.zones; zoneId++ )
    {
        if ( mesh->zones[ zoneId ].is_boundary == true && mesh->zones[ zoneId ].bc_id != 1 )
        {
            int beg = mesh->zones[ zoneId ].start;
            int end = beg + mesh->zones[ zoneId ].num_elems;
            
            for (int bfaceId = beg; bfaceId < end ; bfaceId++)
            {
               double x1,x2,y1,y2, avg_x, avg_y;

                for (int nodeId = 0; nodeId < (MESH_ORDER-1); nodeId++)
                {
                    int node1 = mesh->mBFaces[ bfaceId ].nodes[0];
                    int node2 = mesh->mBFaces[ bfaceId ].nodes[1];
                
                    x1 = mesh->vertices[ node1 ][0]; 	
                    y1 = mesh->vertices[ node1 ][1]; 	
                    
                    x2 = mesh->vertices[ node2 ][0];
                    y2 = mesh->vertices[ node2 ][1];
                    
                    avg_x = x1 + (x2 - x1) * mesh->lookUpArray[ nodeId ];
                    avg_y = y1 + (y2 - y1) * mesh->lookUpArray[ nodeId ];
                    
                    mesh->mBFaces[bfaceId].nodes[ 2 + nodeId ] = nodeNumber; 
                            
                    mesh->vertices[ nodeNumber ][0] = avg_x;
                    mesh->vertices[ nodeNumber ][1] = avg_y;
                    
                    nodeNumber++;
                    
                }
                
            }
            
        }
        
    }
    
    //free( lookUpArray );
    
    return nodeNumber;
}

int getCellNodes ( MeshStruct * mesh, int nodeNumber )
{
    
    // Number of high-order nodes inside each cell:
    int nCellNodes = mesh->sizes.elems * (MESH_ORDER-1) * (MESH_ORDER-1);
    int nodeIdLast;
    int nBegin = 4 + 4 * (MESH_ORDER-1) ;
    // Loop through each cell:
    
    for (int elemId = 0; elemId < mesh->sizes.elems; elemId++)
    {
        
            // Start calculating the position of the nodes insice each cell.
            
            // 1st, so find the corresponding high-order nodes:
            // for instance, for a = mesh->mElem[ elemId ].faces[ 0 ] is opposed to 
            //                   b = mesh->mElem[ elemId ].faces[ 2 ].
            // Get the 1st node of a and the last node of b. They are the ones to be used.
            // Get the 2nd node of a and the last to last node of b. They are the ones.
            // So run through nodes of a up, and through nodes of b down.
            
           // Run through the faces checking whether they're boundaries:
           int faceNumber1isBoundary = 0;
           int faceNumber2isBoundary = 0;
           int faceNumber3isBoundary = 0;
           int faceNumber4isBoundary = 0;
           
           int faceNumber1 = mesh->mElems[ elemId ].faces[ 0 ];
           int faceNumber2 = mesh->mElems[ elemId ].faces[ 1 ];
           int faceNumber3 = mesh->mElems[ elemId ].faces[ 2 ];
           int faceNumber4 = mesh->mElems[ elemId ].faces[ 3 ];

            if ( faceNumber1 < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber1 ].elem_left == elemId )// is boundary
                faceNumber1isBoundary = 1;
            if ( faceNumber2 < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber2 ].elem_left == elemId )// is boundary
                faceNumber2isBoundary = 1;
            if ( faceNumber3 < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber3 ].elem_left == elemId )// is boundary
                faceNumber3isBoundary = 1;
            if ( faceNumber4 < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber4 ].elem_left == elemId )// is boundary
                faceNumber4isBoundary = 1;
                
                
        // FIXME: the (curved) wall nodes should be from 'previous' position - not curved position.
        //        You'll have to calculate them using the 1st order nodes and lookUpArray[ nodeId ]
        //        Otherwise nodes may not even create lines that intersect, actually.
             
        for (int nodeIdVert = 0; nodeIdVert < MESH_ORDER-1; nodeIdVert++)
        {
            // Run in horizontal
            for ( int nodeId = 0; nodeId < MESH_ORDER-1; nodeId++ )
            {
               int node1, node2, node3, node4;
               double x1, y1, x2, y2, x3, y3, x4, y4, x_avg, y_avg, t, vx, vy, wx, wy;
                      
               // Test whether this face should be inverted or not.
               if ( mesh->mFaces[ faceNumber1 ].elem_left == elemId )
               {
                    if ( !faceNumber1isBoundary )
                        node1 = mesh->mFaces[ faceNumber1 ].nodes[ 2 + nodeId ];
                    else
                    {
                        // if it is a boundary, you should not use the projected node. 
                        // Use the previous position instead.
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber1 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber1].nodes[1];
                        int nodePos = ( MESH_ORDER - 2 ) - nodeId;
                        
                        x1 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y1 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                        
                        
                    }   //node1 = mesh->mBFaces[ faceNumber1 ].nodes[ MESH_ORDER - nodeId ]; // bface is inverted!
                }
                else // if I've been through this face before:
                {
                    if ( !faceNumber1isBoundary )
                        node1 = mesh->mFaces[ faceNumber1 ].nodes[ MESH_ORDER - nodeId ];
                    else
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber1 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber1 ].nodes[1];
                        int nodePos = nodeId;

                        x1 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y1 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                        
                    }
                        //node1 = mesh->mBFaces[ faceNumber1 ].nodes[ 2 + nodeId ];                   
                }
                // Opposing face:
                if ( mesh->mFaces[ faceNumber3 ].elem_left == elemId )
                {    
                    if ( !faceNumber3isBoundary )           
                        node3 = mesh->mFaces[ faceNumber3 ].nodes[ MESH_ORDER - nodeId ];
                    else 
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber3 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber3 ].nodes[1];
                        int nodePos = nodeId;

                        x3 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y3 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                        
                    }
                        //node3 = mesh->mBFaces[ faceNumber3 ].nodes[ 2 + nodeId ];
                }
                else
                {
                    if ( !faceNumber3isBoundary )           
                        node3 = mesh->mFaces[ faceNumber3 ].nodes[ 2 + nodeId ];
                    else 
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber3 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber3 ].nodes[1];
                        int nodePos = ( MESH_ORDER - 2 ) - nodeId;


                        x3 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y3 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                    
                    }
                        //node3 = mesh->mBFaces[ faceNumber3 ].nodes[ MESH_ORDER - nodeId ];   
                }
                
                
                // These two should be fixed at nodeIdVert
                if ( mesh->mFaces[ faceNumber2 ].elem_left == elemId )
                {
                    if ( !faceNumber2isBoundary )
                        node2 = mesh->mFaces[ faceNumber2 ].nodes[ 2 + nodeIdVert ];
                    else
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber2 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber2 ].nodes[1];
                        int nodePos = (MESH_ORDER - 2)  - nodeIdVert;

                        
                        x2 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y2 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                        
                    }
                        //node2 = mesh->mBFaces[ faceNumber2 ].nodes[ MESH_ORDER - nodeIdVert ];
                }
                else
                {
                    if ( !faceNumber2isBoundary )
                        node2 = mesh->mFaces[ faceNumber2 ].nodes[ MESH_ORDER - nodeIdVert ];
                    else
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber2 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber2 ].nodes[1];
                        int nodePos = nodeIdVert;
                        
                        x2 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y2 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                    }
                        //node2 = mesh->mBFaces[ faceNumber2 ].nodes[ 2 + nodeIdVert ];
                }
                
                if ( mesh->mFaces[ faceNumber4 ].elem_left == elemId ) 
                {
                    if ( !faceNumber4isBoundary )
                        node4 = mesh->mFaces[ faceNumber4 ].nodes[ MESH_ORDER - nodeIdVert ];
                    else
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber4 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber4 ].nodes[1];
                        int nodePos = nodeIdVert;
                        
                        x4 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y4 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                    }
                       // node4 = mesh->mBFaces[ faceNumber4 ].nodes[ 2 + nodeIdVert ];
                }
                else
                {
                    if ( !faceNumber4isBoundary )
                        node4 = mesh->mFaces[ faceNumber4 ].nodes[ 2 + nodeIdVert ];
                    else
                    {
                        int nodeLeft, nodeRight;
                        nodeLeft = mesh->mBFaces[ faceNumber4 ].nodes[0];
                        nodeRight = mesh->mBFaces[ faceNumber4 ].nodes[1];
                        int nodePos = (MESH_ORDER -2)  - nodeIdVert;
                        
                        x4 = mesh->vertices[ nodeLeft ][0] + ( mesh->vertices[ nodeRight ][0] - mesh->vertices[ nodeLeft ][0] ) * mesh->lookUpArray[ nodePos ];
                        y4 = mesh->vertices[ nodeLeft ][1] + ( mesh->vertices[ nodeRight ][1] - mesh->vertices[ nodeLeft ][1] ) * mesh->lookUpArray[ nodePos ];
                    }
                        //node4 = mesh->mBFaces[ faceNumber4 ].nodes[ MESH_ORDER - nodeIdVert ];
                }
                
                
               // Doing segment intersection here.
         
               // get positions of these nodes:
               if ( !faceNumber1isBoundary )
               {
                   x1 = mesh->vertices[ node1 ][0];
                   y1 = mesh->vertices[ node1 ][1];
               }
               if ( !faceNumber2isBoundary )
               {
                   x2 = mesh->vertices[ node2 ][0];
                   y2 = mesh->vertices[ node2 ][1];
               }
               if ( !faceNumber3isBoundary )
               {
                   x3 = mesh->vertices[ node3 ][0];
                   y3 = mesh->vertices[ node3 ][1];
               }
               if ( !faceNumber4isBoundary )
               {
                   x4 = mesh->vertices[ node4 ][0];
                   y4 = mesh->vertices[ node4 ][1];
               }
               
               double m1, m2;
              
               if ( fabs( x1 - x3 ) < EPS )
               {
                   x_avg = x1;
                   y_avg = y2;
               }
               else
                    m1 = (y1-y3) / (x1-x3);
               
               if( fabs( x4 - x2 ) < EPS )
               {
                   x_avg = x2;
                   y_avg = y1;
               }
               else
                   m2 = (y4-y2) / (x4-x2);
              
              if ( fabs( x1-x3 ) > EPS && fabs( x4-x2 ) > EPS )
              {
                   x_avg = ( m1 * x1 - m2 * x2 + y2 - y1 )/ ( m1 - m2 );
                   y_avg = m1 * (x_avg - x1) + y1; 
              } 
               /*
               vx = x2 - x1;
               vy = y2 - y1;
               
               wx = x4;
               wy = y4;
               
               if ( fabs( x3 * y4 - y3 * x4 ) > EPS )
                    t = ( vx * wy - vy * wx ) / ( x3*y4 - y3*x4 );
               else
                    printf("Division by zero when trying to intersect lines in cell\n");
                
                if ( t > 1.0 || t < 0.0 )
                    printf("Parameter of segment intersection invalid, t = %lf\n", t);
               
               
               // Now, the position of the inner node shall be:
               x_avg = x1 + t * x3 ;
               y_avg = y1 + t * y3 ;
               */
               
               int whichId = nBegin + nodeId + nodeIdVert * (MESH_ORDER-1);
               
               mesh->mElems[ elemId ].nodes[ whichId ] = nodeNumber;
               
               mesh->vertices[ nodeNumber ][0] = x_avg;
               mesh->vertices[ nodeNumber ][1] = y_avg;
               
               nodeNumber++;
                 
            }    
        }
    }
    
    
    //free( lookUpArray );
    
    return nodeNumber;
}

void getConnectivity( MeshStruct * mesh )
{
    
    // Populate mesh->mElems[ elemId ].nodes[4 - end]
    int * beenHere = calloc( mesh->sizes.faces , sizeof(int) );
    
    // Run through elements:
    for (int elemId = 0; elemId < mesh->sizes.elems; elemId++)
    {
        // Run through faces of this element:
        for (int faceId = 0 ; faceId < 4 ; faceId++)
        {
            // This face has the number:
            int faceNumber = mesh->mElems[ elemId ].faces[ faceId ];
            
            // Check whether you passed here before. If so, you have to populate 
            // the nodes in reverse order!
            if ( beenHere[ faceNumber ] )
            {
                // Run through this face's nodes and populate mElems.nodes
                for (int nodeId = 0; nodeId < (MESH_ORDER-1); nodeId++ )
                {
                    if ( faceNumber < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber ].elem_left == elemId )
                        mesh->mElems[ elemId ].nodes[ 4 + faceId * (MESH_ORDER-1) + nodeId ] = mesh->mBFaces[ faceNumber ].nodes[ 2 + nodeId ];
                    else
                    {
                        mesh->mElems[ elemId ].nodes[ 4 + faceId * (MESH_ORDER-1) + nodeId ] = mesh->mFaces[ faceNumber ].nodes[ MESH_ORDER - nodeId ];
                        beenHere[ faceNumber ] = 1;

                    }  
                    
                }
            }
            else 
            {
                // Run through this face's nodes and populate mElems.nodes
                for (int nodeId = 0; nodeId < (MESH_ORDER-1); nodeId++ )
                {
                    if ( faceNumber < mesh->sizes.bfaces && mesh->mBFaces[ faceNumber ].elem_left == elemId )
                        mesh->mElems[ elemId ].nodes[ 4 + faceId * (MESH_ORDER-1) + nodeId ] = mesh->mBFaces[ faceNumber ].nodes[ 2 + nodeId ];
                    else
                    {
                        mesh->mElems[ elemId ].nodes[ 4 + faceId * (MESH_ORDER-1) + nodeId ] = mesh->mFaces[ faceNumber ].nodes[ 2 + nodeId ];
                        beenHere[ faceNumber ] = 1;
                    }   
                }
            }
            

        }
        
    }  
    
    free( beenHere );
    
    // Post-processing step to re-order internal nodes within a cell.
    if ( MESH_ORDER == 3 )
    {
        for (int elemId = 0; elemId < mesh->sizes.elems; elemId++ )
        {
            // swap nodes 14 and 15.
            mesh->mElems[ elemId ].nodes[ 14 ] = mesh->mElems[ elemId ].nodes[ 14 ] ^ mesh->mElems[ elemId ].nodes[ 15 ];
            mesh->mElems[ elemId ].nodes[ 15 ] = mesh->mElems[ elemId ].nodes[ 14 ] ^ mesh->mElems[ elemId ].nodes[ 15 ];
            mesh->mElems[ elemId ].nodes[ 14 ] = mesh->mElems[ elemId ].nodes[ 14 ] ^ mesh->mElems[ elemId ].nodes[ 15 ];
        }
    }
    else if ( MESH_ORDER == 4 )
    {
        for (int elemId = 0; elemId < mesh->sizes.elems; elemId++ )
        {
            // swap lots of nodes inside node:
            // node16 ok
            // 17 goes to 20.
            int temp = mesh->mElems[ elemId ].nodes[20];
            mesh->mElems[ elemId ].nodes[20] = mesh->mElems[ elemId ].nodes[ 17 ];
            
            // 20 goes to 24
            int temp2 = mesh->mElems[ elemId ].nodes[24] ;
            mesh->mElems[ elemId ].nodes[24] = temp;
            
            // 24 goes to 18
            temp = mesh->mElems[ elemId ].nodes[18];
            mesh->mElems[ elemId ].nodes[18] = temp2;
            
            // 18 goes to 17
            temp2 = mesh->mElems[ elemId ].nodes[ 17 ];
            mesh->mElems[ elemId ].nodes[17] = temp;
            
            // 16, 17, 18 done.
            
            // 19 goes to 23.
            temp = mesh->mElems[ elemId ].nodes[ 23 ];
            mesh->mElems[ elemId ].nodes[23] = mesh->mElems[ elemId ].nodes[19];
            
            // 23 goes to 22
            temp2 = mesh->mElems[ elemId ].nodes[22];
            mesh->mElems[ elemId ].nodes[22] = temp;
            
            // 22 goes to 19
            mesh->mElems[ elemId ].nodes[19] = temp2;
            
            // 16, 17, 18, 19, 20, 21, 22, 23, 24 done
        }
    }
    else if ( MESH_ORDER > 4 )
    {
        printf("Re-order not implemented yet for mesh-order = %d. Connectivity will be f*cked up\n", MESH_ORDER);
        exit(0);
    }
    
    return;
    
}





