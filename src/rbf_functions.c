/** DONE: Up to date, we have implemented 
 * a list of nodes that have prescribed 
 * movement, i.e., in this case, boundary
 * nodes (nodeNumber, x,y,z). They are 
 * stored at: mesh->rbf->fixedNodes
 * 15/02/2017
* */

/** DONE: Have implemented a list of nodes 
 * that are going to be moved as a result 
 * of the RBF. These nodes and vertices are 
 * stored at: mesh->rbf->movingNodes in preety
 * much the same way fixedNodes are.
 * */

/** DONE: I have to implement the same list,
 * however, with the previous location of those
 * nodes. For that, I'll have to calculate the 
 * average between 1st order nodes.
 * */


/** Sequence for solving RBF */

/* CREATING NECESSARY DATA {DATA STRUCTURE REQUIREMENTS} */
/** Create a list with x,y,z of points that will have a prescribed movement, before the movement {rbf.fixedNodesBefore}. <OK> */
/** Create a list with x,y,z of points that will not have a prescribed movement, before the movement {rbf.movingNodes}. <OK> */
/** Perform the desired displacement on the nodes with prescribed movement. {This is what curveMesh does} <OK> */
/** Create a list with x,y,z of points that had a prescribed movement, after the movement. {rbf.fixedNodes} <OK> */ 

/* MATRICES FROM THESE DATA */
/** Create a  C matrix, with size (3 + nPrescribed rows, 3 + nPrescribed columns)  <OK> */
/** Create an A matrix, with size ( nNonPrescribed rows, 3 + nPrescribed columns)  <OK> */

/* POPULATING THOSE MATRICES - C MATRIX */
/** The first row and column of C is 1 except for C(0:2, 0:2) = 0; <OK> */
/** C(1,3:end) = xFixedAfterMovement; C(2,3:end) = yFixedAfterMovement <OK> */
/** C(3:end,1) = xFixedAfterMovement; C(3:end,2) = yFixedAfterMovement - NOTE THAT Z CANNOT BE HERE.  <OK>*/
/** Interior of Matrix is: C(i+3,j+3) = phi( fixedAfter(i,xyz), fixedAfter(j,xyz) )  <OK> */

/* POPULATING THOSE MATRICES - A MATRIX */
/** A Matrix: A(:,0) = 1, i.e., first column = 1; A(:,1) = xMoving; A(:,2) = yMoving;  <OK> */
/** Then, A(i, j +3 ) = phi( movingBefore(i,xyz) , fixedAfter(j,xyz) ), FOR i = 1,nPrescribed. j = 1, nNonPrescribed  <OK> */

/* OBTAINING RESULTS */
/** Store the inverse(C) <OK> */
/** Multiply H = A * inverse(C)  <OK> */

/* CALCULATE DISPLACEMENTS */
/** First calculate xDesloc = xFixedAfter - xFixedBefore, yDesloc = yFixedAfter - yFixedBefore  <TODO> */
/** This will produce a vector of displacements:   */
/** U_fixed = zeros( 3+nPrescribed rows, 3 columns); U_fixed(4:end,:) = [ xDesloc(:), yDesloc(:) ];  */
/** To calculate the displacements of the moving nodes:  */
/** U_moving = H * U_fixed  */
/** Now: x_moving(1) = x_moving(1) + U_moving(1) and so on.    */

#include "rbf_functions.h"

#include "native.h"
#include <stdlib.h> // malloc
#include <string.h>
#include <stddef.h>
#include <math.h>


#define ind2d(i,j) i*tam+j


 

void propagateDeformation( MeshStruct * mesh )
{ 
	printf("\n");
	printf("Propagating deformation to the interior ...\n");
	
    /** ALLOCATE MEMORY */
	initializeData( mesh );
	
	/* CREATING NECESSARY DATA */
	getFixedNodes( mesh );
	
	
	getMovingNodes( mesh );
	
	
	/* MATRICES FROM THESE DATA */
	populateMatrices( mesh );
	
	/* CALCULATE DISPLACEMENTS */
	applyDisplacement( mesh );

	/* WRITE RESULTING MESH */
	writeMesh( mesh );
	
	/* FREE UP DATA */
	destroyData( mesh );
	
	printf("\n");
	printf("nFixedNodes  = %d\nnMovingNodes = %d\n", mesh->rbf->nFixedNodes, mesh->rbf->nMovingNodes );
	
	return;
}


void inverse(double* A, int N)
{
    int * IPIV = malloc( N * sizeof(int) ) ;
    int LWORK = N*N;
    double * WORK = malloc( LWORK * sizeof(double) );
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

void initializeData( MeshStruct * mesh )
{
	int nBNodes = MESH_ORDER * mesh->sizes.bfaces;   
    int nNodesFace = (mesh->sizes.faces + mesh->sizes.bfaces) * (MESH_ORDER-1);
    int nNodesCell = 4 + 4 * (MESH_ORDER-1) + (MESH_ORDER-1) * (MESH_ORDER-1);
    int nInternalNodes = mesh->sizes.elems * (MESH_ORDER-1) * (MESH_ORDER-1);
    
    // The total # of nodes is:
    int nNodes = mesh->sizes.nodes +    // normal nodes
                 nNodesFace +           // HO nodes 
                 nInternalNodes;        // Cell nodes
                 
	pFixedNodeStruct fixedNodesBefore;
	fixedNodeStruct fixedNodes;
	movingNodeStruct movingNodes;
	
	mesh->rbf->isSet = calloc( nNodes , sizeof( int ) );
	
	mesh->rbf->fixedNodesBefore = &fixedNodesBefore;
	mesh->rbf->fixedNodesBefore = calloc( nBNodes , sizeof ( fixedNodeStruct ) );
	
	mesh->rbf->fixedNodes = &fixedNodes;
	mesh->rbf->fixedNodes = calloc( nBNodes , sizeof ( fixedNodeStruct ) );
	
	mesh->rbf->movingNodes = &movingNodes;
	mesh->rbf->movingNodes = calloc( nNodes - nBNodes , sizeof ( movingNodeStruct ) );
	
	return;
}


void getFixedNodes( MeshStruct * mesh )
{
    int startingElem, endingElem;
	int nodeNumber = 0, nodeCounter = 0;
	
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
				// collect all three nodes it has.
				
				for (int nodeId = 0; nodeId < 2 + (MESH_ORDER-1); nodeId++) 
				{
					// store the number of each node.
					nodeNumber = mesh->mBFaces[ bFaceNumber ].nodes[ nodeId ] ;
					
					// Now, only enter here if you have not previously entered!
					if ( !mesh->rbf->isSet[ nodeNumber ] )
					{
						mesh->rbf->isSet[ nodeNumber ] = 1; // leave a flag saying you have been here.
						
						// Store data:
						mesh->rbf->fixedNodes[ nodeCounter ].nodeNumber = nodeNumber;
						mesh->rbf->fixedNodes[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0];
						mesh->rbf->fixedNodes[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1];
						mesh->rbf->fixedNodes[ nodeCounter ].z = 0.0;
						
						if (nodeId < 2)
						{
							// If it was a normal node, then it stays at the same position as before.
							mesh->rbf->fixedNodesBefore[ nodeCounter ].nodeNumber = nodeNumber;
							mesh->rbf->fixedNodesBefore[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0];
							mesh->rbf->fixedNodesBefore[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1];
							mesh->rbf->fixedNodesBefore[ nodeCounter ].z = 0.0;
						}
						else if (nodeId >= 2)
						{
							// However, if it is a high-order node:
							// Take the average of nodeId 0 and 1. These were the previous position of projected nodes. 
							int nodeNumber0 = mesh->mBFaces[ bFaceNumber ].nodes[ 0 ];
							int nodeNumber1 = mesh->mBFaces[ bFaceNumber ].nodes[ 1 ];
							
							
							double x_value0 = mesh->vertices[ nodeNumber0 ][0];
							double x_value1 = mesh->vertices[ nodeNumber1 ][0];
							double avg_x = x_value0 + ( x_value1 - x_value0 ) * mesh->lookUpArray[ nodeId - 2 ];		
							
							double y_value0 = mesh->vertices[ nodeNumber0 ][1];
							double y_value1 = mesh->vertices[ nodeNumber1 ][1];
							double avg_y = y_value0 + ( y_value1 - y_value0 ) * mesh->lookUpArray[ nodeId - 2 ];
							
							// Now you finally have the previous position of the projected nodes. Go on and store it.
							mesh->rbf->fixedNodesBefore[ nodeCounter ].nodeNumber = nodeNumber;
							mesh->rbf->fixedNodesBefore[ nodeCounter ].x = avg_x;
							mesh->rbf->fixedNodesBefore[ nodeCounter ].y = avg_y;
							mesh->rbf->fixedNodesBefore[ nodeCounter ].z = 0.0;
							
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

void destroyData( MeshStruct * mesh )
{
	/** FIXME */
	int nM = mesh->rbf->nMovingNodes;
	int nP = mesh->rbf->nFixedNodes; 
	
	for (int i = 0; i < 3 + nP ; i++ )
	{
		free ( mesh->rbf->C[i] );
		free ( mesh->rbf->cInverse[i] );
	}
	for (int i = 0; i < nM ; i++ )
	{
		free ( mesh->rbf->A[i] );
		free ( mesh->rbf->H[i] );
	}
	
	free( mesh->rbf->A );
	free( mesh->rbf->C );
	free( mesh->rbf->cInverse );
	free( mesh->rbf->H );
	
	free( mesh->rbf->isSet );
	free( mesh->rbf->fixedNodes );
	free( mesh->rbf->fixedNodesBefore );
	free( mesh->rbf->movingNodes );
		
	
	return;
}

void getMovingNodes( MeshStruct * mesh )
{
	int nodeNumber = 0, nodeCounter = 0;
	
	// Run through all faces.
	for (int faceNumber = 0; faceNumber < mesh->sizes.faces ; faceNumber++)
	{
		// Also through all nodes of those faces ...
		for (int nodeId = 0; nodeId < 2 + (MESH_ORDER-1); nodeId++) 
		{
			nodeNumber = mesh->mFaces[ faceNumber ].nodes[ nodeId ] ;
			
			// If you haven't stored it, it means it is an interior moving node.
			// Note: if parallelization is being done, make sure getFixedNodes is performed before getMovingNodes function.
			if ( !mesh->rbf->isSet[ nodeNumber ] )
			{
				mesh->rbf->isSet[ nodeNumber ] = 1;
				mesh->rbf->movingNodes[ nodeCounter ].nodeNumber = nodeNumber;
				mesh->rbf->movingNodes[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0];
				mesh->rbf->movingNodes[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1];
				mesh->rbf->movingNodes[ nodeCounter ].z = 0.0;

				nodeCounter++;
				
			} // End if node isSet 
			
		} // End nodeId loop
				
	} // End for faces 
			
	// Run trhough all the elements and get center node.
	for (int elemNumber = 0; elemNumber < mesh->sizes.elems ; elemNumber++ )
	{
		for (int nodeId = 0; nodeId < (MESH_ORDER-1) * (MESH_ORDER-1); nodeId++)
		{
			nodeNumber = mesh->mElems[ elemNumber ].nodes[4 + 4 * (MESH_ORDER-1) + nodeId]; 
			
			if ( !mesh->rbf->isSet[ nodeNumber ] )
			{
				mesh->rbf->isSet[ nodeNumber ] = 1;
				mesh->rbf->movingNodes[ nodeCounter ].nodeNumber = nodeNumber;
				mesh->rbf->movingNodes[ nodeCounter ].x = mesh->vertices[ nodeNumber ][0];
				mesh->rbf->movingNodes[ nodeCounter ].y = mesh->vertices[ nodeNumber ][1];
				mesh->rbf->movingNodes[ nodeCounter ].z = 0.0;
				
				nodeCounter++;

			} // End if node isSet 
		}
	}
	
	// Store number of moving nodes.
	mesh->rbf->nMovingNodes = nodeCounter;
	
	return;
}

void populateMatrices( MeshStruct * mesh )
{
	int nP = mesh->rbf->nFixedNodes;
	int nM = mesh->rbf->nMovingNodes;

    // FIXME - 2D vs 3D	
	mesh->rbf->C = calloc( (3 + nP) , sizeof( double * ) );
	for (int i = 0; i < 3 + nP ; i++ )
		mesh->rbf->C[i] = calloc( (3 + nP) , sizeof( double ) );
	
	mesh->rbf->cInverse = calloc( (3 + nP) , sizeof( double * ) );
	for (int i = 0; i < 3 + nP ; i++ )
		mesh->rbf->cInverse[i] = calloc( (3 + nP) , sizeof( double ) );
	
	mesh->rbf->A = calloc( nM , sizeof( double * ) );
	for (int i = 0; i < nM ; i++ )
		mesh->rbf->A[i] = calloc( (3 + nP) , sizeof( double * ) );
		
	mesh->rbf->H = calloc( nM , sizeof( double * ) );
	for (int i = 0; i < nM ; i++ )
		mesh->rbf->H[i] = calloc( (3 + nP) , sizeof( double * ) );		
	
	// Populate C: 
	// First 3x3 subMatrix
	for (int i = 0; i < 3 ; i++)
		for (int j = 0; j < 3 ; j++ )
			mesh->rbf->C[i][j] = 0.0;
	
	// First row and column	(1)
	for (int i = 3, j = 0; i < 3 + nP; i++ )
		mesh->rbf->C[i][j] = 1.0;
		
	for (int j = 3, i = 0; j < 3 + nP; j++ )
		mesh->rbf->C[i][j] = 1.0;
	
	// Second row and column (x)
	for (int i = 3, j = 1; i < 3 + nP; i++ )
		mesh->rbf->C[i][j] = mesh->rbf->fixedNodesBefore[ i - 3 ].x;
	
	for (int j = 3, i = 1; j < 3 + nP; j++ )
		mesh->rbf->C[i][j] = mesh->rbf->fixedNodesBefore[ j - 3 ].x;
	
	// Third row and column (y)
	for (int i = 3, j = 2; i < 3 + nP; i++ )
		mesh->rbf->C[i][j] = mesh->rbf->fixedNodesBefore[ i - 3 ].y;
	
	for (int j = 3, i = 2; j < 3 + nP; j++ )
		mesh->rbf->C[i][j] = mesh->rbf->fixedNodesBefore[ j - 3 ].y;
	

	// Rest of the matrix positions: phi( fixedNodesBefore , fixedNodesBefore ) 
	for (int i = 3; i < 3 + nP; i++ )
	{
		for (int j = 3 ; j < 3 + nP ; j++ )
		{
			
			double v1[2], v2[2];
			
			// x, y
			v1[0] = mesh->rbf->fixedNodesBefore[ i-3 ].x; v1[1] = mesh->rbf->fixedNodesBefore[ i-3 ].y;
			v2[0] = mesh->rbf->fixedNodesBefore[ j-3 ].x; v2[1] = mesh->rbf->fixedNodesBefore[ j-3 ].y; 
			
			mesh->rbf->C[i][j] = phi ( v1, v2 );
		}
	}
	
	// Now populate matrix A:
	// First column: 1
	for (int i = 0, j = 0; i < nM ; i++)
		mesh->rbf->A[i][j] = 1.0;
	
	// Second column = xMoving
	for (int i = 0, j = 1; i < nM ; i++ )
		mesh->rbf->A[i][j] = mesh->rbf->movingNodes[i].x;
		
	// Third column = yMoving
	for (int i = 0, j = 2; i < nM ; i++ )
		mesh->rbf->A[i][j] = mesh->rbf->movingNodes[i].y;
		
	// Fourth and on: phi ( movingNode, fixedNode)
	for (int i = 0; i < nM ; i++ )
	{
		for (int j = 3; j < 3 + nP; j++)
		{
			double v1[2], v2[2];
			v1[0] = mesh->rbf->movingNodes[i].x; 
			v1[1] = mesh->rbf->movingNodes[i].y;
			
			v2[0] = mesh->rbf->fixedNodesBefore[ j-3 ].x; 
			v2[1] = mesh->rbf->fixedNodesBefore[ j-3 ].y ;
	
			mesh->rbf->A[i][j] = phi ( v1, v2 );
		}
	}
	
	// Create var for macro ind2d(i,j) := i * tam + j
	int tam = 3 + nP;
	
	// Allocate cInverse as 1-D array ( or trip )
	double * cInverseTrip = calloc( tam * tam , sizeof(double) );
	
	for( int i = 0; i < tam ; i++ )
		for (int j = 0; j < tam; j++ )
			cInverseTrip[ ind2d(i,j) ] = mesh->rbf->C[i][j];
		

	// Invert matrix - use Lapack call. *note: matrix inversion is usefull for RBF. 
    inverse( cInverseTrip, tam );
	
	for( int i = 0; i < tam ; i++ )
		for (int j = 0; j < tam; j++ )
			mesh->rbf->cInverse[i][j] = cInverseTrip[ ind2d(i,j) ];
		
	
	/* Calculate H = A * cInverse */
	
	// Allocate memory for matrix multiplication:
	double * H = calloc( nM * (3 + nP ) , sizeof(double) );
	double * A = calloc( nM * (3 + nP ) , sizeof(double) );
	
	// Allocate A as 1-D array:
		for( int i = 0; i < nM ; i++ )
			for (int j = 0; j < 3 + nP; j++ )
				A[ ind2d(i,j) ] = mesh->rbf->A[i][j];
			
	// BLAS-Lapack call for matrix multiplication on the form: C = alpha * A * B + beta * C
	// A, B and C shall be 1-D arrays for vectorization.
	

	//double alpha= 1.0, beta= 0.0;
    double alpha = 1.0, beta = 0.0; 
    char no= 'N', tr= 'T';		// Arrays are in 'transpose mode', as opposed to 'N', normal mode.
    
    (void) no;	// cast to void to silence warning (var may be usefull some other time)
    
    int m = nM; 		// m = nRows of A
    int n = tam;		// n = nColumns of B
    int k = tam;   	    // k = nColumns of A or nRows of B
    
    int lda = k;	// m if 'N' normal mode, k otherwise ('transpose mode')
    int ldb = n;	// k if 'N' normal mode, n otherwise ('transpose mode')
    int ldc = m;	// should be m

	

	dgemm_(&tr, &tr, &m, &n, &k, &alpha, A , &lda, cInverseTrip , &ldb, &beta, H ,&ldc);

	// Put them back into data-structure.
	tam = nM;

	for( int i = 0; i < nM ; i++ )
		for (int j = 0; j < 3 + nP; j++ )
			mesh->rbf->H[i][j] = H[ ind2d(j,i) ];
			
	//matrixMult( nM, 3 + nP, mesh->rbf->A, mesh->rbf->cInverse, mesh->rbf->H ); // Multiplication 'by hand'

	free( cInverseTrip );
	free( H );
	free( A );

	return;
}

// Different radial basis functions:
double phi( double * v1, double * v2 )
{
	
	double finalPhi = 0.0;
	double dx = v1[0] - v2[0];
	double dy = v1[1] - v2[1];
	double dz = 0.0;
	double r = sqrt( dx * dx + dy * dy + dz * dz ); 
	
	
	switch( RBF_TYPE )
	{
		
		case 0: // no RBF at all
	    finalPhi = 0.0;	
		
		case 1:
		// Linear bi-harmonic
		finalPhi = r;
		break;
		
		case 2: 
		// C2 Wendland function:
		if ( r >= 1.0 )
			finalPhi = 0.0;
		else
			finalPhi = ( 4.0 * r + 1.0) * pow( (1.0 - r), 4.0 );
		break;
		
		case 4:
		// C4 Wendland function:
		if ( r >= 1.0 )
			finalPhi = 0.0;
		else
			finalPhi = ( (35.0/3.0) * r * r + 6.0 * r + 1.0) * pow( (1.0 - r), 6.0 );
		break;
		
		
		case 6:
		// C6 Wendland function:
		if ( r >= 1.0 )
			finalPhi = 0.0;
		else
			finalPhi = ( 32.0 * r * r * r + 25.0 * r * r  + 8.0 * r + 1.0) * pow( (1.0 - r), 8.0 );
		break;
		
		case 9:
		// Quadric Biharmonic:
		finalPhi = ( 1.0 + r * r );
		break;
		
		case 10:
		// Inverse Quadric Biharmonic:
		finalPhi = 1.0 / ( 1.0 + r * r );
		break;
		
		case 11:
		finalPhi = 1.0 * exp( - r  ) ;
		break;
		
		case 12:
		finalPhi = r * r * log( r + 0.0000001 );
		break;
		
		
	}
	
	
	return finalPhi;
}

void matrixMult(int nRows, int nColumns, double** matA, double** matB, double** matC) 
{
  int i, j, k;
  // Major-row matrix multiplication (optimized for C)
  // However, dgemm from BLAS is faster :/
   for (i=0; i<nRows; i++)
    for (j=0; j<nColumns; j++)
		matC[i][j] = 0.0;
  
  
   for (i=0; i<nRows; i++)
    for (k=0; k<nColumns; k++)
      for (j=0; j<nColumns; j++)
		matC[i][j] = matC[i][j] + matA[i][k]*matB[k][j];


}    

void applyDisplacement( MeshStruct * mesh )
{
	
	// Generate U_p, the vector with the prescribed displacements.
	int nP = mesh->rbf->nFixedNodes;
	int nM = mesh->rbf->nMovingNodes;
	
	double * Up_x = calloc ( (3 + nP) , sizeof(double) );
	double * Up_y = calloc ( (3 + nP) , sizeof(double) );
	
	double * Um_x = calloc ( nM , sizeof(double) );
	double * Um_y = calloc ( nM , sizeof(double) );
	
	double xp, xpB, yp, ypB;
	
	for (int nodeCounter = 3; nodeCounter < 3 + nP; nodeCounter++ )
	{
		// xFixed:
		xp = mesh->rbf->fixedNodes[ nodeCounter - 3].x;
		// xFixedBefore:
		xpB = mesh->rbf->fixedNodesBefore[ nodeCounter - 3].x;
		
		// yFixed:
		yp = mesh->rbf->fixedNodes[ nodeCounter - 3].y;
		// yFixedBefore:
		ypB = mesh->rbf->fixedNodesBefore[ nodeCounter - 3].y;
		
		// displacement vector:
		Up_x[ nodeCounter ] = xp - xpB;
		Up_y[ nodeCounter ] = yp - ypB;

	}
	
	// Generate U_m, the vector with the new displacements U_m = H * U_p .
	for (int i = 0; i < nM; i++)
	{
		for( int j = 0; j < 3 + nP; j++ )
		{
			Um_x[ i ] = Um_x[ i ] + mesh->rbf->H[i][j] * Up_x[ j ];
			Um_y[ i ] = Um_y[ i ] + mesh->rbf->H[i][j] * Up_y[ j ];
			
		}
	}
	
	// Apply displacement for the moving nodes. The new nodes shall be their previous position + displacement.
	for (int movingNodeId = 0; movingNodeId < nM ; movingNodeId++)
	{
		mesh->rbf->movingNodes[movingNodeId].x = mesh->rbf->movingNodes[movingNodeId].x + Um_x[ movingNodeId ];
		mesh->rbf->movingNodes[movingNodeId].y = mesh->rbf->movingNodes[movingNodeId].y + Um_y[ movingNodeId ];
	}
	
	// Put these new positions back again at data structure mesh->vertices at their correct nodeNumber.
	for (int movingNodeId = 0, nodeNumber; movingNodeId < nM; movingNodeId++ )
	{
		nodeNumber = mesh->rbf->movingNodes[ movingNodeId ].nodeNumber ;
		mesh->vertices[ nodeNumber ][0] = mesh->rbf->movingNodes[ movingNodeId ].x;
		mesh->vertices[ nodeNumber ][1] = mesh->rbf->movingNodes[ movingNodeId ].y;
	}
	
	
	free( Up_x ); free( Up_y );
	free( Um_x ); free( Um_y );
	
	return;
}

void writeMesh ( MeshStruct * mesh )
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
	
	FILE * fw = fopen("MESH_TEST.gmsh","w");
    FILE * fort2 = fopen("fort.2", "w");
	FILE * fort3 = fopen("fort.3", "w");
	FILE * fort5 = fopen("fort.5", "w");
	
    
    
    // List of nodes:
	fprintf(fw, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n", nNodes);
	fprintf(fort2, "%d\n\n", nNodes);
	
	for (int i = 0; i < nNodes; i++ )
	{
		nodeNumber = i;
		x = mesh->vertices[nodeNumber][0];
		y = mesh->vertices[nodeNumber][1];
		z = 0.0;
		
		fprintf(fw, "%d %.15lf %.15lf %.15lf\n", nodeNumber + 1, x, y, z);
        fprintf(fort2, "%d %lf %lf\n", nodeNumber+1, x, y);
	}
	
    
    // Conectivity Table:
	int elemNumber, elemType, elemBType, nTags, bcId , elementaryTag, nNodesPerCell, beg, end;
	
	fprintf(fw, "$EndNodes\n$Elements\n%d\n", nElems);
	fprintf(fort3, "%d\n\n", nElems - mesh->sizes.bfaces );
    fprintf(fort5, "%d\n\n", mesh->sizes.bfaces );

    
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
				fprintf(fort3, "%d", elemNumber+1);
				

				// Now, loop throughout the faces of this cell
				for (int nodeId = 0; nodeId < nNodesPerCell; nodeId++ )
				{
                 	fprintf(fw, " %d", mesh->mElems[ elemNumber ].nodes[ nodeId ] + 1 );
					fprintf(fort3, " %d", mesh->mElems[ elemNumber ].nodes[ nodeId ] + 1 );
				}
				// skip line
				fprintf(fw, "\n");
				fprintf(fort3, "\n");
				
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
				fprintf(fort5, "%d %d", faceNumber, bcId); 
                
				// For each of these b-faces, run throughout their (3) nodes.
				for (int nodeId = 0; nodeId < 2 + (MESH_ORDER-1); nodeId++)
				{	
                    fprintf(fw, " %d", mesh->mBFaces[ bFaceNumber ].nodes[ nodeId ] + 1 );
                    fprintf(fort5, " %d", mesh->mBFaces[ bFaceNumber ].nodes[ nodeId ] + 1 );
				}	
				// Skip line:
				fprintf(fw, "\n");
                fprintf(fort5, "\n");
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
    fclose( fort2 );
    fclose( fort3 );
    fclose( fort5 );
    
	
	return;
}

