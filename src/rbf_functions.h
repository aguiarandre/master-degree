#ifndef RBF_FUNCTIONS_H
#define RBF_FUNCTIONS_H

#include "native.h"
#include <stdlib.h> // malloc
#include <string.h>
#include <stddef.h>
#include <math.h>

// functions to propagate deformation of high order projection:

void initializeData( MeshStruct * mesh );
void propagateDeformation( MeshStruct * mesh );
void getFixedNodes( MeshStruct * mesh );
void getMovingNodes( MeshStruct * mesh );
void populateMatrices( MeshStruct * mesh );
void applyDisplacement( MeshStruct * mesh );
double phi( double * v1, double * v2);
void matrixMult(int nRows, int nColumns, double** matA, double** matB, double** matC);
void matrixMult2(int nRows, int nColumns, double* matA, double* matB, double* matC) ;
void writeMesh( MeshStruct * mesh );
void destroyData( MeshStruct * mesh );


// LU decomoposition of matrix A, in 'array' form.
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

// generate inverse of a matrix given its LU decomposition
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

 void dgemm_(char* TRANSA, char* TRANSB, const int* M,
               const int* N, const int* K, double* alpha, double* A,
               const int* LDA, double* B, const int* LDB, double* beta,
               double* C, const int* LDC);


/*
void cblas_dgemm (const CBLAS_LAYOUT layout, const CBLAS_TRANSPOSE TransA, 
				  const CBLAS_TRANSPOSE TransB, const int M, const int N, 
				  const int K, const double alpha, const double *A, const int lda, 
				  const double *B, const int ldb, const double beta, double *C, const int ldc);
  */             



#endif 
