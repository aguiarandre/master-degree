#ifndef SPLINE_H
#define SPLINE_H

#include "native.h"

typedef struct 
{
    int bSplineId;
    int nPoints;
    int k;
    int m;
    double * knots;
    double * weights;
    
    double * x_ctrl;
    double * y_ctrl;
    double * z_ctrl;
    
    double x;
    double y;
    double z;
    
    // Member Functions here.
    
} BSpline;

// Function Prototype:

double bsplineBasis( double u, int paramsId, int splineId, BSpline * spline, int m );
void reconstructSpline( BSpline * spline, int splineId, double parameter );
bool getParameter( double x, double y, BSpline * spline, int splineId, double * t  );
void splineInit( BSpline * spline );

// HO Features
void curveMesh( MeshStruct * mesh, BSpline * spline );
void initializeCurvedVars( MeshStruct * mesh );

void reOrderBoundaries( MeshStruct * mesh );

int getWallNodes( BSpline * spline, MeshStruct * mesh, int nodeNumber );
int projectNodes( double t1, double t2, int * paramFound, BSpline * spline, int nodeId, MeshStruct * mesh );

int getInnerNodes( MeshStruct * mesh, int nodeNumber );
int getBoundNodes( MeshStruct * mesh, int nodeNumber );
int getCellNodes ( MeshStruct * mesh, int nodeNumber );

void getConnectivity( MeshStruct * mesh );

#endif 
