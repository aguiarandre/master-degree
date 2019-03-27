#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "native.h"
#include "solution.h"
#include "io_cgns.h"
#include "io_bru.h"
#include "BSpline.h"
#include "igesRead.h"
#include "rbf_functions.h"
#include "moveMesh.h"

/* 
 * argv[0] program file name,
 * argv[1] operation mode,
 */
 
   /* program operation mode. Valid values are:
   * 0 - convert CGNS mesh to native format, save to disk (hdf5)
   *     argv[2] input data file name.
   *     argv[3] output data file name.
   * 1 - convert native (hdf5) mesh and solution to CGNS format
   *     argv[2] input mesh file name.
   *     argv[3] solution data file name.
   *     argv[4] output data file name.
   * 2 - convert CGNS to BRU3D format
   *     argv[2] input data file name.
   *     argv[3] output data file name.
   * 3 - convert linear mesh to curved mesh.
   *     argv[2] input mesh file name (.cgns/hdf5)
   *     argv[3] input CAD file name (.iges)
   */
 

int main(int argc, char **argv)
{
	if ( argc < 3 )
	{ 
		printf("Missing Arguments.\n");
		help_datproc();
		exit(1);
	}

  int op_mode = (*argv[1] - '0') % 48;


  MeshStruct mesh;
  mesh_init(&mesh); 
  
  // Create RBF struct to store values.
  rbfStruct rbfPropag;
  mesh.rbf = &rbfPropag;

 // Initialize BSpline
  BSpline * bSpline = (BSpline*) malloc( sizeof(BSpline) *  MAXSPLINES );
  splineInit( bSpline );
	

 
  if ( op_mode == 0 ) {

    io_cgns_read(argv[2], &mesh);

    mesh_create_faces(&mesh);

	mesh_set_bcs(&mesh); // 			mesh.zones[#].bc_id = 0, 1 or 2.
	
	mesh_report(&mesh);




	mesh_write_meta_info(argv[3], &mesh);
	
    mesh_write_data_native(argv[3], &mesh);


  }
  else if ( op_mode == 1 ) {
  
	SolutionStruct sol;
	
	mesh_read_meta_info(argv[2], &mesh);
	
	mesh_read_data_native(argv[2], &mesh);
	
	solution_read(argv[3], &mesh, &sol);
	
	io_cgns_write(argv[4], &mesh, &sol);
	
  }
  else if ( op_mode == 2 ) {

    /* input file creation and BC definition, possibly by GUI. */

    io_cgns_read(argv[2], &mesh);

    //mesh_create_faces(&mesh);

    mesh_report(&mesh);

    io_bru_write(&mesh);

  }
  else if ( op_mode == 3 ) {

    ////////////// process the solver hdf5 and convert it to CGNS. 
    
    // read mesh file from argv[2]
    // read CAD file from argv[3]
    // call library to curve mesh - use as inputs argv[2], argv[3] and argv[4]
    // call function to write final mesh - use as input argv[5]     
    
	io_cgns_read(argv[2], &mesh);

	mesh_create_faces(&mesh);

	mesh_set_bcs(&mesh); // 			mesh.zones[#].bc_id = 0, 1 or 2.

	mesh_report(&mesh);
	
	iges_read(argc, argv[3], &mesh, bSpline);
	
	curveMesh( &mesh, bSpline );    

	propagateDeformation( &mesh );
    
    moveMesh( &mesh );
    
    /* WRITE RESULTING MESH */
	writeMesh( &mesh );
	
    
  }

  else {
    printf("Error: unknown op_mode parameter.\n");
    help_datproc();
    exit(1);
  }

  /* free mesh struct and memory. */
  mesh_finish(&mesh);
}
