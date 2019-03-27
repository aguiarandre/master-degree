#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "native.h"

void hdf5_write_scalar(hid_t file, const char* name, const char* type,
                       const void* data)
{
  hid_t dset_id, dataspace;
  hsize_t dimsf_set[] = {1};
  dataspace = H5Screate_simple(1, dimsf_set, NULL);

  if ( (strcmp(type,"double") == 0) ) {
    dset_id = H5Dcreate(file, name, H5T_NATIVE_DOUBLE, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  else if ( (strcmp(type,"int") == 0) ) {
    dset_id = H5Dcreate(file, name, H5T_NATIVE_INT, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  else {
    printf("Error: hdf5_write_scalar(), unknown data type.\n");
    exit(1);
  }

  H5Sclose(dataspace);
  H5Dclose(dset_id);
}


void hdf5_write_data(hid_t file, const char* name, const char* type,
                     uint* sizes, bool attr_size_bytes, const void* data)
{
  hid_t dset_id, dataspace;

  // Create the dataspace for the dataset.
  hsize_t dimsf[2];
  dimsf[0] = sizes[0];
  dimsf[1] = sizes[1];
  dataspace = H5Screate_simple(2, dimsf, NULL);

  //Create the dataset with default properties and write data
  int item_size;
  if ( (strcmp(type,"double") == 0) ) {
    dset_id = H5Dcreate(file, name, H5T_NATIVE_DOUBLE, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace, H5P_DEFAULT, data);
    // compute size of each element in dataset.
    if ( attr_size_bytes ) {
      item_size = sizeof(double)*sizes[1];
    }
    else {
      item_size = sizes[0];
    }
  }
  else if ( (strcmp(type,"int") == 0) ) {
    dset_id = H5Dcreate(file, name, H5T_NATIVE_INT, dataspace,
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, dataspace, H5P_DEFAULT, data);
    // compute size of each element in dataset.
    if ( attr_size_bytes ) {
      item_size = sizeof(int)*sizes[1];
    }
    else {
      item_size = sizes[0];
    }
  }
  else {
    printf("Error: hdf5_write_data(), unknown data type.\n");
    exit(1);
  }

  H5Sclose(dataspace);
  H5Dclose(dset_id);

  /* attach attributes to dat */

  //open existing data set
  dset_id = H5Dopen(file, name, H5P_DEFAULT);
  //create the data space for the attribute
  hsize_t dims = 1;
  dataspace = H5Screate_simple(1, &dims, NULL);

  //Create an int attribute - size
  hid_t attribute = H5Acreate(dset_id, "size", H5T_NATIVE_INT, dataspace,
      H5P_DEFAULT, H5P_DEFAULT);
  //Write the attribute data.
  H5Awrite(attribute, H5T_NATIVE_INT, &item_size);
  //Close the attribute.
  H5Aclose(attribute);

  //Create an int attribute - dimension
  attribute = H5Acreate(dset_id, "dim", H5T_NATIVE_INT, dataspace,
      H5P_DEFAULT, H5P_DEFAULT);
  //Write the attribute data.
  H5Awrite(attribute, H5T_NATIVE_INT, &sizes[1]);
  H5Aclose(attribute);
  H5Sclose(dataspace);

  //Create an string attribute - type
  dataspace= H5Screate(H5S_SCALAR);
  hid_t atype = H5Tcopy(H5T_C_S1);
  int attlen = strlen(type);
  H5Tset_size(atype, attlen);

  attribute = H5Acreate(dset_id, "type", atype, dataspace,
      H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute, atype, type);
  H5Aclose(attribute);

  //Close the dataspace.
  H5Sclose(dataspace);
  H5Dclose(dset_id);

}


void hdf5_write_vlstring(hid_t file, const char* name, uint* sizes, const void* data)
{
  hid_t dset_id, dataspace;
	hid_t filetype, memtype;

  // Create the dataspace for the dataset.
  hsize_t dimsf[2];
  dimsf[0] = sizes[0];
  dimsf[1] = sizes[1];
  dataspace = H5Screate_simple(2, dimsf, NULL);

	filetype = H5Tcopy (H5T_C_S1);
	H5Tset_size(filetype, H5T_VARIABLE);
	memtype = H5Tcopy (H5T_C_S1);
	H5Tset_size(memtype, H5T_VARIABLE);

	dset_id = H5Dcreate(file, name, filetype, dataspace,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
}


void hdf5_read_scalar(hid_t file, char const *name, char const *type, void* data)
{
  //open existing data set
  hid_t dset_id = H5Dopen(file, name, H5P_DEFAULT);

  char* buffer;
  //initialize data buffer and read data
  if(strcmp(type,"int") == 0)
  {
     buffer = (char *)malloc(sizeof(int));
     H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
     memcpy((void *)data, (void *)buffer,sizeof(int));
  }
  else if(strcmp(type,"float") == 0)
  {
     buffer = (char *)malloc(sizeof(float));
     H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
     memcpy((void *)data, (void *)buffer,sizeof(float));
  }
  else if(strcmp(type,"double") == 0)
  {
     buffer = (char *)malloc(sizeof(double));
     H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
     memcpy((void *)data, (void *)buffer,sizeof(double));
  }
  else {
    printf("Error: hdf5_read_scalar(), unknown data type.\n");
    exit(1);
  }
  
  free(buffer);

  H5Dclose(dset_id);
}


void hdf5_read_data_dim(hid_t file, char const *name, uint* dim)
{
  //open existing data set
  hid_t dset_id = H5Dopen(file, name, H5P_DEFAULT);

  //get OID of the dim attribute
  hid_t attr = H5Aopen(dset_id, "dim", H5P_DEFAULT);
  H5Aread(attr,H5T_NATIVE_INT,dim);
  H5Aclose(attr);
  H5Dclose(dset_id);
}


void hdf5_read_data(hid_t file, char const *name, char const *type, void* data)
{
  hid_t dset_id;
  
  //open existing data set
  dset_id = H5Dopen(file, name, H5P_DEFAULT);

  if(strcmp(type,"int") == 0)
  {
     H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  else if(strcmp(type,"float") == 0)
  {
     H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  else if(strcmp(type,"double") == 0)
  {
     H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  }
  else {
    printf("Error: hdf5_read_data(), unknown data type.\n");
    exit(1);
  }

  H5Dclose(dset_id);
}


void hdf5_read_vlstring(hid_t file, const char* name, uint* sizes, void* data)
{
	hid_t dset_id, space, memtype;

	dset_id = H5Dopen(file, name, H5P_DEFAULT);

	space = H5Dget_space (dset_id);

	// Create the memory datatype.
	memtype = H5Tcopy (H5T_C_S1);
	H5Tset_size (memtype, H5T_VARIABLE);

	// read the data
	H5Dread (dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

	H5Dclose (dset_id);
	H5Sclose (space);
	H5Tclose (memtype);
}
