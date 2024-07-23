////////////////////////////////////////////////////////////////////////////////////////////////////////
// This code is a part of NOAA STAR Ocean Color MSL12 destriping package.                             //
// Developed by Karlis Mikelsons at NOAA/NESDIS/STAR and GST, INC.                 12/12/2014         //
// Please see the README file for credits and scknowledgments.                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef READWRITE_MSL12_HDF5
#define READWRITE_MSL12_HDF5
//defining the file

#include <stdio.h>
#include <string.h>
#include "hdf5.h"
//including packages to use string and hdf5
#define MAX_STR_LEN 1024
// defining max string length for this code



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// templates - allow you to write generic code that can work with different data types.

template <class dtype>
// dtype is a placeholder for a data type that will be specified when the template is instantiated
hid_t get_hdf5_data_type();
// declares a function named get_hdf5_data_type, return type is hid_t

template <>
hid_t get_hdf5_data_type<char>()   { return H5T_NATIVE_CHAR;   }
// specialization of template for char
//returns hid_t type
// when you call get_hdf5_data_type<char>(), it will return the H5T_NATIVE_CHAR value
// "H5T_NATIVE_CHAR is a predefined constant in HDF5 (Hierarchical Data Format) 
//      that represents the native char data type for the platform"

template <>
hid_t get_hdf5_data_type<short>()  { return H5T_NATIVE_SHORT;  }
// specialization of template for short
//returns hid_t type
// when you call get_hdf5_data_type<short>(), it will return the H5T_NATIVE_SHORT value
// "H5T_NATIVE_SHORT is a predefined constant in HDF5 (Hierarchical Data Format) 
//      that represents the native short data type for the platform"

template <>
hid_t get_hdf5_data_type<int>()    { return H5T_NATIVE_INT;    }
// specialization of template for int

template <>
hid_t get_hdf5_data_type<long>()   { return H5T_NATIVE_LONG;   }
// specialization of template for long

template <>
hid_t get_hdf5_data_type<float>()  { return H5T_NATIVE_FLOAT;  }
// specialization of template for float

template <>
hid_t get_hdf5_data_type<double>() { return H5T_NATIVE_DOUBLE; }
// specialization of template for double


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class dtype>
int read_msl12_hdf5(dtype ** buffer, int * nx, int * ny, char * dataset_name, char * filename, int readwrite = 0) {

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // this finction reads or writes (parameter readwrite = 0 or 1, respectively)
    // a two-dimensional array of data 
    // identified by string sds_name in hdf4 file
    // from/to hdf file with name given by string filename
    //
    // on return, nx and ny are the dimensions of the two-dimensional array read/written
    // with ny being the major index and nx the minor index
    // so that entry (ix,iy) corresponds to location buffer[0][iy*nx+ix]
    // for (0 <= ix < nx) and (0 <= iy < ny)
    //
    // in case of read:
    //   this function allocates buffer where data are read, and points buffer[0] to this buffer;
    //   is it up to caller to deallocate this buffer later
    //
    // in case of write:
    //   buffer[0] points to data that should be written when this function is called;
    //   the buffer of data remains intact
    //
    // on success, the return value is 0;
    // nonzero value means problem of some kind
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    hid_t   file_id, dataset, dataset_factors, dataspace, H5T_DATA_TYPE;
    herr_t  hdferr;
    int     info, rank_dset, i, iprint = 0;
    unsigned long long   dims[2], dimsizes[2], maxdimsizes[2];
    char    dataset_full_name[1024];

    // look for datasets in these 4 locations
    char    loc1[] = "/AOD/";
    char    loc2[] = "/Lrc/";
    char    loc3[] = "/Lt/";
    char    loc4[] = "/Rrs/";
    char    *dataset_locations[4] = { loc1, loc2, loc3, loc4 };
    int     idataset_locations, ndataset_locations = 4;


    if(iprint > 0) printf("read_msl12_hdf5\n");     

    H5T_DATA_TYPE = get_hdf5_data_type<dtype>();
    if(iprint>0) printf("H5T_DATA_TYPE = %i\n", H5T_DATA_TYPE); 

    // open instance of HDF5 library for use
    hdferr = H5open();
    if(hdferr!=0) { 
        printf("Cannot initialize HDF5 library!\n"); 
        return -1; 
    }

    // open file for read/write
    if(readwrite==0)  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    else              file_id = H5Fopen(filename, H5F_ACC_RDWR,   H5P_DEFAULT);
    if(file_id<0) { 
        printf("Cannot open HDF5 file %s!\n", filename); 
        return -1; 
    }
    dims[1] = 2;
    dims[2] = 0;

    // try to locate the dataset
    htri_t ispathvalid = 0;
    for(idataset_locations = 0; idataset_locations<ndataset_locations; idataset_locations++){
        sprintf(dataset_full_name, "%s%s\0", dataset_locations[idataset_locations], dataset_name);
//        if(iprint>0) printf("Trying %s\n", dataset_full_name);

        // see if we can find dataset_full_name
        ispathvalid = H5Lexists(file_id, dataset_full_name, H5P_DEFAULT);
        if(ispathvalid==1) {
//            if(iprint>0) printf("Path valid, %i\n", ispathvalid);
            break;
        }
//        else {
//            if(iprint>0) printf("Path %s not valid %i\n", dataset_full_name, ispathvalid);
//        }
    }

    // if dataset nowhere found, return
    if(ispathvalid!=1) {

        // close file
        hdferr = H5Fclose(file_id);
        if(hdferr<0) {
            printf("Cannot close HDF5 file %s!\n", filename);
            return -1;
        }

        // close hdf5 library instance
        hdferr = H5close();
        if(hdferr<0) {
            printf("Cannot close HDF5 library!\n");
            return -1;
        }
        
        return 1;
    }


    // try to open dataset dataset_full_name
    dataset = H5Dopen(file_id, dataset_full_name, H5P_DEFAULT);
    if(dataset<0) { 
        printf("Cannot open HDF5 dataset %s!\n", dataset_full_name); 

        // close file
        hdferr = H5Fclose(file_id);
        if(hdferr<0) {
            printf("Cannot close HDF5 file %s!\n", filename);
            return -1;
        }

        // close hdf5 library instance
        hdferr = H5close();
        if(hdferr<0) {
            printf("Cannot close HDF5 library!\n");
            return -1;
        }
        
        return -1; 
    }

    // get the dataspace for the dataset
    dataspace = H5Dget_space(dataset);
    /* H5Dget_space() makes a copy of the dataspace of the dataset specified by dset_id. 
    The function returns an identifier for the new copy of the dataspace.
    */
    if(dataspace<0) {
         printf("Cannot open HDF5 dataspace for dataset %s!\n", dataset_full_name); 
         return -1; 
    }

    // H5Sget_simple_extent_ndims - Determines the dimensionality of a dataspace.
    // get the number of dimensions and make sure it is 2
    rank_dset = H5Sget_simple_extent_ndims(dataspace);
    if(rank_dset!=2) { 
        printf("Unexpected rank of dataspace %i expected 2\n", rank_dset); 
        return -1; 
    }

    // get the size of each dimentsion
    rank_dset = H5Sget_simple_extent_dims(dataspace, dimsizes, maxdimsizes);
    if(rank_dset<0)  {
        printf("Cannot get dataspace dimensions!\n");
        return -1; 
    }
    if(iprint>0) printf("dimsizes = %i  %i\n", dimsizes[0], dimsizes[1]);

    *ny = dimsizes[0];
    *nx = dimsizes[1];
    if(iprint>0) printf("nx = %i ny = %i\n", *nx, *ny);
    

    // check if we need to read or write data
    if(readwrite==0) {

        // allocate data buffer for data to be read into
        *buffer = (dtype *) malloc(dimsizes[0]*dimsizes[1]*sizeof(dtype));
        if(*buffer==NULL) {
            printf("Cannot allocate memory\n"); 
            return -1; 
        }

        // read data into buffer
        // H5Dread - reads a dataset, specified by its identifier dset_id, from the file into an application memory buffer *buffer
        hdferr = H5Dread( dataset, H5T_DATA_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
        if(hdferr<0) { 
            printf("Cannot read data to hdf!\n");
            return -1; 
        }
    } 
    else {

        // write data from buffer into hdf5 file
        hdferr = H5Dwrite(dataset, H5T_DATA_TYPE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *buffer);
        if(hdferr<0) { 
            printf("Cannot write data to hdf!\n");
            return -1; 
        }
    } //  if(readwrite==0)

    // close everything:
    // close dataset
    hdferr = H5Dclose(dataset);
    if(hdferr<0) { 
        printf("Cannot close HDF5 dataset %s!\n", dataset_full_name);
        return -1; 
    }

    // close file
    hdferr = H5Fclose(file_id);
    if(hdferr<0) {
        printf("Cannot close HDF5 file %s!\n", filename);
        return -1;
    }

    // close hdf5 library instance
    hdferr = H5close();
    if(hdferr<0) {
        printf("Cannot close HDF5 library!\n");
        return -1;
    }

    return 0;

};



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class dtype>
int       write_msl12_hdf5(dtype ** buffer, int * nx, int * ny, char * dataset_name, char * filename) {
    return read_msl12_hdf5(         buffer,       nx,       ny,        dataset_name,        filename, /* readwrite = */ 1 );
};



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
int write_msl12_hdf5_attr_strarr(char * filename, char * attr_name, int nstr, char ** str_arr){

    hid_t   file_id, dataset, attr_id, type_id;
    herr_t  hdferr;
    int     i, iprint = 1;   /* set iprint > 0 for debug printout */
    int     retval = 0;

    char attrbuff[MAX_STR_LEN], attrbuff_old[MAX_STR_LEN];
    for (i=0; i<MAX_STR_LEN; i++) { attrbuff[i] = 0; attrbuff_old[i] = 0; }

    char attrFieldStr[] = "/\0";

    if(nstr<=0) return 0; // nothing to do

    hdferr = H5open();
    if(hdferr!=0) {
        printf("Cannot initialize HDF5 library!\n");
        return -1;
    }

    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    if(file_id<0) { 
        printf("Cannot open HDF5 file %s!\n", filename);
        return -1;
    }

    if(H5Aexists(file_id, attr_name)>0) {

        // attribute already exists
        if(iprint>0) printf("Attribute %s     found in file %s\n", attr_name, filename);

        // open attribute 
        attr_id = H5Aopen(file_id, attr_name, H5P_DEFAULT );
        if(attr_id<0) { printf("Cannot open a HDF5 attribute\n"); return -1; }

        // get attribute type
        type_id = H5Aget_type(attr_id);
        if(type_id<0) { printf("Cannot get a HDF5 attribute type (H5Aget_type)\n"); return -1; }        

        // read attribute 
        hdferr = H5Aread( attr_id, type_id, attrbuff_old );
        if(hdferr<0) {

            // attribute found, but cannot read it - problem
            printf("Cannot read attr %s with H5Aread\n", attr_name);  
            return -1; 

        }
        else {

            // print warning message for attributes already destriped
            printf("WARNING: Following products were already destriped: %s = %s\n", attr_name, attrbuff_old); 

        }

        retval = 1;

        hdferr = H5Adelete( file_id, attr_name ); 
        if(hdferr<0) {

            // cannot delete the old attribute to write a new one
            printf("Cannot delete old HDF5 attribute (H5Adelete)!\n");
            return -1;

        }

    } // if attribute already exists ...

 
    // modify attribute for destriped products
    // concatenate strings in array str_arr
    int istr, nchar = 0;
    attrbuff[0] = '\0';
    for(istr=0; istr<nstr; istr++){
        if(istr>0) strcat(attrbuff," \0");
        strcat(attrbuff, str_arr[istr]);
    }
    nchar = strlen(attrbuff);
    strcat(attrbuff, "\0");
    if(iprint>0) { printf("%s\n", attrbuff); fflush(stdout); }

    // create a new attribute type based on length of the string
    type_id = H5Tcopy(H5T_C_S1);
    if(type_id<0)  {
        printf("Cannot create a new datatype!\n");
        return -1;
    }

    // adjust the size of datatype
    hdferr = H5Tset_size(type_id, nchar+1 /* 256 H5T_VARIABLE */);
    if(hdferr<0) {
        printf("Cannot set HDF5 type size (H5Tset_size)!\n");
         return -1;
    }

    // create a new dataspace
    hid_t space_id = H5Screate(H5S_SCALAR);
    if(space_id<0) {
        printf("Cannot create a new dataspace!\n");
        return -1;
    }

    // attribute does not exist, create it
    attr_id = H5Acreate(file_id, attr_name, type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
    if(attr_id<0) {
        printf("Cannot create a HDF5 attribute\n");
        return -1;
    }

    // write the new attribute to the hdf file
    hdferr = H5Awrite( attr_id, type_id, (const void *) attrbuff );
    if(hdferr==-1) {
        printf("Cannot write attribute!\n");
        return -1;
    }

    // close everything
    // close the attribute
    hdferr = H5Aclose(attr_id);
    if(hdferr<0) {
        printf("Cannot close HDF5 attribute %s!\n", attr_name);
        return -1;
    }

    // close dataspace
    hdferr = H5Sclose(space_id);
    if(hdferr<0) {
        printf("Cannot close HDF5 dataspace (H5Sclose) %s!\n", attr_name);
        return -1;
    }

    // close datatype
    hdferr = H5Tclose(type_id);
    if(hdferr<0) {
        printf("Cannot close HDF5 datatype (H5Tclose) %s!\n", attr_name);
        return -1;
    }

    // close the file
    hdferr = H5Fclose(file_id);
    if(hdferr<0) {
        printf("Cannot close HDF5 file %s!\n", filename);
        return -1;
    }

    // close library use
    hdferr = H5close();
    if(hdferr<0) {
        printf("Cannot close HDF5 library!\n");
        return -1;
    }

    return retval;

};



#endif
