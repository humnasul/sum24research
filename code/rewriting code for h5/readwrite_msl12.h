////////////////////////////////////////////////////////////////////////////////////////////////////////
// This code is a part of NOAA STAR Ocean Color MSL12 destriping package.                             //
// Developed by Karlis Mikelsons at NOAA/NESDIS/STAR and GST, INC.                 12/12/2014         //
// Please see the README file for credits and scknowledgments.                                        //
////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef READWRITE_MSL12
#define READWRITE_MSL12

#include <string.h>
#include "readwrite_msl12_hdf4.h"
#include "readwrite_msl12_hdf5.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class dtype>
int read_msl12(dtype ** buffer, int * nx, int * ny, char * dataset_name, char * filename, int readwrite = 0) {

   int ilen = strlen(filename);
   if( (ilen>=4) &&
       (filename[ilen-4]=='.') &&
       (filename[ilen-3]=='h') &&
       (filename[ilen-2]=='d') &&
       (filename[ilen-1]=='f') )
   {
       return read_msl12_hdf4(buffer, nx, ny, dataset_name, filename, readwrite);
   }
   else if( (ilen>=3) &&
            (filename[ilen-3]=='.') &&
            (filename[ilen-2]=='n') &&
            (filename[ilen-1]=='c') )
   {
       return read_msl12_hdf5(buffer, nx, ny, dataset_name, filename, readwrite);
   }
   else if( (ilen>=3) &&
            (filename[ilen-3]=='.') &&
            (filename[ilen-2]=='h') &&
            (filename[ilen-1]=='5') )
   {
       return read_msl12_hdf5(buffer, nx, ny, dataset_name, filename, readwrite);
   }
   else return -1; 

};


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class dtype>
int       write_msl12(dtype ** buffer, int * nx, int * ny, char * dataset_name, char * filename) {
    return read_msl12(         buffer,       nx,       ny,        dataset_name,        filename, /* readwrite = */ 1 );
};



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
int write_msl12_attr_strarr(char * filename, char * attr_name, int nstr, char ** str_arr){

   int ilen = strlen(filename);
   if( (ilen>=4) &&
       (filename[ilen-4]=='.') &&
       (filename[ilen-3]=='h') &&
       (filename[ilen-2]=='d') &&
       (filename[ilen-1]=='f') )
   {
       return write_msl12_hdf4_attr_strarr(filename, attr_name, nstr, str_arr);
   }

   else if( (ilen>=3) &&
            (filename[ilen-3]=='.') &&
            (filename[ilen-2]=='n') &&
            (filename[ilen-1]=='c') )
   {
       return write_msl12_hdf5_attr_strarr(filename, attr_name, nstr, str_arr);
   }
   else return -1; 

};


#endif
