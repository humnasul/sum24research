///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//   This package [destripe_msl12] is designed as a standalone postprocessing tool                       //
//   for MSL12 Ocean Color satellite data processing. It performs destriping of the                      //
//   water leaving radiances calculated by MSL12 and recalculates the derived                            //
//   quantities: chlorophyll-a concentration (chlor-a) and diffuse light attenuation                     //
//   coefficient at 490nm (kd490).                                                                       //
//                                                                                                       //
//   This code was developed in NOAA/NESDIS/STAR Ocean Color group by Karlis Mikelsons,                  //
//   and is partially based on an earlier version of destriping algorithm for use with the               //
//   Sea Surface Temperature data written by Marouan Bouali.                                             //
//                                                                                                       //
//   Please acknowledge any use of these codes in your presentations                                     //
//   and publications by citing the following publications:                                              //
//                                                                                                       //
//   1. K. Mikelsons, M. Wang, L. Jiang, and M. Bouali,                                                  //
//      "Destriping algorithm for improved satellite-derived                                             //
//      ocean color product imagery," Opt. Express  22, 28058-28070 (2014).                              //
//                                                                                                       //
//   2. M. Bouali and A. Ignatov, "Adaptive Reduction of Striping                                        //
//      for Improved Sea Surface Temperature Imagery                                                     //
//      from Suomi National Polar-Orbiting Partnership (S-NPP)                                           //
//      Visible Infrared Imaging Radiometer Suite (VIIRS)",                                              //
//      J. Atmos. Oceanic Technol., 31, 150–163 (2014).                                                  //
//                                                                                                       //
//   Please report any bugs to Karlis.Mikelsons@noaa.gov                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mfhdf.h"
#include "readwrite_msl12.h"
#include "l2_flags.h"
#include "get_chl.c"
#include "get_k490.c"
#include "destripe.h"
#include "fill_restore_viirs_bowtie.c"
#include "get_nif_ndf.c"
#include "write_plotscript.c"

#define HMAX 128000
#define NTHREADS 8

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {

    int *l2_flags1 = NULL;
    short *buffer0 = NULL;
    short *buffer1 = NULL;

    float **bufferf1 = NULL;
    float **bufferf2 = NULL;
    int   **binary_M = NULL;
    float **bowtie   = NULL;
    int histogram[HMAX];

    int nx0, ny0, nx1, ny1, nx2, ny2, nya, nyb;
    int is, ns, ix, iy, status1, i, j, n, nthreads;

    int Ndet, Niter;
    int Ndet_arr[40], Niter_arr[40];
    float Qmin, Qmax, Thresh_x, Thresh_y, NEdQ, r1, scalefact = 1000.0;
    float Qmin_arr[40], Qmax_arr[40], Tx_arr[40], Ty_arr[40], NEdQ_arr[40];
    float Qplotmin_chlor, Qplotmax_chlor, Qplotmin_kd490, Qplotmax_kd490, Qplotmin[40], Qplotmax[40];

    int is_nLw_445 = 0, is_nLw_488 = 0, is_nLw_555 = 0, is_nLw_672 = 0;
    int is_nLw_443 = 0, is_nLw_551 = 0, is_nLw_645 = 0, is_nLw_667 = 0;

    float Fo[16];
    float **nLws, *chlor;
    short *kd490, *kdpar;
    float kd490f, kdparf;
    int  nattrstr = 0;
    char ** attrstr_arr = NULL, **bandnames = NULL;

    clock_t start_ticks, end_ticks;
    start_ticks = clock();
    time_t  start_time, end_time;
    start_time = time(NULL);

    attrstr_arr = allocate_2d_c(64,64);
    bandnames   = allocate_2d_c(64,64);
    if( (attrstr_arr==NULL) || (bandnames==NULL) ) { 
        printf("ERROR: Cannot allocate memory\n");
        return -1; 
    }

    // check if enough command line arguments
    if(argc<3) { 
        printf("Not enough arguments!\nUsage:\n");
        printf(" %s msl12_produced_L2.hdf destriping_param_file.txt [prev=msl12_prev_L2.hdf] [next=msl12_next_L2.hdf] [-plots]\n", argv[0]);
        printf("where:\n");
        printf("msl12_produced_L2.hdf      is a file containing the granule to be destriped\n");
        printf("destriping_param_file.txt  is destriping parameter file (differs for VIIRS and MODIS)\n");
        printf("prev=msl12_prev_L2.hdf     denotes optional preceding granule to be used for better destriping near lower boundary\n");
        printf("next=msl12_next_L2.hdf     denotes optional following granule to be used for better destriping near upper boundary\n");
        printf("-plots                     if present, generate a script to plot destriped quantities with OCDAPS\n");
        printf("Note: granules denoted by prev and next options need to contain original (NOT destriped) msl12 data\n");
        return 0;
    }

    printf("Destriping %s\n", argv[1]);
    
    char * pprevfile = NULL;
    char * pnextfile = NULL;
    bool plotscript = false;

    // process command line arguments
    for(i=3; i<argc; i++) {
        if( (argv[i][0]=='p') && 
            (argv[i][1]=='r') &&
            (argv[i][2]=='e') &&
            (argv[i][3]=='v') &&
            (argv[i][4]=='=') )
        {
            pprevfile=&(argv[i][5]);
        }

        if( (argv[i][0]=='n') && 
            (argv[i][1]=='e') &&
            (argv[i][2]=='x') &&
            (argv[i][3]=='t') &&
            (argv[i][4]=='=') ) 
        {
            pnextfile=&(argv[i][5]);
        }

        if( (argv[i][0]=='-') && 
            (argv[i][1]=='p') &&
            (argv[i][2]=='l') &&
            (argv[i][3]=='o') &&
            (argv[i][4]=='t') &&
            (argv[i][5]=='s') ) 
        {
            plotscript = true;
        }

    } // for(i=3; i<argc; i++) 

    if(pprevfile!=NULL) printf("Using prevfile = %s\n", pprevfile);
    if(pnextfile!=NULL) printf("Using nextfile = %s\n", pnextfile);

    // set number of threads
    nthreads = omp_get_max_threads();
    if(nthreads>8) nthreads = 8; // limit number of threads to 8 - more than that is very unlikely to improve performance
    omp_set_num_threads(nthreads);
    // printf("maxthreads =  %i\n", omp_get_max_threads());


    // find if this is VIIRS or MODIS file
    // search for the first letter of filename
    int pos = 0;
    for(i=(strlen(argv[1])-2); i>=0; i--){
        if(argv[1][i]=='/') {
            pos = i+1;
            break;
        }
    }

    // set up Fo data needed for chlor and kd490
    if(argv[1][pos] == 'V') { // VIIRS
        printf("Using VIIRS bands for chlor and kd490\n");
        Fo[1] = 190.26;
        Fo[2] = 198.87;
        Fo[3] = 184.23;
        Fo[4] = 150.45;
    } 
    else { // MODIS
        printf("Using MODIS bands for chlor and kd490\n");
        Fo[1] = 187.30;
        Fo[3] = 194.85;
        Fo[5] = 186.58;
        Fo[6] = 183.95;
        Fo[7] = 158.59;
        Fo[8] = 152.20;
    }

    ///////////////////////////////////////////////////////////////////////
    //////////////////// read parameters  /////////////////////////////////
    
    // open parameter file
    FILE *fp = fopen(argv[2],"r");
    if(fp==NULL) {  printf("ERROR: Cannot open parameter file\n");  return -1; }
 
    // read parameters
    is = 0;
    for(i=0; i<40; i++){

        j = fscanf(fp,"%s %i %i %f %f %f %f %f\n", bandnames[is], 
                   &(Ndet_arr[is]), &(Niter_arr[is]), &(NEdQ_arr[is]),
                   &(Tx_arr[is]), &(Ty_arr[is]), &(Qmin_arr[is]), &(Qmax_arr[is]));

        if(j==8) { is++; } else break;
    }
    ns = is;

    // close param. file
    fclose(fp);
    ////////////////// done reading parameters  /////////////////////////// 
    ///////////////////////////////////////////////////////////////////////


    printf("Reading l2_flags\n");
    status1 = read_msl12(&l2_flags1, &nx1, &ny1, "l2_flags\0", argv[1]);
    if(status1!=0) { printf("ERROR: Cannot read flag data! %i\n", 10*status1);  return -1; }
    printf("Data read\n");

    // allocate space for nLw data  - later needed  for kd490 and chlor
    // as well as temp. 2D arrays needed for destriping 
    // use dimensions from flags
    nLws     = allocate_2d_f(ny1*nx1, 16);
    // also, allocate other data buffers
    bufferf1 = allocate_2d_f(ny1+160, nx1);
    bufferf2 = allocate_2d_f(ny1+160, nx1);
    binary_M = allocate_2d_i(ny1+160, nx1);
    bowtie   = allocate_2d_f(ny1+160, nx1);
    buffer1  = (short*) malloc((ny1+160)*nx1*sizeof(short));
    if( (nLws==NULL) || (bufferf1==NULL) || (bufferf2==NULL) || (binary_M==NULL) || (bowtie==NULL) || (buffer1==NULL) ) {
        printf("ERROR: Cannot allocate data \n"); 
        return -1;
    }
    printf("Memory allocated\n");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///   loop over all possible bands              /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(is=0; is<ns; is++) {
        printf("======================================================================================\n");
        printf("Variable %s\n", bandnames[is]);
 
        NEdQ     = NEdQ_arr[is];
        Thresh_x = Tx_arr[is];
        Thresh_y = Ty_arr[is];
        Qmin     = Qmin_arr[is];
        Qmax     = Qmax_arr[is];

        Ndet     = Ndet_arr[is];
        Niter    = Niter_arr[is];

        printf("%s %5i %5i %f %f %f %f %f\n", bandnames[is], (Ndet_arr[is]), (Niter_arr[is]), (NEdQ_arr[is]),
                                               (Tx_arr[is]), (Ty_arr[is]), (Qmin_arr[is]), (Qmax_arr[is]));

        // read the data for previous adjacent granule, if supplied (for better treatment of boundary)
        nya = 0;
        if(pprevfile!=NULL) {
            printf("Using previous granule %s for better destriping near lower boundary\n", pprevfile);
            status1 = read_msl12(&buffer0, &nx0, &ny0, bandnames[is], pprevfile);
            if(status1!=0) { 
                printf("Band %s not found in file %s\n", bandnames[is], pprevfile); 
            }
            else {
                nya = 80;
                if(ny0<80) nya = ny0;
#pragma omp parallel for 
                for(iy = 0; iy<nx0*nya; iy++) buffer1[iy] = buffer0[nx0*(ny0 - nya) + iy];

                free(buffer0);
            }
        }
        ny0 = nya; // index of where the lower boundary rows end, and where the current granule data will start


        // read the data of granule to destripe
        status1 = read_msl12(&buffer0, &nx1, &ny1, bandnames[is], argv[1]);
        if(status1!=0) {
            printf("Band %s not found in file %s\n", bandnames[is], argv[1]);
            continue; // if this band not found, move on to next band
        }
#pragma omp parallel for 
        for(iy=0; iy<nx1*ny1; iy++) buffer1[nx0*ny0 + iy] = buffer0[iy];
        free(buffer0);


        // read the data for next adjacent granule, if supplied (for better treatment of boundary)
        nyb = 0;
        if(pnextfile!=NULL) {
            printf("Using   next   granule %s for better destriping near upper boundary\n", pnextfile);
            status1 = read_msl12(&buffer0, &nx2, &ny2, bandnames[is], pnextfile);
            if(status1!=0) { 
                printf("Band %s not found in file %s\n", bandnames[is], pnextfile); 
            }
            else {
                nyb = 80;
                if(ny2<80) nyb = ny2;
#pragma omp parallel for 
                for(iy = 0; iy<nx2*nyb; iy++) buffer1[nx1*(nya+ny1) + iy] = buffer0[iy];
                free(buffer0);
            }
        }
        ny2 = ny0 + ny1 + nyb; // total number of rows to destripe

        printf("Data read\n");


        // move -1000 to -32766
#pragma omp parallel for 
        for(iy=0; iy<nx1*ny2; iy++) { if(buffer1[iy]==-1000) buffer1[iy] = -32766; }

        // scale floating point data to physical values
#pragma omp parallel for 
        for(iy=0; iy<nx1*ny2; iy++) { bufferf1[0][iy] = buffer1[iy]*0.001; }

        // fill VIIRS bowtie area with interpolated values
        if(argv[1][pos] == 'V') { // VIIRS
            fill_viirs_bowtie( bufferf1, bowtie, nx1, ny2, -1.0, 20.0, -100.0);
        }

        /////////////////////////////////////////////////////////////////////////////////////
        // loop over all pixels and find the histogram for x derivatives
        // using integer array makes histograms easier
        n = 0;
        for(i=0;i<HMAX;i++){ histogram[i] = 0; }
        for(iy=0; iy<ny1; iy++) {
            for(ix=0; ix<(nx1-1); ix++){

                if( (buffer1[(ny0+iy)*nx1 + ix    ] < Qmin*scalefact) || 
                    (buffer1[(ny0+iy)*nx1 + ix    ] > Qmax*scalefact) || 
                    (l2_flags1[iy*nx1 + ix    ] & HIGLINT) )  continue;

                if( (buffer1[(ny0+iy)*nx1 + ix + 1] < Qmin*scalefact) ||
                    (buffer1[(ny0+iy)*nx1 + ix + 1] > Qmax*scalefact) ||
                    (l2_flags1[iy*nx1 + ix + 1] & HIGLINT) )  continue;

                i =  abs( buffer1[(ny0+iy)*nx1 + ix ] - buffer1[(ny0+iy)*nx1 + ix + 1] );
                if(i<HMAX) {
                    histogram[i]++;
                    n++;
                }
            } // ix
        } // iy

        j = 0;
        // for(i=0;i<HMAX;i++){ printf("%i %i\n", i,  histogram[i]); }
        for(i=0;i<HMAX;i++){ j += histogram[i]; if(j>0.99*n) break; }
        r1 =  i * 1.2 / scalefact;
        if( Thresh_x > r1 )  Thresh_x = r1;


        /////////////////////////////////////////////////////////////////////////////////////
        // loop over all pixels and find the histogram for y derivatives
        // using integer array makes histograms easier
        n = 0;
        for(i=0;i<HMAX;i++){ histogram[i] = 0; }
        for(iy=0; iy<(ny1-1); iy++) {
            for(ix=0; ix<nx1; ix++){

                if( (buffer1[(ny0+iy  )*nx1 + ix] < Qmin*scalefact) ||
                    (buffer1[(ny0+iy  )*nx1 + ix] > Qmax*scalefact) || 
                    (l2_flags1[iy*nx1 + ix    ] & HIGLINT) )  continue;

                if( (buffer1[(ny0+iy+1)*nx1 + ix] < Qmin*scalefact) ||
                    (buffer1[(ny0+iy+1)*nx1 + ix] > Qmax*scalefact) || 
                    (l2_flags1[(iy+1)*nx1 + ix] & HIGLINT) )  continue;

                i =  abs( buffer1[(ny0+iy)*nx1 + ix] - buffer1[(ny0+iy+1)*nx1 + ix] );
                if(i<HMAX) {
                    histogram[i]++;
                    n++;
                } 
            } // ix
        } // iy

        j = 0;
        // for(i=0;i<HMAX;i++){ printf("%i %i\n", i,  histogram[i]); }
        for(i=0;i<HMAX;i++){ j += histogram[i]; if(j>0.99*n) break; }
        r1 =  i * 1.2 / scalefact;
        if( Thresh_y > r1 )  Thresh_y = r1;

        /////////////////////////////////////////////////////////////////////////////////////
        // loop over all pixels and find the histogram for values to get the plotting range
        n = 0;
        for(i=0;i<HMAX;i++){ histogram[i] = 0; }
        for(iy=0; iy<ny1; iy++) {
            for(ix=0; ix<nx1; ix++){
                if( (buffer1[(ny0+iy)*nx1 + ix]<Qmin*scalefact) || 
                    (buffer1[(ny0+iy)*nx1 + ix]>Qmax*scalefact) 
               /* ||(l2_flags1[iy*nx1 + ix    ] & HIGLINT) */ )   continue;

                i =  buffer1[(ny0+iy)*nx1 + ix ] + 1000;
                if((i>=0)&&(i<HMAX)) {
                    histogram[i]++;
                    n++;
                }
            } // ix
        } // iy 

        j = 0;
        for(i=0;i<HMAX;i++){ j += histogram[i]; if(j>0.01*n) break; }
        Qplotmin[is] =  (i-1000)/scalefact;

        j = 0;
        for(i=0;i<HMAX;i++){ j += histogram[i]; if(j>0.99*n) break; }
        Qplotmax[is] = (i-1000)/scalefact;

        if(Qplotmax[is]-Qplotmin[is]<0.02) { Qplotmax[is] = Qplotmin[is] + 0.02; }
        /////////////////////////////////////////////////////////////////////////////////////
    
        printf("Thresh_x = %f Thresh_y = %f\n", Thresh_x, Thresh_y);


        // we use negative NEdT value to trigger automatic estimation of NEdT by destriping algorithm
        status1 = destripe_main_frame(bufferf1, bufferf2, binary_M, Ndet, Niter, -NEdQ, Thresh_x, Thresh_y, Qmin, Qmax, nx1, ny2);
        if(status1!=0) { printf("ERROR: Destriping failure %i!\n", 100*status1); }
        printf("Destriping done\n");

        // calculate NIF and NDF
        float nif, ndf;
        int nwf;
        get_nif_ndf( &(bufferf1[ny0]), &(bufferf2[ny0]), &(binary_M[ny0]), nx1, ny1, &nif, &ndf, &nwf);
        //    printf("NIF = %2.6f\n", nif);
        //    printf("NDF = %2.6f\n", ndf);
        //    printf("NWF = %i\n", nwf);
        printf("  Band    NIF        NDF      domain_size\n %s  %2.6f  %2.6f   %i\n", bandnames[is], nif, ndf, nwf);

        // restore bowtie areas in case of VIIRS
        if(argv[1][pos] == 'V') { // VIIRS
            restore_viirs_bowtie( bufferf2, bowtie, nx1, ny2, -100.0);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(argv[1][pos] == 'V') { // VIIRS
            if(strcmp(bandnames[is],"nLw_445\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][1] = bufferf2[ny0][i]; printf("nLw_445 saved\n"); is_nLw_445 = 1; }
            if(strcmp(bandnames[is],"nLw_488\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][2] = bufferf2[ny0][i]; printf("nLw_488 saved\n"); is_nLw_488 = 1; }
            if(strcmp(bandnames[is],"nLw_555\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][3] = bufferf2[ny0][i]; printf("nLw_555 saved\n"); is_nLw_555 = 1; }
            if(strcmp(bandnames[is],"nLw_672\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][4] = bufferf2[ny0][i]; printf("nLw_672 saved\n"); is_nLw_672 = 1; }
        }
        else {  // modis
            if(strcmp(bandnames[is],"nLw_443\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][1] = bufferf2[ny0][i]; printf("nLw_443 saved\n"); is_nLw_443 = 1; }
            if(strcmp(bandnames[is],"nLw_488\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][3] = bufferf2[ny0][i]; printf("nLw_488 saved\n"); is_nLw_488 = 1; }
            if(strcmp(bandnames[is],"nLw_551\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][5] = bufferf2[ny0][i]; printf("nLw_551 saved\n"); is_nLw_551 = 1; }
            if(strcmp(bandnames[is],"nLw_555\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][6] = bufferf2[ny0][i]; printf("nLw_555 saved\n"); is_nLw_555 = 1; }
            if(strcmp(bandnames[is],"nLw_645\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][7] = bufferf2[ny0][i]; printf("nLw_645 saved\n"); is_nLw_645 = 1; }
            if(strcmp(bandnames[is],"nLw_667\0")==0){ for(i=0; i<nx1*ny1; i++)  nLws[i][8] = bufferf2[ny0][i]; printf("nLw_667 saved\n"); is_nLw_667 = 1; }
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // scale data back to short integer values
#pragma omp parallel for 
        for(iy=0; iy<nx1*ny1; iy++) { buffer1[iy] = (short) round(bufferf2[ny0][iy]*scalefact); }

        status1 = write_msl12(&buffer1, &nx1, &ny1, bandnames[is], argv[1]);
        if(status1!=0) { printf("ERROR: Cannot write nLw data %i\n", 10*status1); }
        sprintf(attrstr_arr[nattrstr], "%s\0",  bandnames[is]); nattrstr++; // add this band to the list of bands destriped

        printf("======================================================================================\n");
    } // for(is=0; is<ns; is++)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///   loop over all possible bands complete     /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    free(buffer1);
    free(l2_flags1);
    free(bufferf1[0]); free(bufferf1);
    free(bufferf2[0]); free(bufferf2);
    free(binary_M[0]); free(binary_M);
    free(bowtie[0]);   free(bowtie);

    /////////////////////////////////////////////////////////////////////
    // calculate chl and kd490
    chlor = (float *) malloc(ny1*nx1*sizeof(float));
    kd490 = (short *) malloc(ny1*nx1*sizeof(short));
    kdpar = (short *) malloc(ny1*nx1*sizeof(short));
    if( (chlor==NULL) || (kd490==NULL) || (kdpar==NULL) ) {
        printf("ERROR: Cannot allocate memory\n");
        return -1;
    }
  
    int getchlor = 0, getkd490 = 0;

    if(argv[1][pos] == 'V') { // VIIRS
        if( is_nLw_445 && is_nLw_488 && is_nLw_555 ) {
            getchlor = 1;
            for(i=0; i<nx1*ny1; i++){
                chlor[i] = get_chl_oc3_viirs(nLws[i], Fo);
            }  
        }

        if( is_nLw_488 && is_nLw_555 && is_nLw_672 ) {
            getkd490 = 1;
            for(i=0;i<nx1*ny1;i++){
                // kd490[i] = (short) (get_k490_noaa_viirs_pix(nLws[i], Fo)/2.0E-4);
                get_k490_noaa_viirs_pix( nLws[i], Fo, &kd490f, &kdparf );
                kd490[i] = (short) (kd490f/2.0E-4);
                kdpar[i] = (short) (kdparf/2.0E-4);
            }  
        }

    }
    else { // MODIS
        if( is_nLw_443 && is_nLw_488 && is_nLw_551 && is_nLw_555 ) {
            getchlor = 1;
            for(i=0;i<nx1*ny1;i++){
                chlor[i] = get_chl_oc3_modis(nLws[i], Fo);
            }
        }
        if( is_nLw_488 && is_nLw_551 && is_nLw_555 && is_nLw_645 && is_nLw_667 ) {
            getkd490 = 1;
            for(i=0;i<nx1*ny1;i++){
                // kd490[i] = (short) (get_k490_noaa_modis_pix(nLws[i], Fo)/2.0E-4);
                get_k490_noaa_modis_pix( nLws[i], Fo, &kd490f, &kdparf );
                kd490[i] = (short) (kd490f/2.0E-4);
                kdpar[i] = (short) (kdparf/2.0E-4);
            }
        }  
    }

    // if chlor_a was updated, save it to hdf file
    if(getchlor) {

        // write calculated chlor to hdf file
        status1 = write_msl12(&chlor, &nx1, &ny1, "chlor_a\0",   argv[1]);
        if(status1!=0) {
            printf("ERROR: Cannot write chlor data %i\n", 10*status1); 
            printf("(Are chlor data present in %s?)\n", argv[1]);
        }
        else {
            // add this band to the list of bands destriped
            sprintf(attrstr_arr[nattrstr], "%s\0", "chlor_a\0"); 
            nattrstr++; 
            printf("chlor_a recalculated and saved\n");
        }

    }
    else {
        printf("One or more bands needed for chlor_a calculation was not processed\n");
        printf("chlor_a was not updated\n");
    }

    // if kd490 was updated, save it to hdf file
    if(getkd490) {

        // write calculated kd490 to hdf file
        status1 = write_msl12(&kd490, &nx1, &ny1, "k490_noaa\0", argv[1]);
        if(status1!=0) {
            printf("ERROR: Cannot write kd490 data %i\n", 10*status1); 
            printf("(Are kd490 data present in %s?)\n", argv[1]);
        }
        else {
            // add this band to the list of bands destriped
            sprintf(attrstr_arr[nattrstr], "%s\0", "k490_noaa\0");
            nattrstr++;
            printf("kd490 recalculated and saved\n");
        }

        // write calculated kdpar to hdf file
        status1 = write_msl12(&kdpar, &nx1, &ny1, "kpar_noaa\0", argv[1]);
        if(status1!=0) { 
            printf("ERROR: Cannot write kdpar data %i\n", 10*status1); 
            printf("(Are kdpar data present in %s?)\n", argv[1]);
        }
        else {
            // add this band to the list of bands destriped
            sprintf(attrstr_arr[nattrstr], "%s\0", "kpar_noaa\0");
            nattrstr++;
            printf("kdpar recalculated and saved\n");
        }

    }
    else {
        printf("One or more bands needed for kd490 and kdpar calculation was not processed\n");
        printf("kd490 and kdpar were not updated\n");
    }
    printf("======================================================================================\n");
    // calculate chl and kd490 complete
    /////////////////////////////////////////////////////////////////////

    // write destriping attribute to hdf file
    char attr_name[] = "DESTRIPED_PRODUCTS\0";
    status1 = write_msl12_attr_strarr(argv[1], attr_name, nattrstr, attrstr_arr);
    if(status1<0) { printf("ERROR: Cannot write attribute in hdf file %i\n", 40*status1); }
    else           { printf("DESTRIPED_PRODUCTS attribute written in hdf file\n"); }
    free(attrstr_arr[0]); free(attrstr_arr);
    printf("======================================================================================\n");


    /////////////////////////////////////////////////////////////////////
    // write a script file for plotting all images
    /////////////////////////////////////////////////////////////////////
    if(plotscript) {
        write_plotscript(ns, nx1, ny1, chlor, kd490, argv[1], bandnames, Qplotmin, Qplotmax);
    }

    free(chlor); 
    free(kd490);
    free(kdpar);
    free(nLws[0]); free(nLws);
    free(bandnames[0]); free(bandnames);
    /////////////////////////////////////////////////////////////////////
  
    end_ticks = clock();
    end_time = time(NULL);
    printf("threads = %i   cpu_time_used = %i (sec)   real_time_used = %i (sec)\n", nthreads,(int) ((end_ticks - start_ticks)/(1.0*CLOCKS_PER_SEC)), (int) (end_time-start_time) ) ;
    return 0;
};
