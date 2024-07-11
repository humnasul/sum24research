#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define HMAX 128000

int write_plotscript(int ns, int nx1, int ny1, float * chlor, short * kd490, char * basename, char ** bandnames, float * Qplotmin, float * Qplotmax)
{

    int histogram[HMAX];
    int i, j, n, ix, iy, is;
    float Qplotmin_chlor, Qplotmax_chlor, Qplotmin_kd490, Qplotmax_kd490;


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all pixels and find the histogram for values of chlor to get the plotting range
    n = 0;
    for(i=0;i<HMAX;i++){ histogram[i] = 0; }
    for(iy=0; iy<nx1*ny1; iy++) {
        i = (int) round(100.0*chlor[iy]);
        if((i>=0)&&(i<HMAX)) {
            histogram[i]++;
            n++;
        }
    }
    j = 0;
    for(i=0;i<HMAX;i++){
        j += histogram[i];
        if(j>0.01*n){ break; }
    }
    Qplotmin_chlor = 0.01*i;
    j = 0;
    for(i=0;i<HMAX;i++){
        j += histogram[i];
        if(j>0.99*n){ break; }
    }
    Qplotmax_chlor = 0.01*i;
    if(Qplotmax_chlor-Qplotmin_chlor<0.05) { Qplotmax_chlor = Qplotmin_chlor + 0.05; }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loop over all pixels and find the histogram for values of chlor to get the plotting range
    n = 0;
    for(i=0;i<HMAX;i++){ histogram[i] = 0; }
    for(iy=0; iy<nx1*ny1; iy++) {
        i = (int) kd490[iy];
        if((i>=0)&&(i<HMAX)) {
            histogram[i]++;
            n++;
        }
    }
    j = 0;
    for(i=0;i<HMAX;i++){
        j += histogram[i];
        if(j>0.01*n){ break; }
    }
    Qplotmin_kd490 = i*2.0E-4;
    j = 0;
    for(i=0;i<HMAX;i++){
        j += histogram[i];
        if(j>0.99*n){ break; }
    }
    Qplotmax_kd490 = i*2.0E-4;
    if(Qplotmax_kd490-Qplotmin_kd490<0.05) { Qplotmax_kd490 = Qplotmin_kd490 + 0.05; }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////
    // write a seadas script file for plotting all images
    /////////////////////////////////////////////////////////////////////

    // create a file name for plot script file
    char plotscriptfile[256];
    sprintf(plotscriptfile, "%s.plot.txt.destr.txt\0", basename);

    // create and open plotscript file
    FILE * fp = NULL;
    fp = fopen(plotscriptfile,"w");
    if(fp==NULL) {
        printf("Cannot write plot script file \n");
        return -1;
    }
  
    fprintf(fp,"l2file=\'%s\'\n\n", basename); 
    fprintf(fp,"l2file1=l2file"); 
    fprintf(fp,"+\'.destr.hdf\'\n"); 
 
    fprintf(fp,"\nload, l2file1, prod_name=[\'l2_flags\',\'chlor_a\',\'k490_noaa\'"); 
    for(is=0; is<ns; is++){
        fprintf(fp,",\'%s\'", bandnames[is]); 
    }
    fprintf(fp,"]\n\n");

    fprintf(fp,"xmin=%i\n",1);
    fprintf(fp,"xmax=%i\n",nx1);
    fprintf(fp,"ymin=%i\n",1);
    fprintf(fp,"ymax=%i\n",ny1);

    // higlint - overlay higlint flag over chlor to show where glint is hign
    fprintf(fp,"loadpal,\'$OCDAPS/config/color_luts/standard/02-standard_chl.lut\'\n");
    fprintf(fp,"mband_cmd, cmd_array=[\'result=b2\',\'iy=ARRAY_BIT(r1,[4],jcnt)\',\'if (jcnt GT 0) then result(iy) = 640.0\'], bandname=\'higlint\'\n");
    fprintf(fp,"display, fbuf=%i, band_no=%i, smin=0.01, smax=640.0, stype=\'LOG\'\n", ns+4, ns+4);
    fprintf(fp,"ofile = l2file + \'_higlint.png.destr.png\'\n");
    fprintf(fp,"out, ofile, /display, /cbar, region=[xmin,xmax,ymin,ymax]\n\n");
 
    // chlor
    fprintf(fp,"loadpal,\'$OCDAPS/config/color_luts/standard/02-standard_chl.lut\'\n");
    fprintf(fp,"display, fbuf=2, band_no=2, smin=%f, smax=%f, stype=\'LOG\'\n", Qplotmin_chlor, Qplotmax_chlor);
    fprintf(fp,"ofile = l2file + \'_chlor.png.destr.png\'\n");
    fprintf(fp,"out, ofile, /display, /cbar, region=[xmin,xmax,ymin,ymax]\n\n");

    // kd490
    fprintf(fp,"loadpal,\'$OCDAPS/config/color_luts/standard/05-standard_k490.lut\'\n");
    fprintf(fp,"display, fbuf=3, band_no=3, smin=%f, smax=%f, stype=\'LOG\'\n", Qplotmin_kd490, Qplotmax_kd490);
    fprintf(fp,"ofile = l2file + \'_k490.png.destr.png\'\n");
    fprintf(fp,"out, ofile, /display, /cbar, region=[xmin,xmax,ymin,ymax]\n\n");

    // nLw
    for(is=0; is<ns; is++){
        // destriped bands (nLw)
        fprintf(fp,"loadpal,\'$OCDAPS/config/color_luts/standard/03-standard_sst.lut\'\n");
        fprintf(fp,"display, fbuf=%i, band_no=%i, smin=%f, smax=%f, stype=\'LINEAR\'\n", is+4, is+4, Qplotmin[is], Qplotmax[is]);
        fprintf(fp,"ofile = l2file + \'_%s.png.destr.png\'\n", bandnames[is]);
        fprintf(fp,"out, ofile, /display, /cbar, region=[xmin,xmax,ymin,ymax]\n\n");
    }

    fprintf(fp,"clear_up\n");
    fclose(fp);
    /////////////////////////////////////////////////////////////////////

    return 0;
};

