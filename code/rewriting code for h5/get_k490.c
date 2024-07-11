/*---------------------------------------------------------------------*/
/* get_k_490.c -  diffuse attenuation coefficient for MSl12.           */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     k_490 - diffuse attenuation coefficient, 1 value per pixel.     */
/*                                                                     */
/* Algorithm provided by: J. Mueller                                   */
/* OCTS/POLDER coefficents: S. Bailey, 16 July 2001                    */
/* Written by: B. Franz, SAIC-GSC, SIMBIOS Project, October 1999       */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
// #include "l12_proto.h"



/* Added on 5/30/2013 by lidej */
void get_k490_noaa_modis_pix(float nLw[], float Fo[], float *k490, float *kpar)
{
    const int idx488 = 3, idx551 = 5, idx555 = 6, idx645 = 7, idx667 = 8;
    const float a = 0.1853, b = -1.349, Kw = 0;
    const float badval = -1, maxval = 6.4;
    float kd_clear, kd_turbid, kpar_clear, kpar_turbid, w;
    float Rrs488, Rrs645, Rrs667, rrs488, rrs645, rrs667;

    if (nLw[idx488] <= 0.0 || nLw[idx551] <= 0.0 ) {
      *k490 = badval;
      *kpar = badval;
    }
    else {
        kd_clear = badval;
        kd_turbid = badval;
        Rrs488 = nLw[idx488]/Fo[idx488];
        Rrs645 = nLw[idx645]/Fo[idx645];
        Rrs667 = nLw[idx667]/Fo[idx667];
        rrs488 = 4*Rrs488/(0.52+1.7*Rrs488);
        rrs645 = 4*Rrs645/(0.52+1.7*Rrs645);
        rrs667 = 4*Rrs667/(0.52+1.7*Rrs667);

        if (nLw[idx551] < 20.0 )
            kd_clear = a*pow(nLw[idx488]/nLw[idx551],b) + Kw;
        else 
            kd_clear = a*pow(nLw[idx488]/nLw[idx555],b) + Kw;

        kd_clear = MAX(MIN(kd_clear, maxval), 0);
        kpar_clear = MAX(0,0.0864+0.8*kd_clear-0.00137/kd_clear);

        /* use 667 band first */
        if (nLw[idx667] < 10.0 && nLw[idx667] > 0.0)
            kd_turbid = 2.697e-4/rrs488+1.045*rrs667/rrs488+4.18*(7e-4+2.7135*rrs667)* (1-0.52*exp(-2.533e-3/rrs488-9.817*rrs667/rrs488));
        else 
            kd_turbid = -9.785e-4/rrs488+0.8321*rrs645/rrs488+4.18*(-2.54e-3 +2.1598*rrs645)*(1-0.52*exp(9.19e-3/rrs488-7.81*rrs645/rrs488));

        kd_turbid = MAX(MIN(kd_turbid, maxval), 0);
        kpar_turbid = 0.8045*pow(kd_turbid,0.917);

        /* calculate weigh w */
        w = -1.175 + 4.512*rrs645/rrs488;
        if (Rrs645/Rrs488 < 0.26064) w = 0.0;
        if (Rrs645/Rrs488 > 0.4821)  w = 1.0;
        *k490 = (1.0 - w)*kd_clear   + w*kd_turbid;
        *kpar = (1.0 - w)*kpar_clear + w*kpar_turbid;
    }

    if ( (nLw[idx488] > 20.0) || (nLw[idx555] > 20.0) || (*k490 < 0.0) ) {
        *k490 = badval;
        *kpar = badval;
    }

    return;
};


/* Added on 5/30/2013 by lidej */
void get_k490_noaa_modis_pix_(float nLw[], float Fo[], float *k490, float *kpar)
{
    get_k490_noaa_modis_pix(nLw, Fo, k490, kpar);
    return;
};


/* Added on 5/30/2013 by lidej */
void get_k490_noaa_viirs_pix(float nLw[], float Fo[], float *k490, float *kpar)
{
    const int idx486 = 2, idx551 = 3, idx671 = 4;
    const float a = 0.1853, b = -1.349, Kw = 0;
    const float badval = -1, maxval = 6.4;
    float  kd_clear, kd_turbid, kpar_clear, kpar_turbid, w;
    float  rrs486, rrs671, Rrs486, Rrs671;

    if ( (nLw[idx486] <= 0.0) || (nLw[idx551] <= 0.0) ) {
        *k490 = badval;
        *kpar = badval;
    }
    else {
        kd_clear = badval;
        kd_turbid = badval;
        Rrs486 = nLw[idx486]/Fo[idx486];
        Rrs671 = nLw[idx671]/Fo[idx671];
        rrs486 = 4*Rrs486/(0.52+1.7*Rrs486);
        rrs671 = 4*Rrs671/(0.52+1.7*Rrs671);

        kd_clear = a*pow(nLw[idx486]/nLw[idx551],b) + Kw;
        kd_clear = MAX(MIN(kd_clear,maxval),0);
        kpar_clear = MAX(0,0.0864+0.8*kd_clear-0.00137/kd_clear);
        kd_turbid = 2.697e-4/rrs486
                  + 1.045*rrs671/rrs486
                  + 4.18*(7e-4 + 2.7135*rrs671) * (1-0.52*exp(-2.533e-3/rrs486-9.817*rrs671/rrs486));

        kd_turbid = MAX(MIN(kd_turbid,maxval),0);
        kpar_turbid = 0.8045*pow(kd_turbid,0.917);

        /* calculate weigh w */
        w = -1.175 + 4.512*rrs671/rrs486;
        if (Rrs671/Rrs486 < 0.26064) w = 0.0;
        if (Rrs671/Rrs486 > 0.4821)  w = 1.0;
        *k490 = (1.0 - w)*kd_clear   + w*kd_turbid;
        *kpar = (1.0 - w)*kpar_clear + w*kpar_turbid;
    }

    if ( (nLw[idx486] > 20.0) || (nLw[idx551] > 20.0) || (*k490 < 0.0) ) {
        *k490 = badval;
        *kpar = badval;
    }

    return;
};


/* Added on 5/30/2013 by lidej */
void get_k490_noaa_viirs_pix_(float nLw[], float Fo[], float *k490, float *kpar)
{
    get_k490_noaa_viirs_pix(nLw, Fo, k490, kpar);
    return;
};


