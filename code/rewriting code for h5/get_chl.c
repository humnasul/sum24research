/*
Humna Sultan
Edited to suit Rrs and OC-SMART input
removed all instances of nLw / Fo and replaced with Rrs
- in the case that code could not be replaced, it was removed
- OC-SMART file only uses Rrs and nLw input does not exist
*/

#include <stdlib.h>
#include <math.h>
//#include "l12_proto.h"

#ifndef MAX
#define MAX(A,B)    ((A) > (B) ? (A) : (B))  /* Greater of (A,B) */
#endif

#ifndef MIN
#define MIN(A,B)    ((A) < (B) ? (A) : (B))  /* Lesser  of (A,B) */
#endif

#define CHL_BAD   -1.0

static float pi = 3.141592654;
static float chlmin = 0.0;
static float chlmax = 100.0;
static float chlbad = CHL_BAD;

/////////////////////////////////////////////////////////////////////////////////////////
// calculate chl_a for modis
/////////////////////////////////////////////////////////////////////////////////////////
float get_chl_oc3_modis(float Rrs[])
{
    float rat;
    float minRrs;
    float chl = chlbad;
    //float Rrs551_555;
    float Rrs443 = Rrs[1];
    float Rrs488 = Rrs[3];
    // Rrs551 and Rrs555 variables removed

    //    static float a[] = {0.283,-2.753,1.457,0.659,-1.403};
    static float a[] = {0.2424,-2.5828,1.7057,-0.3415,-0.8818};	// updated v6 coefficient - lidej.2013.10.25

    minRrs = MIN(Rrs443,Rrs488);

    /* For cases where  Rrs551 is saturated, Rrs555 is used 3/7/07*/ 
    // code for Rrs551 and Rrs555 removed

    // We require Rrs551 to be positive, and we require that if any band
    // goes negative, it must occur in order of wavelength 

    if ( /* (Rrs551_555 > 0.0) && */ Rrs488 > 0.0 ) {
        rat = MAX(Rrs443,Rrs488); // removed dividing by Rrs551_555

        /* Fail if ratio is unphysical (Rat=0.21 -> Chl=640) */
        if (rat > 0.21) { 
            rat = log10(rat);
            chl = (float) pow(10.0,(a[0]+rat*(a[1]+rat*(a[2]+rat*(a[3]+rat*a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
        }
    }

    /*
    if ( (nLw[3] > 20.0) || (nLw[6] > 20.0) ) {
        chl = -1;
    } */
   // removed for now

    /*set the value to be -2 for negative water leaving radiance*/ 
    if ( /*  (Rrs551_555 < 0.0) || */ MAX(Rrs443,Rrs488) < 0.0 ) {
        chl=-1;
    }
    
    return chl;
};

/////////////////////////////////////////////////////////////////////////////////////////
// calculate chl_a for viirs
//  REMOVED
/////////////////////////////////////////////////////////////////////////////////////////

