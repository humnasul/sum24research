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
    float Rrs551_555;
    float Rrs443 = Rrs[1];
    float Rrs488 = Rrs[3];
    float Rrs551 = Rrs[5];
    float Rrs555 = Rrs[6];

    //    static float a[] = {0.283,-2.753,1.457,0.659,-1.403};
    static float a[] = {0.2424,-2.5828,1.7057,-0.3415,-0.8818};	// updated v6 coefficient - lidej.2013.10.25

    minRrs = MIN(Rrs443,Rrs488);

    /* For cases where  Rrs551 is saturated, Rrs555 is used 3/7/07*/ 
    if (nLw[5] < 20.0) {
        Rrs551_555 = Rrs551;
    }
    else {
        Rrs551_555 = Rrs555;
    }

    // We require Rrs551 to be positive, and we require that if any band
    // goes negative, it must occur in order of wavelength 
    if ( (Rrs551_555 > 0.0) && (Rrs488 > 0.0) ) {
        rat = MAX(Rrs443,Rrs488)/Rrs551_555;

        /* Fail if ratio is unphysical (Rat=0.21 -> Chl=640) */
        if (rat > 0.21) { 
            rat = log10(rat);
            chl = (float) pow(10.0,(a[0]+rat*(a[1]+rat*(a[2]+rat*(a[3]+rat*a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
        }
    }

    if ( (nLw[3] > 20.0) || (nLw[6] > 20.0) ) {
        chl = -1;
    }

    /*set the value to be -2 for negative water leaving radiance*/ 
    if ( (Rrs551_555 < 0.0) || (MAX(Rrs443,Rrs488) < 0.0) ) {
        chl=-1;
    }
    
    return chl;
};

/////////////////////////////////////////////////////////////////////////////////////////
// calculate chl_a for viirs
/////////////////////////////////////////////////////////////////////////////////////////
float get_chl_oc3_viirs(float nLw[], float Fo[])
{

    float rat;
    float chl = chlbad;

    float Rrs445 = nLw[1]/Fo[1];
    float Rrs488 = nLw[2]/Fo[2];
    float Rrs555 = nLw[3]/Fo[3];

    static float a[] = {0.2228,-2.4683,1.5867,-0.4275,-0.7768};

    // We require Rrs555 to be positive, and we require that if any band
    //   goes negative, it must occur in order of wavelength 
    if ( (Rrs555 > 0.0) && (Rrs488 > 0.0) ) {
        rat = MAX(Rrs445,Rrs488)/Rrs555;

        /* Fail if ratio is unphysical (Rat=0.21 -> Chl=640) */
        if (rat > 0.21) {
            rat = log10(rat);
            chl = (float) pow(10.0,(a[0]+rat*(a[1]+rat*(a[2]+rat*(a[3]+rat*a[4])))));
            chl = (chl > chlmin ? chl : chlmin);
        }
    }

    if ( (nLw[2] > 20.0) || (nLw[3] > 20.0) ) {
        chl=CHL_BAD;
    }

    /*set the value to be -2 for negative water leaving radiance*/
    if ( (Rrs555 < 0.0) || (MAX(Rrs445,Rrs488) < 0.0) ) {
        chl=CHL_BAD;
    }

    return (chl);
};

