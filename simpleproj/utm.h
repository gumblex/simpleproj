#include "sp_common.h"

// More exact: Poder/Engsager
typedef struct PoderEngsager
{
    /* parameters from PJ */
    double  a;       /* semimajor axis (radius if eccentricity==0) */
    double  es;      /* first  eccentricity squared */
    double  k0;      /* General scaling factor - e.g. the 0.9996 of UTM */
    double  lam0;    /* central meridian */
    double  phi0;    /* central parallel */
    double  x0;      /* false easting */
    double  y0;      /* false northing  */

    double    Qn;     /* Merid. quad., scaled to the projection */
    double    Zb;     /* Radius vector in polar coord. systems  */
    double    cgb[6]; /* Constants for Gauss -> Geo lat */
    double    cbg[6]; /* Constants for Geo lat -> Gauss */
    double    utg[6]; /* Constants for transv. merc. -> geo */
    double    gtu[6]; /* Constants for geo -> transv. merc. */
} UTM;

PJ_XY utm_forward (PJ_LONLAT coords, UTM *P);

PJ_LONLAT utm_inverse (PJ_XY xy, UTM *P);

int utm_get_zone(double lon);

int utm_setup(UTM *P, int zone, int south);
