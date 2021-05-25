#include <math.h>
#include "utm.h"

/* Adapted from PROJ: src/projections/tmerc.cpp

 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
*/

/* Constant for "exact" transverse mercator */
#define PROJ_ETMERC_ORDER 6

inline static double adjlon (double lon) {
    /* Let lon slightly overshoot, to avoid spurious sign switching at the date line */
    if (fabs (lon) < M_PI + 1e-12)
        return lon;

    /* adjust to 0..2pi range */
    lon += M_PI;

    /* remove integral # of 'revolutions'*/
    lon -= M_TWOPI * floor(lon / M_TWOPI);

    /* adjust back to -pi..pi range */
    lon -= M_PI;

    return lon;
}


/*****************************************************************************/
//
//                  Exact Transverse Mercator functions
//
//
// The code in this file is largly based upon procedures:
//
// Written by: Knud Poder and Karsten Engsager
//
// Based on math from: R.Koenig and K.H. Weise, "Mathematische
// Grundlagen der hoeheren Geodaesie und Kartographie,
// Springer-Verlag, Berlin/Goettingen" Heidelberg, 1951.
//
// Modified and used here by permission of Reference Networks
// Division, Kort og Matrikelstyrelsen (KMS), Copenhagen, Denmark
//
/*****************************************************************************/

/* Helper functions for "exact" transverse mercator */
inline
static double gatg(const double *p1, int len_p1, double B, double cos_2B, double sin_2B) {
    double h = 0, h1, h2 = 0;

    const double two_cos_2B = 2*cos_2B;
    const double* p = p1 + len_p1;
    h1 = *--p;
    while (p - p1) {
        h = -h2 + two_cos_2B*h1 + *--p;
        h2 = h1;
        h1 = h;
    }
    return (B + h*sin_2B);
}

/* Complex Clenshaw summation */
inline
static double clenS(const double *a, int size,
                    double sin_arg_r, double cos_arg_r,
                    double sinh_arg_i, double cosh_arg_i,
                    double *R, double *I) {
    double      r, i, hr, hr1, hr2, hi, hi1, hi2;

    /* arguments */
    const double* p = a + size;
    r          =  2*cos_arg_r*cosh_arg_i;
    i          = -2*sin_arg_r*sinh_arg_i;

    /* summation loop */
    hi1 = hr1 = hi = 0;
    hr = *--p;
    for (; a - p;) {
        hr2 = hr1;
        hi2 = hi1;
        hr1 = hr;
        hi1 = hi;
        hr  = -hr2 + r*hr1 - i*hi1 + *--p;
        hi  = -hi2 + i*hr1 + r*hi1;
    }

    r   = sin_arg_r*cosh_arg_i;
    i   = cos_arg_r*sinh_arg_i;
    *R  = r*hr - i*hi;
    *I  = r*hi + i*hr;
    return *R;
}


/* Real Clenshaw summation */
static double clens(const double *a, int size, double arg_r) {
    double      r, hr, hr1, hr2, cos_arg_r;

    const double* p = a + size;
    cos_arg_r  = cos(arg_r);
    r          =  2*cos_arg_r;

    /* summation loop */
    hr1 = 0;
    hr = *--p;
    for (; a - p;) {
        hr2 = hr1;
        hr1 = hr;
        hr  = -hr2 + r*hr1 + *--p;
    }
    return sin (arg_r)*hr;
}

/* Ellipsoidal, forward */
static PJ_XY exact_e_fwd (PJ_LP lp, UTM *P) {
    PJ_XY xy = {0.0,0.0};

    /* ell. LAT, LNG -> Gaussian LAT, LNG */
    double Cn  = gatg (P->cbg, PROJ_ETMERC_ORDER, lp.phi, cos(2*lp.phi), sin(2*lp.phi));
    /* Gaussian LAT, LNG -> compl. sph. LAT */
    const double sin_Cn = sin (Cn);
    const double cos_Cn = cos (Cn);
    const double sin_Ce = sin (lp.lam);
    const double cos_Ce = cos (lp.lam);

    const double cos_Cn_cos_Ce = cos_Cn*cos_Ce;
    Cn     = atan2 (sin_Cn, cos_Cn_cos_Ce);

    const double inv_denom_tan_Ce = 1. / hypot (sin_Cn, cos_Cn_cos_Ce);
    const double tan_Ce = sin_Ce*cos_Cn * inv_denom_tan_Ce;
#if 0
    // Variant of the above: found not to be measurably faster
    const double sin_Ce_cos_Cn = sin_Ce*cos_Cn;
    const double denom = sqrt(1 - sin_Ce_cos_Cn * sin_Ce_cos_Cn);
    const double tan_Ce = sin_Ce_cos_Cn / denom;
#endif

    /* compl. sph. N, E -> ell. norm. N, E */
    double Ce  = asinh ( tan_Ce );     /* Replaces: Ce  = log(tan(FORTPI + Ce*0.5)); */

/*
 *  Non-optimized version:
 *  const double sin_arg_r  = sin(2*Cn);
 *  const double cos_arg_r  = cos(2*Cn);
 *
 *  Given:
 *      sin(2 * Cn) = 2 sin(Cn) cos(Cn)
 *          sin(atan(y)) = y / sqrt(1 + y^2)
 *          cos(atan(y)) = 1 / sqrt(1 + y^2)
 *      ==> sin(2 * Cn) = 2 tan_Cn / (1 + tan_Cn^2)
 *
 *      cos(2 * Cn) = 2cos^2(Cn) - 1
 *                  = 2 / (1 + tan_Cn^2) - 1
 */
    const double two_inv_denom_tan_Ce = 2 * inv_denom_tan_Ce;
    const double two_inv_denom_tan_Ce_square = two_inv_denom_tan_Ce * inv_denom_tan_Ce;
    const double tmp_r = cos_Cn_cos_Ce * two_inv_denom_tan_Ce_square;
    const double sin_arg_r  = sin_Cn * tmp_r;
    const double cos_arg_r  = cos_Cn_cos_Ce * tmp_r - 1;

/*
 *  Non-optimized version:
 *  const double sinh_arg_i = sinh(2*Ce);
 *  const double cosh_arg_i = cosh(2*Ce);
 *
 *  Given
 *      sinh(2 * Ce) = 2 sinh(Ce) cosh(Ce)
 *          sinh(asinh(y)) = y
 *          cosh(asinh(y)) = sqrt(1 + y^2)
 *      ==> sinh(2 * Ce) = 2 tan_Ce sqrt(1 + tan_Ce^2)
 *
 *      cosh(2 * Ce) = 2cosh^2(Ce) - 1
 *                   = 2 * (1 + tan_Ce^2) - 1
 *
 * and 1+tan_Ce^2 = 1 + sin_Ce^2 * cos_Cn^2 / (sin_Cn^2 + cos_Cn^2 * cos_Ce^2)
 *                = (sin_Cn^2 + cos_Cn^2 * cos_Ce^2 + sin_Ce^2 * cos_Cn^2) / (sin_Cn^2 + cos_Cn^2 * cos_Ce^2)
 *                = 1. / (sin_Cn^2 + cos_Cn^2 * cos_Ce^2)
 *                = inv_denom_tan_Ce^2
 *
 */
    const double sinh_arg_i = tan_Ce * two_inv_denom_tan_Ce;
    const double cosh_arg_i = two_inv_denom_tan_Ce_square - 1;

    double dCn, dCe;
    Cn += clenS (P->gtu, PROJ_ETMERC_ORDER,
                 sin_arg_r, cos_arg_r, sinh_arg_i, cosh_arg_i,
                 &dCn, &dCe);
    Ce += dCe;
    if (fabs (Ce) <= 2.623395162778) {
        xy.y  = P->Qn * Cn + P->Zb;  /* Northing */
        xy.x  = P->Qn * Ce;          /* Easting  */
    } else {
        //proj_errno_set(P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        xy.x = xy.y = HUGE_VAL;
    }
    return xy;
}


/* Ellipsoidal, inverse */
static PJ_LP exact_e_inv (PJ_XY xy, UTM *P) {
    PJ_LP lp = {0.0,0.0};

    /* normalize N, E */
    double Cn = (xy.y - P->Zb)/P->Qn;
    double Ce = xy.x/P->Qn;

    if (fabs(Ce) <= 2.623395162778) { /* 150 degrees */
        /* norm. N, E -> compl. sph. LAT, LNG */
        const double sin_arg_r  = sin(2*Cn);
        const double cos_arg_r  = cos(2*Cn);

        //const double sinh_arg_i = sinh(2*Ce);
        //const double cosh_arg_i = cosh(2*Ce);
        const double exp_2_Ce = exp(2*Ce);
        const double half_inv_exp_2_Ce = 0.5 / exp_2_Ce;
        const double sinh_arg_i = 0.5 * exp_2_Ce - half_inv_exp_2_Ce;
        const double cosh_arg_i = 0.5 * exp_2_Ce + half_inv_exp_2_Ce;

        double dCn_ignored, dCe;
        Cn += clenS(P->utg, PROJ_ETMERC_ORDER,
                    sin_arg_r, cos_arg_r, sinh_arg_i, cosh_arg_i,
                    &dCn_ignored, &dCe);
        Ce += dCe;

        /* compl. sph. LAT -> Gaussian LAT, LNG */
        const double sin_Cn = sin (Cn);
        const double cos_Cn = cos (Cn);

#if 0
        // Non-optimized version:
        double sin_Ce, cos_Ce;
        Ce = atan (sinh (Ce));  // Replaces: Ce = 2*(atan(exp(Ce)) - FORTPI);
        sin_Ce = sin (Ce);
        cos_Ce = cos (Ce);
        Ce     = atan2 (sin_Ce, cos_Ce*cos_Cn);
        Cn     = atan2 (sin_Cn*cos_Ce,  hypot (sin_Ce, cos_Ce*cos_Cn));
#else
/*
 *      One can divide both member of Ce = atan2(...) by cos_Ce, which gives:
 *      Ce     = atan2 (tan_Ce, cos_Cn) = atan2(sinh(Ce), cos_Cn)
 *
 *      and the same for Cn = atan2(...)
 *      Cn     = atan2 (sin_Cn, hypot (sin_Ce, cos_Ce*cos_Cn)/cos_Ce)
 *             = atan2 (sin_Cn, hypot (sin_Ce/cos_Ce, cos_Cn))
 *             = atan2 (sin_Cn, hypot (tan_Ce, cos_Cn))
 *             = atan2 (sin_Cn, hypot (sinhCe, cos_Cn))
 */
        const double sinhCe = sinh (Ce);
        Ce     = atan2 (sinhCe, cos_Cn);
        const double modulus_Ce = hypot (sinhCe, cos_Cn);
        Cn     = atan2 (sin_Cn, modulus_Ce);
#endif

        /* Gaussian LAT, LNG -> ell. LAT, LNG */

        // Optimization of the computation of cos(2*Cn) and sin(2*Cn)
        const double tmp = 2 * modulus_Ce / (sinhCe * sinhCe + 1);
        const double sin_2_Cn  = sin_Cn * tmp;
        const double cos_2_Cn  = tmp * modulus_Ce - 1.;
        //const double cos_2_Cn = cos(2 * Cn);
        //const double sin_2_Cn = sin(2 * Cn);

        lp.phi = gatg (P->cgb,  PROJ_ETMERC_ORDER, Cn, cos_2_Cn, sin_2_Cn);
        lp.lam = Ce;
    }
    else {
        //proj_errno_set(P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        lp.phi = lp.lam = HUGE_VAL;
    }
    return lp;
}


int utm_get_zone(double lon) {
    long zone = lround((floor ((adjlon (lon * M_PI / 180) + M_PI) * 30. / M_PI)));
    if (zone < 0)
        zone = 0;
    else if (zone >= 60)
        zone = 59;
    return zone + 1;
}


int utm_setup(UTM *P, int zone, int south) {
    if (zone > 0 && zone <= 60)
        --zone;
    else
        return 1;

    P->a = 6378137.0;
    P->es = 0.006694379990;

    P->y0 = south ? 10000000. : 0.;
    P->x0 = 500000.;

    P->lam0 = (zone + .5) * M_PI / 30. - M_PI;
    P->k0 = 0.9996;
    P->phi0 = 0.;

    /* third flattening */
    const double n = 0.00167922038638370426;
    double np = n;

    /* COEF. OF TRIG SERIES GEO <-> GAUSS */
    /* cgb := Gaussian -> Geodetic, KW p190 - 191 (61) - (62) */
    /* cbg := Geodetic -> Gaussian, KW p186 - 187 (51) - (52) */
    /* PROJ_ETMERC_ORDER = 6th degree : Engsager and Poder: ICC2007 */

    P->cgb[0] = n*( 2 + n*(-2/3.0  + n*(-2      + n*(116/45.0 + n*(26/45.0 +
                n*(-2854/675.0 ))))));
    P->cbg[0] = n*(-2 + n*( 2/3.0  + n*( 4/3.0  + n*(-82/45.0 + n*(32/45.0 +
                n*( 4642/4725.0))))));
    np     *= n;
    P->cgb[1] = np*(7/3.0 + n*( -8/5.0  + n*(-227/45.0 + n*(2704/315.0 +
                n*( 2323/945.0)))));
    P->cbg[1] = np*(5/3.0 + n*(-16/15.0 + n*( -13/9.0  + n*( 904/315.0 +
                n*(-1522/945.0)))));
    np     *= n;
    /* n^5 coeff corrected from 1262/105 -> -1262/105 */
    P->cgb[2] = np*( 56/15.0  + n*(-136/35.0 + n*(-1262/105.0 +
                n*( 73814/2835.0))));
    P->cbg[2] = np*(-26/15.0  + n*(  34/21.0 + n*(    8/5.0   +
                n*(-12686/2835.0))));
    np     *= n;
    /* n^5 coeff corrected from 322/35 -> 332/35 */
    P->cgb[3] = np*(4279/630.0 + n*(-332/35.0 + n*(-399572/14175.0)));
    P->cbg[3] = np*(1237/630.0 + n*( -12/5.0  + n*( -24832/14175.0)));
    np     *= n;
    P->cgb[4] = np*(4174/315.0 + n*(-144838/6237.0 ));
    P->cbg[4] = np*(-734/315.0 + n*( 109598/31185.0));
    np     *= n;
    P->cgb[5] = np*(601676/22275.0 );
    P->cbg[5] = np*(444337/155925.0);

    /* Constants of the projections */
    /* Transverse Mercator (UTM, ITM, etc) */
    np = n*n;
    /* Norm. mer. quad, K&W p.50 (96), p.19 (38b), p.5 (2) */
    P->Qn = P->k0/(1 + n) * (1 + np*(1/4.0 + np*(1/64.0 + np/256.0)));
    /* coef of trig series */
    /* utg := ell. N, E -> sph. N, E,  KW p194 (65) */
    /* gtu := sph. N, E -> ell. N, E,  KW p196 (69) */
    P->utg[0] = n*(-0.5  + n*( 2/3.0 + n*(-37/96.0 + n*( 1/360.0 +
                n*(  81/512.0 + n*(-96199/604800.0))))));
    P->gtu[0] = n*( 0.5  + n*(-2/3.0 + n*(  5/16.0 + n*(41/180.0 +
                n*(-127/288.0 + n*(  7891/37800.0 ))))));
    P->utg[1] = np*(-1/48.0 + n*(-1/15.0 + n*(437/1440.0 + n*(-46/105.0 +
                n*( 1118711/3870720.0)))));
    P->gtu[1] = np*(13/48.0 + n*(-3/5.0  + n*(557/1440.0 + n*(281/630.0 +
                n*(-1983433/1935360.0)))));
    np      *= n;
    P->utg[2] = np*(-17/480.0 + n*(  37/840.0 + n*(  209/4480.0  +
                n*( -5569/90720.0 ))));
    P->gtu[2] = np*( 61/240.0 + n*(-103/140.0 + n*(15061/26880.0 +
                n*(167603/181440.0))));
    np      *= n;
    P->utg[3] = np*(-4397/161280.0 + n*(  11/504.0 + n*( 830251/7257600.0)));
    P->gtu[3] = np*(49561/161280.0 + n*(-179/168.0 + n*(6601661/7257600.0)));
    np     *= n;
    P->utg[4] = np*(-4583/161280.0 + n*(  108847/3991680.0));
    P->gtu[4] = np*(34729/80640.0  + n*(-3418889/1995840.0));
    np     *= n;
    P->utg[5] = np*(-20648693/638668800.0);
    P->gtu[5] = np*(212378941/319334400.0);

    /* Gaussian latitude value of the origin latitude */
    const double Z = gatg (P->cbg, PROJ_ETMERC_ORDER, P->phi0, cos(2*P->phi0), sin(2*P->phi0));

    /* Origin northing minus true northing at the origin latitude */
    /* i.e. true northing = N - P->Zb                         */
    P->Zb  = - P->Qn*(Z + clens(P->gtu, PROJ_ETMERC_ORDER, 2*Z));

    return 0;
}


PJ_XY utm_forward (PJ_LONLAT coords, UTM *P) {
    PJ_LP p_lp = {coords.lon * DEG_TO_RAD - P->lam0, coords.lat * DEG_TO_RAD};
    PJ_XY result = exact_e_fwd(p_lp, P);
    result.x = result.x * P->a + P->x0;
    result.y = result.y * P->a + P->y0;
    return result;
}


PJ_LONLAT utm_inverse (PJ_XY xy, UTM *P) {
    PJ_XY p_xy = {(xy.x - P->x0) / P->a, (xy.y - P->y0) / P->a};
    PJ_LP p_lp = exact_e_inv(p_xy, P);
    PJ_LONLAT result = {(p_lp.lam + P->lam0) * RAD_TO_DEG, p_lp.phi * RAD_TO_DEG};
    return result;
}
