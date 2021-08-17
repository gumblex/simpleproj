#include <math.h>
#include <float.h>
#include "sp_common.h"

#define UNUSED(x) (void)(x)

// Krasovsky 1940 ellipsoid
#define GCJ_A 6378245
// f = 1/298.3; e^2 = 2*f - f**2
#define GCJ_EE 0.00669342162296594323

// Baidu's artificial deviations
#define BD_DLAT 0.0060
#define BD_DLON 0.0065

// Baidu Mercator uses EPSG:7008, Clarke 1866
// PROJ:
// +proj=merc +a=6378206.4 +b=6356583.8 +lat_ts=0.0 +lon_0=0.0
// +x_0=0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs
#define BD_A 6378206.4
#define BD_E 0.0822718542230032587696


static inline int sanity_in_china_p(PJ_LONLAT *coords) {
    return ((0.8293 <= coords->lat) && (coords->lat <= 55.8271) &&
            (72.004 <= coords->lon) && (coords->lon <= 137.8347));
}


static inline int sanity_in_china_p_gcj(PJ_LONLAT *coords) {
    return ((0.83056 <= coords->lat) && (coords->lat <= 55.8296) &&
            (72.0077 <= coords->lon) && (coords->lon <= 137.8437));
}


static inline int sanity_in_china_p_bd(PJ_LONLAT *coords) {
    return ((0.83676 <= coords->lat) && (coords->lat <= 55.8355) &&
            (72.0142 <= coords->lon) && (coords->lon <= 137.8503));
}


PJ_LONLAT prcoords_wgs_gcj(PJ_LONLAT wgs, int check_china) {
    if (check_china && !sanity_in_china_p(&wgs)) {
        return wgs;
    }

    double x = wgs.lon - 105;
    double y = wgs.lat - 35;

    /*
    These distortion functions accept (x = lon - 105, y = lat - 35).
    They return distortions in terms of arc lengths, in meters.

    In other words, you can pretty much figure out how much you will be off
    from WGS-84 just through evaulating them...

    For example, at the (mapped) center of China (105E, 35N), you get a
    default deviation of <300, -100> meters.
    */

    double x_pi = x * M_PI, y_pi = y * M_PI, sq_x = sqrt(fabs(x));

    double dLat_m = (-100 + 2 * x + 3 * y + 0.2 * y * y + 0.1 * x * y +
        0.2 * sq_x + (
            2 * sin(x_pi * 6) + 2 * sin(x_pi * 2) +
            2 * sin(y_pi) + 4 * sin(y_pi / 3) +
            16 * sin(y_pi / 12) + 32 * sin(y_pi / 30)
        ) * 20 / 3);
    double dLon_m = (300 + x + 2 * y + 0.1 * x * x + 0.1 * x * y +
        0.1 * sq_x + (
            2 * sin(x_pi * 6) + 2 * sin(x_pi * 2) +
            2 * sin(x_pi) + 4 * sin(x_pi / 3) +
            15 * sin(x_pi / 12) + 30 * sin(x_pi / 30)
        ) * 20 / 3);

    double radLat = wgs.lat * DEG_TO_RAD;
    double magic = 1 - GCJ_EE * pow(sin(radLat), 2); // just a common expr

    // [[:en:Latitude#Length_of_a_degree_of_latitude]]
    double lat_deg_arclen = DEG_TO_RAD * (GCJ_A * (1 - GCJ_EE)) / pow(magic, 1.5);
    // [[:en:Longitude#Length_of_a_degree_of_longitude]]
    double lon_deg_arclen = DEG_TO_RAD * (GCJ_A * cos(radLat) / sqrt(magic));

    // The screwers pack their deviations into degrees and disappear.
    // Note how they are mixing WGS-84 and Krasovsky 1940 ellipsoids here...
    PJ_LONLAT result = {
        wgs.lon + (dLon_m / lon_deg_arclen),
        wgs.lat + (dLat_m / lat_deg_arclen)
    };
    return result;
}


PJ_LONLAT prcoords_gcj_wgs(PJ_LONLAT gcj, int check_china) {
    if (check_china && !sanity_in_china_p_gcj(&gcj)) return gcj;

    /* rev_transform_rough; accuracy ~2e-6 deg (meter-level) */
    PJ_LONLAT diff = prcoords_wgs_gcj(gcj, 0);
    diff.lon -= gcj.lon;
    diff.lat -= gcj.lat;
    PJ_LONLAT result = {
        gcj.lon - diff.lon,
        gcj.lat - diff.lat
    };
    return result;
}


PJ_LONLAT prcoords_gcj_bd(PJ_LONLAT gcj, int check_china) {
    UNUSED(check_china);
    double x = gcj.lon;
    double y = gcj.lat;

    // trivia: pycoordtrans actually describes how these values are calculated
    double r = sqrt(x * x + y * y) + 0.00002 * sin(y * DEG_TO_RAD * 3000);
    double theta = atan2(y, x) + 0.000003 * cos(x * DEG_TO_RAD * 3000);

    // Hard-coded default deviations again!
    PJ_LONLAT result = {
        r * cos(theta) + BD_DLON,
        r * sin(theta) + BD_DLAT
    };
    return result;
}


// Yes, we can implement a "precise" one too.
PJ_LONLAT prcoords_bd_gcj(PJ_LONLAT bd, int check_china) {
    // accuracy ~1e-7 deg (decimeter-level; exceeds usual data accuracy)
    UNUSED(check_china);
    double x = bd.lon - BD_DLON;
    double y = bd.lat - BD_DLAT;

    // trivia: pycoordtrans actually describes how these values are calculated
    double r = sqrt(x * x + y * y) - 0.00002 * sin(y * DEG_TO_RAD * 3000);
    double theta = atan2(y, x) - 0.000003 * cos(x * DEG_TO_RAD * 3000);

    PJ_LONLAT result = {r * cos(theta), r * sin(theta)};
    return result;
}


PJ_LONLAT prcoords_bd_wgs(PJ_LONLAT bd, int check_china) {
    if (check_china && !sanity_in_china_p_bd(&bd)) return bd;
    return prcoords_gcj_wgs(prcoords_bd_gcj(bd, 0), 0);
}


PJ_LONLAT prcoords_wgs_bd(PJ_LONLAT bd, int check_china) {
    return prcoords_gcj_bd(prcoords_wgs_gcj(bd, check_china), check_china);
}


/* generic "exact function" factory, Caijun 2014
 * gcj: 4 calls to prcoords_wgs_gcj; ~0.1mm acc */
#define __precise_conv(bad, wgs, fwd, rev) do {  \
    wgs = rev(bad, 0);  \
    PJ_LONLAT diff = {99, 99};  \
    /* Wait till we hit fixed point or get bored */  \
    unsigned i = 0;  \
    while ((i < 10) && hypot(diff.lon, diff.lat) > DBL_EPSILON) {  \
        diff = fwd(wgs, 0);  \
        diff.lon -= bad.lon; diff.lat -= bad.lat;  \
        wgs.lon -= diff.lon; wgs.lat -= diff.lat;  \
        i += 1;  \
    }  \
} while (0)


/*
Precise functions using caijun 2014 method

Why "bored"? Because they usually exceed source data accuracy -- the
original GCJ implementation contains noise from a linear-modulo PRNG,
and Baidu seems to do similar things with their API too.
*/

PJ_LONLAT prcoords_gcj_wgs_exact(PJ_LONLAT gcj, int check_china) {
    if (check_china && !sanity_in_china_p_gcj(&gcj)) return gcj;
    PJ_LONLAT result;
    __precise_conv(gcj, result, prcoords_wgs_gcj, prcoords_gcj_wgs);
    return result;
}

PJ_LONLAT prcoords_bd_gcj_exact(PJ_LONLAT bd, int check_china) {
    if (check_china && !sanity_in_china_p_bd(&bd)) return bd;
    PJ_LONLAT result;
    __precise_conv(bd, result, prcoords_gcj_bd, prcoords_bd_gcj);
    return result;
}

PJ_LONLAT prcoords_bd_wgs_exact(PJ_LONLAT bd, int check_china) {
    if (check_china && !sanity_in_china_p_bd(&bd)) return bd;
    PJ_LONLAT result;
    __precise_conv(bd, result, prcoords_wgs_bd, prcoords_bd_wgs);
    return result;
}


/****************************************************************************
* Convert tau' = sinh(psi) = tan(chi) to tau = tan(phi).  The code is taken
* from GeographicLib::Math::tauf(taup, e).
*
* Here
*   phi = geographic latitude (radians)
* psi is the isometric latitude
*   psi = asinh(tan(phi)) - e * atanh(e * sin(phi))
*       = asinh(tan(chi))
* chi is the conformal latitude
*
* The representation of latitudes via their tangents, tan(phi) and tan(chi),
* maintains full *relative* accuracy close to latitude = 0 and +/- pi/2.
* This is sometimes important, e.g., to compute the scale of the transverse
* Mercator projection which involves cos(phi)/cos(chi) tan(phi)
*
* From Karney (2011), Eq. 7,
*
*   tau' = sinh(psi) = sinh(asinh(tan(phi)) - e * atanh(e * sin(phi)))
*        = tan(phi) * cosh(e * atanh(e * sin(phi))) -
*          sec(phi) * sinh(e * atanh(e * sin(phi)))
*        = tau * sqrt(1 + sigma^2) - sqrt(1 + tau^2) * sigma
* where
*   sigma = sinh(e * atanh( e * tau / sqrt(1 + tau^2) ))
*
* For e small,
*
*    tau' = (1 - e^2) * tau
*
* The relation tau'(tau) can therefore by reliably inverted by Newton's
* method with
*
*    tau = tau' / (1 - e^2)
*
* as an initial guess.  Newton's method requires dtau'/dtau.  Noting that
*
*   dsigma/dtau = e^2 * sqrt(1 + sigma^2) /
*                 (sqrt(1 + tau^2) * (1 + (1 - e^2) * tau^2))
*   d(sqrt(1 + tau^2))/dtau = tau / sqrt(1 + tau^2)
*
* we have
*
*   dtau'/dtau = (1 - e^2) * sqrt(1 + tau'^2) * sqrt(1 + tau^2) /
*                (1 + (1 - e^2) * tau^2)
*
* This works fine unless tau^2 and tau'^2 overflows.  This may be partially
* cured by writing, e.g., sqrt(1 + tau^2) as hypot(1, tau).  However, nan
* will still be generated with tau' = inf, since (inf - inf) will appear in
* the Newton iteration.
*
* If we note that for sufficiently large |tau|, i.e., |tau| >= 2/sqrt(eps),
* sqrt(1 + tau^2) = |tau| and
*
*   tau' = exp(- e * atanh(e)) * tau
*
* So
*
*   tau = exp(e * atanh(e)) * tau'
*
* can be returned unless |tau| >= 2/sqrt(eps); this then avoids overflow
* problems for large tau' and returns the correct result for tau' = +/-inf
* and nan.
*
* Newton's method usually take 2 iterations to converge to double precision
* accuracy (for WGS84 flattening).  However only 1 iteration is needed for
* |chi| < 3.35 deg.  In addition, only 1 iteration is needed for |chi| >
* 89.18 deg (tau' > 70), if tau = exp(e * atanh(e)) * tau' is used as the
* starting guess.
****************************************************************************/

static inline double pj_sinhpsi2tanphi(const double taup, const double e) {
    static const int numit = 5;
    // min iterations = 1, max iterations = 2; mean = 1.954
    static const double rooteps = sqrt(DBL_EPSILON);
    static const double tol = rooteps / 10; // the criterion for Newton's method
    static const double tmax = 2 / rooteps; // threshold for large arg limit exact
    const double e2m = 1 - e * e;
    const double abs_taup = fabs(taup);
    const double stol = tol * (1.0 < abs_taup ? abs_taup : 1.0);
    // The initial guess.  70 corresponds to chi = 89.18 deg (see above)
    double tau = abs_taup > 70 ? taup * exp(e * atanh(e)) : taup / e2m;
    if (!(fabs(tau) < tmax))      // handles +/-inf and nan and e = 1
        return tau;
    // If we need to deal with e > 1, then we could include:
    // if (e2m < 0) return std::numeric_limits<double>::quiet_NaN();
    int i = numit;
    for (; i; --i) {
        double tau1 = sqrt(1 + tau * tau);
        double sig = sinh( e * atanh(e * tau / tau1) );
        double taupa = sqrt(1 + sig * sig) * tau - sig * tau1;
        double dtau = ( (taup - taupa) * (1 + e2m * (tau * tau)) /
                        (e2m * tau1 * sqrt(1 + taupa * taupa)) );
        tau += dtau;
        if (!(fabs(dtau) >= stol))  // backwards test to allow nans to succeed.
            break;
    }
    // if (i == 0)
        // proj_context_errno_set(ctx, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    return tau;
}


PJ_XY prcoords_bd_merc_fwd(PJ_LONLAT bd) {
    double phi = bd.lat * DEG_TO_RAD;
    // Instead of calling tan and sin, call sin and cos which the compiler
    // optimizes to a single call to sincos.
    double sphi = sin(phi);
    double cphi = cos(phi);
    PJ_XY result = {
        bd.lon * DEG_TO_RAD * BD_A,
        (asinh(sphi/cphi) - BD_E * atanh(BD_E * sphi)) * BD_A
    };
    return result;
}


PJ_LONLAT prcoords_bd_merc_inv(PJ_XY point) {
    PJ_LONLAT result = {
        point.x * RAD_TO_DEG / BD_A,
        atan(pj_sinhpsi2tanphi(sinh(point.y / BD_A), BD_E)) * RAD_TO_DEG
    };
    return result;
}
