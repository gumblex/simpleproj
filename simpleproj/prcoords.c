#include <math.h>
#include "sp_common.h"

#define UNUSED(x) (void)(x)

// Krasovsky 1940 ellipsoid
#define GCJ_A 6378245
// f = 1/298.3; e^2 = 2*f - f**2
#define GCJ_EE 0.00669342162296594323

// Epsilon to use for "exact" iterations.
#define PRC_EPS 1e-5

// Baidu's artificial deviations
#define BD_DLAT 0.0060
#define BD_DLON 0.0065

// Mean Earth Radius
#define EARTH_R 6371000


static inline int sanity_in_china_p(PJ_LONLAT *coords) {
    return ((0.8293 <= coords->lat) && (coords->lat <= 55.8271) &&
            (72.004 <= coords->lon) && (coords->lon <= 137.8347));
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

    double dLat_m = (-100 + 2 * x + 3 * y + 0.2 * y * y + 0.1 * x * y +
        0.2 * sqrt(fabs(x)) + (
            2 * sin(x * 6 * M_PI) + 2 * sin(x * 2 * M_PI) +
            2 * sin(y * M_PI) + 4 * sin(y / 3 * M_PI) +
            16 * sin(y / 12 * M_PI) + 32 * sin(y / 30 * M_PI)
        ) * 20 / 3);
    double dLon_m = (300 + x + 2 * y + 0.1 * x * x + 0.1 * x * y +
        0.1 * sqrt(fabs(x)) + (
            2 * sin(x * 6 * M_PI) + 2 * sin(x * 2 * M_PI) +
            2 * sin(x * M_PI) + 4 * sin(x / 3 * M_PI) +
            15 * sin(x / 12 * M_PI) + 30 * sin(x / 30 * M_PI)
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
    /* rev_transform_rough; accuracy ~2e-6 deg (meter-level) */
    PJ_LONLAT diff = prcoords_wgs_gcj(gcj, check_china);
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
    return prcoords_gcj_wgs(prcoords_bd_gcj(bd, check_china), check_china);
}


PJ_LONLAT prcoords_wgs_bd(PJ_LONLAT bd, int check_china) {
    return prcoords_gcj_bd(prcoords_wgs_gcj(bd, check_china), check_china);
}


/* generic "bored function" factory, Caijun 2014
 * gcj: 4 calls to prcoords_wgs_gcj; ~0.1mm acc */
#define __precise_conv(bad, wgs, fwd, rev) do {  \
    wgs = rev(bad, 0);  \
    PJ_LONLAT diff = {99, 99};  \
    /* Wait till we hit fixed point or get bored */  \
    unsigned i = 0;  \
    while ((i < 10) && hypot(diff.lon, diff.lat) > PRC_EPS) {  \
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

PJ_LONLAT prcoords_gcj_wgs_bored(PJ_LONLAT gcj, int check_china) {
    if (check_china && !sanity_in_china_p(&gcj)) return gcj;
    PJ_LONLAT result;
    __precise_conv(gcj, result, prcoords_wgs_gcj, prcoords_gcj_wgs);
    return result;
}

PJ_LONLAT prcoords_bd_gcj_bored(PJ_LONLAT bd, int check_china) {
    if (check_china && !sanity_in_china_p(&bd)) return bd;
    PJ_LONLAT result;
    __precise_conv(bd, result, prcoords_gcj_bd, prcoords_bd_gcj);
    return result;
}

PJ_LONLAT prcoords_bd_wgs_bored(PJ_LONLAT bd, int check_china) {
    if (check_china && !sanity_in_china_p(&bd)) return bd;
    PJ_LONLAT result;
    __precise_conv(bd, result, prcoords_wgs_bd, prcoords_bd_wgs);
    return result;
}
