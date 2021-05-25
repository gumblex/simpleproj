#include <math.h>
#include "geodesic.h"
#include "distance.h"

double sp_distance(double lon1, double lat1, double lon2, double lat2) {
    const double a = 6378137.0;
    const double rf = 298.257223563;
    double s12;
    struct geod_geodesic gd;
    geod_init(&gd, a, 1 / rf);
    geod_inverse(&gd, lat1, lon1, lat2, lon2, &s12, 0, 0);
    return s12;
}
