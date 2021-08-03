import os
import sys
from cffi import FFI

if sys.platform == "win32":
    # fix mingw builds
    import distutils.cygwinccompiler
    distutils.cygwinccompiler.get_msvcr = lambda: []

ffibuilder = FFI()

ffibuilder.cdef(r"""
typedef struct { double   x,   y; }  PJ_XY;
typedef struct { double lon, lat; }  PJ_LONLAT;

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

PJ_LONLAT prcoords_wgs_gcj(PJ_LONLAT wgs, int check_china);
PJ_LONLAT prcoords_gcj_wgs(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_gcj_bd(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_bd_gcj(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_bd_wgs(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_wgs_bd(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_gcj_wgs_exact(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_bd_gcj_exact(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_bd_wgs_exact(PJ_LONLAT bd, int check_china);
PJ_XY prcoords_bd_merc_fwd(PJ_LONLAT bd);
PJ_LONLAT prcoords_bd_merc_inv(PJ_XY point);

double sp_distance(double lon1, double lat1, double lon2, double lat2);
""")

_curdir = os.path.dirname(__file__)
_sources_path = [
    os.path.join(_curdir, s) for s in
    ('utm.c', 'prcoords.c', 'distance.c', 'geodesic.c')
]

ffibuilder.set_source("simpleproj._simpleproj", '''
#include "distance.h"
#include "utm.h"
#include "prcoords.h"
''', sources=_sources_path, libraries=['m'], extra_compile_args=['-I' + _curdir])

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
