#ifndef SIMPLEPROD_COMMON_H
#define SIMPLEPROD_COMMON_H

#ifndef M_PI
#define M_PI            3.141592653589793238462643383
#endif
#define M_TWOPI         6.283185307179586476925286767   /* 2*pi */
#define RAD_TO_DEG    57.29577951308232087679815481
#define DEG_TO_RAD   0.01745329251994329576923690768

typedef struct { double   x,   y; }  PJ_XY;
typedef struct { double lam, phi; }  PJ_LP;
typedef struct { double lon, lat; }  PJ_LONLAT;

#endif
