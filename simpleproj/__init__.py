#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
from ._simpleproj import ffi, lib

EARTH_EQUATORIAL_RADIUS = 6378137


class UTM:
    """UTM projection converter."""
    def __init__(self, zone, south=False):
        self.zone = zone
        self.south = south
        self._ctx = ffi.new("UTM *")
        result = lib.utm_setup(self._ctx, zone, bool(south))
        if result:
            raise ValueError("invalid zone %d" % zone)

    def forward(self, lon, lat):
        """Convert (lon, lat) to (x, y)."""
        result = lib.utm_forward({'lon': lon, 'lat': lat}, self._ctx)
        return (result.x, result.y)

    def inverse(self, x, y):
        """Convert (x, y) to (lon, lat)."""
        result = lib.utm_inverse({'x': x, 'y': y}, self._ctx)
        return (result.lon, result.lat)

    @classmethod
    def get_zone(self, lon):
        """Get UTM Zone by longitude."""
        return lib.utm_get_zone(lon)


def wgs_gcj(lon, lat, check_china=True):
    """Convert WGS84 coordinates to GCJ02."""
    result = lib.prcoords_wgs_gcj({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_wgs(lon, lat, check_china=True):
    """Convert GCJ02 coordinates to WGS84 (fast estimation)."""
    result = lib.prcoords_gcj_wgs({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_bd(lon, lat):
    """Convert GCJ02 coordinates to BD09."""
    result = lib.prcoords_gcj_bd({'lon': lon, 'lat': lat}, 0)
    return (result.lon, result.lat)


def bd_gcj(lon, lat):
    """Convert BD09 coordinates to GCJ02 (fast estimation)."""
    result = lib.prcoords_bd_gcj({'lon': lon, 'lat': lat}, 0)
    return (result.lon, result.lat)


def bd_wgs(lon, lat, check_china=True):
    """Convert BD09 coordinates to WGS84 (fast estimation)."""
    result = lib.prcoords_bd_wgs({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def wgs_bd(lon, lat, check_china=True):
    """Convert WGS84 coordinates to BD09."""
    result = lib.prcoords_wgs_bd({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_wgs_exact(lon, lat, check_china=True):
    """Convert GCJ02 coordinates to WGS84 accurately."""
    result = lib.prcoords_gcj_wgs_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def bd_gcj_exact(lon, lat, check_china=True):
    """Convert BD09 coordinates to GCJ02 accurately."""
    result = lib.prcoords_bd_gcj_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def bd_wgs_exact(lon, lat, check_china=True):
    """Convert BD09 coordinates to WGS84 accurately."""
    result = lib.prcoords_bd_wgs_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


gcj_wgs_bored = gcj_wgs_exact
bd_gcj_bored = bd_gcj_exact
bd_wgs_bored = bd_wgs_exact


def bd_merc_fwd(lon, lat):
    """Convert BD09 (lon, lat) to Baidu Mercator projection."""
    result = lib.prcoords_bd_merc_fwd({'lon': lon, 'lat': lat})
    return (result.x, result.y)


def bd_merc_inv(x, y):
    """Convert Baidu Mercator projection (x, y) to BD09 (lon, lat)."""
    result = lib.prcoords_bd_merc_inv({'x': x, 'y': y})
    return (result.lon, result.lat)


def earth_distance(lon1, lat1, lon2, lat2):
    return lib.sp_distance(lon1, lat1, lon2, lat2)


def from3857_to4326(x, y):
    lon = math.degrees(x / EARTH_EQUATORIAL_RADIUS)
    lat = math.degrees(math.atan(math.sinh(y/EARTH_EQUATORIAL_RADIUS)))
    return (lon, lat)


def from4326_to3857(lon, lat):
    x = math.radians(lon) * EARTH_EQUATORIAL_RADIUS
    y = math.asinh(math.tan(math.radians(lat))) * EARTH_EQUATORIAL_RADIUS
    return (x, y)

