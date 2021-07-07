#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
from ._simpleproj import ffi, lib

EARTH_EQUATORIAL_RADIUS = 6378137


class UTM:
    def __init__(self, zone, south=False):
        self.zone = zone
        self.south = south
        self._ctx = ffi.new("UTM *")
        result = lib.utm_setup(self._ctx, zone, bool(south))
        if result:
            raise ValueError("invalid zone %d" % zone)

    def forward(self, lon, lat):
        result = lib.utm_forward({'lon': lon, 'lat': lat}, self._ctx)
        return (result.x, result.y)

    def inverse(self, x, y):
        result = lib.utm_inverse({'x': x, 'y': y}, self._ctx)
        return (result.lon, result.lat)

    @classmethod
    def get_zone(self, lon):
        return lib.utm_get_zone(lon)


def wgs_gcj(lon, lat, check_china=True):
    result = lib.prcoords_wgs_gcj({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_wgs(lon, lat, check_china=True):
    result = lib.prcoords_gcj_wgs({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_bd(lon, lat):
    result = lib.prcoords_gcj_bd({'lon': lon, 'lat': lat}, 0)
    return (result.lon, result.lat)


def bd_gcj(lon, lat):
    result = lib.prcoords_bd_gcj({'lon': lon, 'lat': lat}, 0)
    return (result.lon, result.lat)


def bd_wgs(lon, lat, check_china=True):
    result = lib.prcoords_bd_wgs({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def wgs_bd(lon, lat, check_china=True):
    result = lib.prcoords_wgs_bd({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def gcj_wgs_exact(lon, lat, check_china=True):
    result = lib.prcoords_gcj_wgs_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def bd_gcj_exact(lon, lat, check_china=True):
    result = lib.prcoords_bd_gcj_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


def bd_wgs_exact(lon, lat, check_china=True):
    result = lib.prcoords_bd_wgs_exact({'lon': lon, 'lat': lat}, check_china)
    return (result.lon, result.lat)


gcj_wgs_bored = gcj_wgs_exact
bd_gcj_bored = bd_gcj_exact
bd_wgs_bored = bd_wgs_exact


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

