#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

setup(
    name='simpleproj',
    version='0.1',
    description='Simple geodetic coordinate transformations',
    author='gumblex',
    license='MIT',
    keywords='gis',
    setup_requires=['cffi>=1.11.5'],
    packages=['simpleproj'],
    cffi_modules=[
        'simpleproj/_sp_build.py:ffibuilder',
    ],
    install_requires=['cffi>=1.11.5'],
    platforms='any',
)
