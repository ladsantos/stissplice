#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(
    name="stissplice",
    version="1.0.4",
    author="Leonardo dos Santos",
    author_email="ldsantos@stsci.edu",
    packages=["stissplice"],
    url="https://github.com/ladsantos/stissplice",
    license="MIT",
    description="Splicer of STIS Echelle Spectra from Hubble Space Telescope",
    install_requires=[line.strip() for line in
                      open('requirements.txt', 'r').readlines()],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ]
)
