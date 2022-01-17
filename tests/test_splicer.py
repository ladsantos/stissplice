#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains useful tests to splice Echelle spectra from the Hubble
Space Telescope.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
from wesh import splicer
import matplotlib.pyplot as plt
from astropy.io import fits


# Test the entire splicing pipeline
def test_pipeline():
    dataset = 'oblh01040'
    prefix = '../data/'
    wavelength, flux, uncertainty = splicer.splice_pipeline(dataset, prefix)
    print(wavelength)
