#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains useful tests to splice Echelle spectra from the Hubble
Space Telescope.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from stissplice import splicer
import numpy as np


# Test the entire splicing pipeline
def test_pipeline(precision_threshold=1E-6):
    dataset = 'oblh01040'
    prefix = '../data/'

    spectrum_table = splicer.splice_pipeline(dataset, prefix,
                                             weight='sensitivity')

    truth = np.loadtxt('../docs/source/spliced_spectrum_truth.dat', skiprows=1)
    wl_truth = truth[:, 0]
    f_truth = truth[:, 1]
    u_truth = truth[:, 2]

    f_scale = 1E-12
    wl_diff = abs(np.sum(spectrum_table['WAVELENGTH'].data - wl_truth))
    f_diff = abs(np.sum(spectrum_table['FLUX'].data - f_truth)) / f_scale
    u_diff = abs(np.sum(spectrum_table['ERROR'].data - u_truth)) / f_scale
    diff = wl_diff + f_diff + u_diff

    assert(diff < precision_threshold)
