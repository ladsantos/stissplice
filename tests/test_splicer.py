#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains useful tests to splice Echelle spectra from the Hubble
Space Telescope.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from wesh import splicer, tools


# Test the entire splicing pipeline
def test_pipeline(precision_threshold=1E-6):
    dataset = 'oblh01040'
    prefix = '../data/'

    spectrum_table = splicer.splice_pipeline(dataset, prefix)
    i0 = tools.nearest_index(spectrum_table['WAVELENGTH'].data, 2310)
    test = spectrum_table['FLUX'][i0]
    assert(abs(test - 2.460101E-12) / test < precision_threshold)
