#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains useful tools to splice Echelle spectra from the Hubble
Space Telescope.
"""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import numpy as np
from wesh import tools
from astropy.io import fits
from astropy.table import Table


__all__ = ["read_spectrum", "find_overlap_all", "find_overlap_pair",
           "find_overlap_trio", "merge_overlap", "splice", "splice_pipeline"]


# Read the Echelle spectrum based on dataset name and a prefix for file location
def read_spectrum(filename, prefix, truncate_edge_left=None,
                  truncate_edge_right=None):
    """
    This is a fairly straightforward function to read the spectrum from a `x1d`
    FITS file.

    Parameters
    ----------
    filename (``str``):
        Name of the ``*_x1d.fits`` file containing the spectrum.

    prefix (``str``):
        Path to the ``*_x1d.fits`` file containing the spectrum.

    truncate_edge_left (``int``, optional):
        Set the number of low-resolution pixels at the left edge of the detector
        where the spectra should be truncated. If ``None``, then no truncation
        is applied. Default is ``None``.

    truncate_edge_right (``int``, optional):
        Set the number of low-resolution pixels at the right edge of the
        detector where the spectra should be truncated. If ``None``, then no
        truncation is applied. Default is ``None``.

    Returns
    -------
    spectrum (``list``):
        List of all the orders contained in the Echelle spectrum and their
        respective fluxes.
    """
    with fits.open(prefix + '%s_x1d.fits' % filename) as hdu:
        header = hdu[0].header
        data = hdu['SCI'].data
    optical_element = header['OPT_ELEM']
    if optical_element[0] != 'E':
        raise TypeError("This is not an Echelle spectrum.")
    wavelength = data['WAVELENGTH']
    flux = data['FLUX']
    uncertainty = data['ERROR']
    data_quality = data['DQ']
    n_orders = len(wavelength)

    tel = truncate_edge_left
    if truncate_edge_right is not None:
        ter = -truncate_edge_right
    else:
        ter = truncate_edge_right

    # We index the spectral regions as a `dict` in order to avoid confusion with
    # too many numerical indexes. Also, since the orders are in reverse order,
    # we index them in the opposite way
    spectrum = [{'wavelength': wavelength[-i - 1][tel:ter],
                 'flux': flux[-i - 1][tel:ter],
                 'uncertainty': uncertainty[-i - 1][tel:ter],
                 'data_quality': data_quality[-i - 1][tel:ter]}
                for i in range(n_orders)]
    return spectrum


# Identify overlaps in the whole spectrum
def find_overlap_all(spectrum):
    """
    Find and return the overlapping sections of the Echelle spectrum.

    Parameters
    ----------
    spectrum (``list``):
        List of dictionaries containing the orders of the Echelle spectrum.
        It should resemble the output of ``read_spectrum()``.

    Returns
    -------
    unique_sections (``list``):
        List containing the unique sections of the spectrum.

    overlap_pair_sections (``list``):
        List containing the overlapping pairs of the spectrum.

    overlap_trio_sections (``list``):
        List containing the overlapping trios of the spectrum.
    """
    n_orders = len(spectrum)

    # Identify the wavelength borders of each order
    borders = []
    for order in spectrum:
        borders.append([min(order['wavelength']), max(order['wavelength'])])

    unique_sections = []
    first_trio = []
    second_trio = []
    third_trio = []
    first_pair = []
    second_pair = []

    # The following code is hacky and not pretty to look at, but it works. Sorry
    # for the mess!

    # First we deal with the first orders
    order = spectrum[0]
    wl = order['wavelength']
    idx = [tools.nearest_index(wl, borders[1][0]),
           tools.nearest_index(wl, borders[2][0])]

    # There is always a unique section here
    unique_idx = np.arange(0, idx[0], 1)
    unique_sections.append({'wavelength': order['wavelength'][unique_idx],
                            'flux': order['flux'][unique_idx],
                            'uncertainty': order['uncertainty'][unique_idx],
                            'data_quality': order['data_quality'][unique_idx]})

    # There is also a pair overlap. But since it's with the next order, it's
    # considered a "second pair"
    overlap_01 = np.arange(idx[0], idx[1], 1)
    second_pair.append({'wavelength': order['wavelength'][overlap_01],
                        'flux': order['flux'][overlap_01],
                        'uncertainty': order['uncertainty'][overlap_01],
                        'data_quality': order['data_quality'][overlap_01]})

    if idx[1] < 1023:
        # There is a trio overlap
        overlap_012 = np.arange(idx[1], 1023, 1)
        third_trio.append(
            {'wavelength': order['wavelength'][overlap_012],
             'flux': order['flux'][overlap_012],
             'uncertainty': order['uncertainty'][overlap_012],
             'data_quality': order['data_quality'][overlap_012]}
        )
    else:
        pass

    # Now the second order
    order = spectrum[1]
    wl = order['wavelength']
    idx = [tools.nearest_index(wl, borders[2][0]),
           tools.nearest_index(wl, borders[0][1]),
           tools.nearest_index(wl, borders[3][0])]

    # There are two pairs, potentially one or two trios, and potentially a
    # unique section
    if idx[1] > idx[0]:
        # There is a trio overlap
        overlap_01 = np.arange(0, idx[0], 1)
        overlap_012 = np.arange(idx[0], idx[1], 1)
        unique_1 = None

        if idx[2] < 1023:
            # There is another trio overlap
            overlap_12 = np.arange(idx[1], idx[2], 1)
            overlap_123 = np.arange(idx[2], 1023, 1)
        else:
            overlap_123 = None
            overlap_12 = np.arange(idx[0], idx[2], 1)

    else:
        # No trio overlap
        overlap_01 = np.arange(0, idx[1], 1)
        overlap_012 = None
        unique_1 = np.arange(idx[1], idx[0], 1)
        overlap_123 = None
        overlap_12 = np.arange(idx[0], idx[2], 1)

    # Add the to the lists
    if unique_1 is not None:
        unique_sections.append(
            {'wavelength': order['wavelength'][unique_1],
             'flux': order['flux'][unique_1],
             'uncertainty': order['uncertainty'][unique_1],
             'data_quality': order['data_quality'][unique_1]}
        )
    else:
        pass

    first_pair.append({'wavelength': order['wavelength'][overlap_01],
                       'flux': order['flux'][overlap_01],
                       'uncertainty': order['uncertainty'][overlap_01],
                       'data_quality': order['data_quality'][overlap_01]})
    second_pair.append({'wavelength': order['wavelength'][overlap_12],
                        'flux': order['flux'][overlap_12],
                        'uncertainty': order['uncertainty'][overlap_12],
                        'data_quality': order['data_quality'][overlap_12]})

    if overlap_012 is not None:
        second_trio.append({'wavelength': order['wavelength'][overlap_012],
                            'flux': order['flux'][overlap_012],
                            'uncertainty': order['uncertainty'][overlap_012],
                            'data_quality': order['data_quality'][overlap_012]})
    else:
        pass

    if overlap_123 is not None:
        third_trio.append({'wavelength': order['wavelength'][overlap_123],
                           'flux': order['flux'][overlap_123],
                           'uncertainty': order['uncertainty'][overlap_123],
                           'data_quality': order['data_quality'][overlap_123]})

    # Now we deal with the third to third before last orders in a loop
    for i in range(n_orders - 4):
        order = spectrum[i + 2]
        wl = order['wavelength']
        idx = [tools.nearest_index(wl, borders[i][1]),
               tools.nearest_index(wl, borders[i + 3][0]),
               tools.nearest_index(wl, borders[i + 1][1]),
               tools.nearest_index(wl, borders[i + 4][0])]
        if idx[0] > 0:
            overlap_idx_012 = np.arange(0, idx[0], 1)
        else:
            overlap_idx_012 = None
        if idx[2] < idx[1]:
            overlap_idx_12 = np.arange(idx[0], idx[2], 1)
            overlap_idx_123 = None
            unique_idx_2 = np.arange(idx[2], idx[1], 1)
            overlap_idx_23 = np.arange(idx[1], idx[3], 1)
        else:
            overlap_idx_12 = np.arange(idx[0], idx[1], 1)
            overlap_idx_123 = np.arange(idx[1], idx[2], 1)
            unique_idx_2 = None
            overlap_idx_23 = np.arange(idx[2], idx[3], 1)
        if idx[3] < 1023:
            overlap_idx_234 = np.arange(idx[3], 1023, 1)
        else:
            overlap_idx_234 = None

        first_pair.append(
            {'wavelength': order['wavelength'][overlap_idx_12],
             'flux': order['flux'][overlap_idx_12],
             'uncertainty': order['uncertainty'][overlap_idx_12],
             'data_quality': order['data_quality'][overlap_idx_12]}
        )
        second_pair.append(
            {'wavelength': order['wavelength'][overlap_idx_23],
             'flux': order['flux'][overlap_idx_23],
             'uncertainty': order['uncertainty'][overlap_idx_23],
             'data_quality': order['data_quality'][overlap_idx_23]}
        )

        if overlap_idx_012 is not None:
            first_trio.append(
                {'wavelength': order['wavelength'][overlap_idx_012],
                 'flux': order['flux'][overlap_idx_012],
                 'uncertainty': order['uncertainty'][overlap_idx_012],
                 'data_quality': order['data_quality'][overlap_idx_012]}
            )
        else:
            pass

        if overlap_idx_123 is not None:
            second_trio.append(
                {'wavelength': order['wavelength'][overlap_idx_123],
                 'flux': order['flux'][overlap_idx_123],
                 'uncertainty': order['uncertainty'][overlap_idx_123],
                 'data_quality': order['data_quality'][overlap_idx_123]}
            )
        else:
            pass

        if overlap_idx_234 is not None:
            third_trio.append(
                {'wavelength': order['wavelength'][overlap_idx_234],
                 'flux': order['flux'][overlap_idx_234],
                 'uncertainty': order['uncertainty'][overlap_idx_234],
                 'data_quality': order['data_quality'][overlap_idx_234]}
            )
        else:
            pass

        if unique_idx_2 is not None:
            unique_sections.append(
                {'wavelength': order['wavelength'][unique_idx_2],
                 'flux': order['flux'][unique_idx_2],
                 'uncertainty': order['uncertainty'][unique_idx_2],
                 'data_quality': order['data_quality'][unique_idx_2]}
            )
        else:
            pass

    # Now we deal with the last orders. Almost there!
    order = spectrum[-2]
    wl = order['wavelength']
    idx = [tools.nearest_index(wl, borders[-4][1]),
           tools.nearest_index(wl, borders[-1][0]),
           tools.nearest_index(wl, borders[-3][1])]

    # There are two pairs, potentially one or two trios, and potentially a
    # unique section
    if idx[0] > 0:
        # There is a trio overlap
        overlap_012 = np.arange(0, idx[0], 1)
        if idx[2] > idx[1]:
            # There is another trio overlap
            overlap_123 = np.arange(idx[1], idx[2], 1)
            overlap_12 = np.arange(idx[0], idx[1], 1)
            overlap_23 = np.arange(idx[2], 1023, 1)
            unique_2 = None
        else:
            overlap_123 = None
            unique_2 = np.arange(idx[2], idx[1], 1)
            overlap_12 = np.arange(idx[0], idx[2], 1)
            overlap_23 = np.arange(idx[1], 1023, 1)
    else:
        overlap_012 = None
        overlap_123 = None
        unique_2 = np.arange(idx[2], idx[1], 1)
        overlap_12 = np.arange(idx[0], idx[2], 1)
        overlap_23 = np.arange(idx[1], 1023, 1)

    # Add the to the lists
    if unique_2 is not None:
        unique_sections.append(
            {'wavelength': order['wavelength'][unique_2],
             'flux': order['flux'][unique_2],
             'uncertainty': order['uncertainty'][unique_2],
             'data_quality': order['data_quality'][unique_2]}
        )
    else:
        pass

    first_pair.append({'wavelength': order['wavelength'][overlap_12],
                       'flux': order['flux'][overlap_12],
                       'uncertainty': order['uncertainty'][overlap_12],
                       'data_quality': order['data_quality'][overlap_12]})
    second_pair.append({'wavelength': order['wavelength'][overlap_23],
                        'flux': order['flux'][overlap_23],
                        'uncertainty': order['uncertainty'][overlap_23],
                        'data_quality': order['data_quality'][overlap_23]})

    if overlap_012 is not None:
        first_trio.append({'wavelength': order['wavelength'][overlap_012],
                           'flux': order['flux'][overlap_012],
                           'uncertainty': order['uncertainty'][overlap_012],
                           'data_quality': order['data_quality'][overlap_012]})
    else:
        pass

    if overlap_123 is not None:
        second_trio.append({'wavelength': order['wavelength'][overlap_123],
                            'flux': order['flux'][overlap_123],
                            'uncertainty': order['uncertainty'][overlap_123],
                            'data_quality': order['data_quality'][overlap_123]})

    # Finally deal with the last order
    order = spectrum[-1]
    wl = order['wavelength']
    idx = [tools.nearest_index(wl, borders[-3][1]),
           tools.nearest_index(wl, borders[-2][1])]

    # There is always a unique section here
    unique_idx = np.arange(idx[1], 1023, 1)
    unique_sections.append(
        {'wavelength': order['wavelength'][unique_idx],
         'flux': order['flux'][unique_idx],
         'uncertainty': order['uncertainty'][unique_idx],
         'data_quality': order['data_quality'][unique_idx]}
    )

    # There is also a pair overlap. But since it's with the previous order, it's
    # considered a "first pair"
    overlap_23 = np.arange(idx[0], idx[1], 1)
    first_pair.append(
        {'wavelength': order['wavelength'][overlap_23],
         'flux': order['flux'][overlap_23],
         'uncertainty': order['uncertainty'][overlap_23],
         'data_quality': order['data_quality'][overlap_23]}
    )

    if idx[0] > 0:
        # There is a trio overlap
        overlap_123 = np.arange(0, idx[0], 1)
        first_trio.append(
            {'wavelength': order['wavelength'][overlap_123],
             'flux': order['flux'][overlap_123],
             'uncertainty': order['uncertainty'][overlap_123],
             'data_quality': order['data_quality'][overlap_123]}
        )
    else:
        pass

    # With all that done, we assemble the overlap sections into a large list
    overlap_pair_sections = []
    overlap_trio_sections = []
    n_pairs = len(first_pair)
    n_trios = len(first_trio)
    for i in range(n_pairs):
        overlap_pair_sections.append([first_pair[i], second_pair[i]])
    if n_trios > 0:
        for i in range(n_trios):
            overlap_trio_sections.append([first_trio[i], second_trio[i],
                                          third_trio[i]])

    return unique_sections, overlap_pair_sections, overlap_trio_sections


# Identify overlapping regions in each order
def find_overlap_pair(order_pair):
    """
    Find and return the overlapping sections of the Echelle spectrum.

    Parameters
    ----------
    order_pair (Sequence):
        The pair of spectral sections containing two ``dict`` objects, one for
        each order. Can be either an array, list, or sequence.

    Returns
    -------
    unique_sections (``list``):
        List of unique sections in for the first and second spectral orders in a
        pair, respectively. "First" and "second" here do not refer to the actual
        order numbers of the spectrum, but rather in the pair of neighboring
        orders.

    overlap_sections (``list``):
        List of overlap sections in for the first order overlap with
        the second, and the second order overlap with the first, respectively.
        "First" and "second" here do not refer to the actual order numbers of
        the spectrum, but rather in the pair of neighboring orders.
    """
    order_0, order_1 = order_pair

    # Identify the wavelength values in the borders of each order
    borders_0 = \
        np.array([min(order_0['wavelength']), max(order_0['wavelength'])])
    borders_1 = \
        np.array([min(order_1['wavelength']), max(order_1['wavelength'])])

    # Identify the indexes where the orders overlap
    i0 = tools.nearest_index(order_0['wavelength'], borders_1[0])
    i1 = tools.nearest_index(order_1['wavelength'], borders_0[1])
    len_order_0 = len(order_0['wavelength'])
    len_order_1 = len(order_1['wavelength'])
    overlap_index_0 = np.arange(i0, len_order_0, 1)
    overlap_index_1 = np.arange(0, i1, 1)

    # Break down the order pair into four sections: two are unique spectral
    # sections, and the other two are the overlapping spectral sections
    overlap_0 = {'wavelength': order_0['wavelength'][overlap_index_0],
                 'flux': order_0['flux'][overlap_index_0],
                 'uncertainty': order_0['uncertainty'][overlap_index_0],
                 'data_quality': order_0['data_quality'][overlap_index_0]
                 }
    overlap_1 = {'wavelength': order_1['wavelength'][overlap_index_1],
                 'flux': order_1['flux'][overlap_index_1],
                 'uncertainty': order_1['uncertainty'][overlap_index_1],
                 'data_quality': order_1['data_quality'][overlap_index_1]
                 }
    unique_0 = {'wavelength': order_0['wavelength'][:i0],
                'flux': order_0['flux'][:i0],
                'uncertainty': order_0['uncertainty'][:i0],
                'data_quality': order_0['data_quality'][:i0]
                }
    unique_1 = {'wavelength': order_1['wavelength'][i1:],
                'flux': order_1['flux'][i1:],
                'uncertainty': order_1['uncertainty'][i1:],
                'data_quality': order_1['data_quality'][i1:]
                }

    unique_sections = [unique_0, unique_1]
    overlap_sections = [overlap_0, overlap_1]

    return unique_sections, overlap_sections


# Find overlap between three sections
def find_overlap_trio(order_trio):
    """
    Find and return the overlapping sections of a trio of orders in the Echelle
    spectrum.

    Parameters
    ----------
    order_trio (Sequence):
        The trio of spectral sections containing two ``dict`` objects, one for
        each order. Can be either an array, list, or sequence.

    Returns
    -------
    unique_sections (``list``):
        List of unique sections in for the first, second and third spectral
        orders in a trio, respectively. "First", "second" and "third"
        here do not refer to the actual order numbers of the spectrum, but
        rather in the trio of neighboring orders.

    overlap_sections (``list``):
        List of overlap sections in for [0] the first order overlap with
        the second and third, [1] the first order overlap with the second, [2]
        the second order overlap with the third, [3] the second order overlap
        with the first and third, [4] the second order overlap with the first,
        [5] the third order overlap with the second, [6] the third order overlap
        with the first and second, respectively. "First", "second" and "third"
        here do not refer to the actual order numbers of the spectrum, but
        rather in the trio of neighboring orders.
    """
    order_0, order_1, order_2 = order_trio

    # Identify the wavelength values in the borders of each order
    borders_0 = \
        np.array([min(order_0['wavelength']), max(order_0['wavelength'])])
    borders_1 = \
        np.array([min(order_1['wavelength']), max(order_1['wavelength'])])
    borders_2 = \
        np.array([min(order_2['wavelength']), max(order_2['wavelength'])])

    # Identify the indexes where the orders overlap. Things are a bit convoluted
    # here but we need to identify where there are overlaps between 3 orders and
    # also between 2 orders
    i00 = tools.nearest_index(order_0['wavelength'], borders_1[0])
    i01 = tools.nearest_index(order_0['wavelength'], borders_2[0])
    # If i01 is lower than the length of order_0, then there is a trio overlap
    len_order_0 = len(order_0['wavelength'])
    len_order_1 = len(order_1['wavelength'])
    if i01 < len_order_0 - 1:
        i10 = tools.nearest_index(order_1['wavelength'], borders_2[0])
        i11 = tools.nearest_index(order_1['wavelength'], borders_0[1])
    else:
        i10 = tools.nearest_index(order_1['wavelength'], borders_0[1])
        i11 = tools.nearest_index(order_1['wavelength'], borders_2[0])
    i20 = tools.nearest_index(order_2['wavelength'], borders_0[1])
    i21 = tools.nearest_index(order_2['wavelength'], borders_1[1])

    # If i01 is lower than the length of order_0, then there is a trio overlap
    if i01 < len_order_0 - 1:
        overlap_index_0_1 = np.arange(i00, i01, 1)  # Indexes in order zero
        # corresponding to the overlap with order 1
        overlap_index_0_12 = np.arange(i01, len_order_0, 1)  # Indexes in order
        # zero corresponding to the overlap with orders 1 and 2
        overlap_index_1_0 = np.arange(0, i10, 1)  # Indexes in order 1
        # corresponding to the overlap with order 0
        overlap_index_1_02 = np.arange(i10, i11, 1)  # Indexes in order 1
        # corresponding to the overlap with orders 0 and 2
        overlap_index_1_2 = np.arange(i11, len_order_1, 1)  # Indexes in order 1
        # corresponding to the overlap with order 2
        overlap_index_2_01 = np.arange(0, i20, 1)  # Indexes in order 2
        # corresponding to the overlap with orders 1 and 2
        overlap_index_2_1 = np.arange(i20, i21, 1)  # Indexes in order 2
        # corresponding to the overlap with order 1

        # Break down the order trio into ten sections: three are unique spectral
        # sections, and the other seven are the overlapping spectral sections
        overlap_0_12 = {'wavelength': order_0['wavelength'][overlap_index_0_12],
                        'flux': order_0['flux'][overlap_index_0_12],
                        'uncertainty': order_0['uncertainty'][
                            overlap_index_0_12],
                        'data_quality': order_0['data_quality'][
                            overlap_index_0_12]
                        }
        overlap_1_02 = {'wavelength': order_1['wavelength'][overlap_index_1_02],
                        'flux': order_1['flux'][overlap_index_1_02],
                        'uncertainty': order_1['uncertainty'][
                            overlap_index_1_02],
                        'data_quality': order_1['data_quality'][
                            overlap_index_1_02]
                        }
        overlap_2_01 = {'wavelength': order_2['wavelength'][overlap_index_2_01],
                        'flux': order_2['flux'][overlap_index_2_01],
                        'uncertainty': order_2['uncertainty'][
                            overlap_index_2_01],
                        'data_quality': order_2['data_quality'][
                            overlap_index_2_01]
                        }
        unique_1 = None

    # Otherwise, if i01 is 1023, then there is no trio overlap
    else:
        overlap_index_0_1 = np.arange(i00, len_order_0, 1)  # Indexes in order
        # zero corresponding to the overlap with order 1
        overlap_index_1_0 = np.arange(0, i10, 1)  # Indexes in order zero
        # corresponding to the overlap with order 0
        overlap_index_1_2 = np.arange(i11, len_order_1, 1)  # Indexes in order 1
        # corresponding to the overlap with order 2
        overlap_index_2_1 = np.arange(0, i21, 1)  # Indexes in order 2
        # corresponding to the overlap with order 1

        # Break down the order trio into ten sections: three are unique spectral
        # sections, and the other seven are the overlapping spectral sections
        overlap_0_12 = None
        overlap_1_02 = None
        overlap_2_01 = None

        unique_1 = {'wavelength': order_1['wavelength'][i10:i11],
                    'flux': order_1['flux'][i10:i11],
                    'uncertainty': order_1['uncertainty'][i10:i11],
                    'data_quality': order_1['data_quality'][i10:i11]
                    }

    # The overlap_0_1, overlap_1_0, overlap_1_2, overlap_2_1, unique_0 and
    # unique_2 are the same whether there is a trio overlap or not
    overlap_0_1 = {'wavelength': order_0['wavelength'][overlap_index_0_1],
                   'flux': order_0['flux'][overlap_index_0_1],
                   'uncertainty': order_0['uncertainty'][overlap_index_0_1],
                   'data_quality': order_0['data_quality'][
                       overlap_index_0_1]
                   }
    overlap_1_0 = {'wavelength': order_1['wavelength'][overlap_index_1_0],
                   'flux': order_1['flux'][overlap_index_1_0],
                   'uncertainty': order_1['uncertainty'][overlap_index_1_0],
                   'data_quality': order_1['data_quality'][
                       overlap_index_1_0]
                   }
    overlap_1_2 = {'wavelength': order_1['wavelength'][overlap_index_1_2],
                   'flux': order_1['flux'][overlap_index_1_2],
                   'uncertainty': order_1['uncertainty'][overlap_index_1_2],
                   'data_quality': order_1['data_quality'][
                       overlap_index_1_2]
                   }
    overlap_2_1 = {'wavelength': order_2['wavelength'][overlap_index_2_1],
                   'flux': order_2['flux'][overlap_index_2_1],
                   'uncertainty': order_2['uncertainty'][overlap_index_2_1],
                   'data_quality': order_2['data_quality'][
                       overlap_index_2_1]
                   }
    unique_0 = {'wavelength': order_0['wavelength'][:i00],
                'flux': order_0['flux'][:i00],
                'uncertainty': order_0['uncertainty'][i00:],
                'data_quality': order_0['data_quality'][i00:]
                }

    unique_2 = {'wavelength': order_2['wavelength'][i21:],
                'flux': order_2['flux'][i21:],
                'uncertainty': order_2['uncertainty'][i21:],
                'data_quality': order_2['data_quality'][i21:]
                }

    unique_sections = [unique_0, unique_1, unique_2]
    overlap_sections = [overlap_0_1, overlap_0_12, overlap_1_0,
                        overlap_1_02, overlap_1_2, overlap_2_01,
                        overlap_2_1]

    return unique_sections, overlap_sections


# Merge overlapping sections
def merge_overlap(overlap_sections,
                  acceptable_dq_flags=(0, 64, 128, 1024, 2048)):
    """
    Merges overlapping spectral regions. The basic workflow of this function
    is to interpolate the sections into a common wavelength table and calculate
    the weighted mean flux for each wavelength bin. If the fluxes are
    inconsistent between each other, the code can use the flux with higher SNR
    instead of the mean. If there are still outlier fluxes (compared to
    neighboring pixels), the code uses the flux from the lower SNR section
    instead. Co-added (merged) pixels will have their DQ flag set to `32768` if
    they are the result of combining good pixels (according to the list of
    acceptable flags). Their DQ flag will be set to `65536` if the combined
    pixels do not have an acceptable DQ flag.
    
    Parameters
    ----------
    overlap_sections (``list``):
        List of dictionaries containing the overlapping spectra of neighboring
        orders.

    acceptable_dq_flags (array-like, optional):
        Data-quality flags that are acceptable when co-adding overlapping
        spectra. The default values are (0, 64, 128, 1024, 2048), which
        correspond to: 0 = regular pixel, 64 = vignetted pixel, 128 = pixel in
        overscan region, 1024 = small blemish, 2048 = more than 30% of
        background pixels rejected by sigma-clipping in the data reduction.

    Returns
    -------
    overlap_merged (``dict``):
        Dictionary containing the merged overlapping spectrum.
    """
    n_overlaps = len(overlap_sections)

    # First we need to determine which spectrum has a lower SNR
    avg_snr = np.array([np.mean(ok['flux'] / ok['uncertainty'])
                        for ok in overlap_sections])

    # We interpolate the higher-SNR spectra to the wavelength bins of the higher
    # SNR spectrum.
    min_snr_idx = np.where(avg_snr == max(avg_snr))[0][0]
    overlap_ref = overlap_sections.pop(min_snr_idx)

    f_interp = []
    err_interp = []
    for i in range(n_overlaps - 1):
        f_interp.append(np.interp(overlap_ref['wavelength'],
                                  overlap_sections[i]['wavelength'],
                                  overlap_sections[i]['flux']))
        err_interp.append(np.interp(overlap_ref['wavelength'],
                          overlap_sections[i]['wavelength'],
                          overlap_sections[i]['uncertainty']))
    f_interp = np.array(f_interp)
    err_interp = np.array(err_interp)

    # Merge the spectra. We will take the weighted averaged, with weights equal
    # to the inverse of the uncertainties squared multiplied by a scale factor
    # to avoid numerical overflows.
    scale = 1E-20
    weights_interp = (1 / err_interp) ** 2 * scale
    weights_ref = (1 / overlap_ref['uncertainty']) ** 2 * scale

    # Here we deal with the data-quality flags. We only accept flags that are
    # listed in `acceptable_dq_flags`. Let's initialize the dq flag arrays
    dq_ref = overlap_ref['data_quality']

    # We create a new data-quality array filled with 32768, which is what we
    # establish as the flag for co-added pixels
    dq_merge = np.ones_like(dq_ref, dtype=int) * 32768

    # The interpolated dq flag array is a bit more involved. First we
    # interpolate the original array to the new wavelength grid, and then we
    # round all the values to the nearest integer
    dq_interp = []
    for i in range(n_overlaps - 1):
        dq_interp.append(np.rint(np.interp(overlap_ref['wavelength'],
                                           overlap_sections[i]['wavelength'],
                                           overlap_sections[i]['data_quality']))
                         )
    dq_interp = np.array(dq_interp)
    # However this does not guarantee the interpolated and rounded dq values are
    # valid dq flags. Since the interpolation occurs at very small wavelength
    # shifts, for now we assume that all dq flags will be valid. This may be
    # changed in the future.
    # We start assuming that all the dq weights are zero
    dq_weights_ref = np.zeros_like(dq_ref)
    dq_weights_interp = np.zeros_like(dq_interp)
    # And then for each acceptable dq, if the element of the dq array is one
    # of the acceptable flags, we set its dq weight to one
    for adq in acceptable_dq_flags:
        dq_weights_ref[np.where(dq_ref == adq)] = 1
        dq_weights_interp[np.where(dq_interp == adq)] = 1

    # Now we need to verify if we are setting the dq weighting to zero in both
    # the reference and the interpolated dqs. If this is the case, we will
    # set their weights to one and then flag these pixels
    sum_dq_weights = np.copy(dq_weights_ref + np.sum(dq_weights_interp, axis=0))
    dq_weights_ref[sum_dq_weights < 1] = 1
    for i in range(n_overlaps - 1):
        dq_weights_interp[i][sum_dq_weights < 1] = 1
    dq_merge[sum_dq_weights < 1] = 65536

    # And then we multiply the original weights by the dq weights
    weights_interp *= dq_weights_interp
    weights_ref *= dq_weights_ref

    # This following array will be important later
    sum_weights = np.sum(weights_interp, axis=0) + weights_ref

    # Finally co-add the overlaps
    wl_merge = np.copy(overlap_ref['wavelength'])

    f_merge = np.zeros_like(overlap_ref['flux'])
    err_merge = np.zeros_like(overlap_ref['uncertainty'])
    for i in range(n_overlaps - 1):
        f_merge += f_interp[i] * weights_interp[i]
        err_merge += err_interp[i] ** 2 * weights_interp[i] ** 2
    f_merge += overlap_ref['flux'] * weights_ref
    err_merge += overlap_ref['uncertainty'] ** 2 * weights_ref ** 2
    f_merge = f_merge / sum_weights
    err_merge = err_merge ** 0.5 / sum_weights

    overlap_merged = {'wavelength': wl_merge, 'flux': f_merge,
                      'uncertainty': err_merge, 'data_quality': dq_merge}

    return overlap_merged


# Splice the spectra
def splice(unique_spectra_list, merged_pair_list, merged_trio_list):
    """
    Concatenate the unique and the (merged) overlapping spectra.

    Parameters
    ----------
    unique_spectra_list (``list``):
        List of unique spectra.

    merged_pair_list (``list``):
        List of merged overlapping pair spectra.

    merged_trio_list (``list``):
        List of merged overlapping trio spectra.

    Returns
    -------
    spliced_wavelength (``numpy.ndarray``):
        Array containing the wavelengths in the entire spectrum.

    spliced_flux (``numpy.ndarray``):
        Array containing the fluxes in the entire spectrum.

    spliced_uncertainty (``numpy.ndarray``):
        Array containing the flux uncertainties in the entire spectrum.
    """
    n_pair_overlap = len(merged_pair_list)
    all_spectra = []

    # We always start with the first unique spectrum
    all_spectra.append(unique_spectra_list[0])
    k = 1
    for i in range(n_pair_overlap):
        all_spectra.append(merged_pair_list[i])
        try:
            all_spectra.append(merged_trio_list[i])
            pass
        except IndexError:
            all_spectra.append(unique_spectra_list[k])
            k += 1
            pass
    spliced_wavelength = \
        np.concatenate([spectrum['wavelength'] for spectrum in all_spectra])
    spliced_flux = \
        np.concatenate([spectrum['flux'] for spectrum in all_spectra])
    spliced_uncertainty = \
        np.concatenate([spectrum['uncertainty'] for spectrum in all_spectra])
    spliced_data_quality = \
        np.concatenate([spectrum['data_quality'] for spectrum in all_spectra])

    return spliced_wavelength, spliced_flux, spliced_uncertainty, \
        spliced_data_quality


# The splice pipeline does everything
def splice_pipeline(dataset, prefix='./', update_fits=False, output_file=None,
                    truncate_edge_left=None, truncate_edge_right=None,
                    acceptable_dq_flags=(0, 64, 128, 1024, 2048)):
    """
    The main workhorse of the package. This pipeline performs all the steps
    necessary to merge overlapping spectral sections and splice them with the
    unique sections.

    Parameters
    ----------
    dataset (``str``):
        String that identifies the dataset (example: ``'oblh01040'``).

    prefix (``str``):
        Path to the ``*_x1d.fits`` file containing the spectrum.

    update_fits (``bool``, optional):
        THIS FEATURE HAS NOT BEEN TESTED YET. Use carefully, since it can modify
        fits files permanently. Parameter that decides whether to update the
        ``*_x1d.fits`` file with a new extension containing the spliced
        spectrum.

    output_file (``str`` or ``None``, optional):
        String containing the location to save the output spectrum as an ascii
        file. If ``None``, no output file is saved and the code returns an
        Astropy Table instead. Default is ``None``.

    acceptable_dq_flags (array-like, optional):
        Data-quality flags that are acceptable when co-adding overlapping
        spectra. The default values are (0, 64, 128, 1024, 2048), which
        correspond to: 0 = regular pixel, 64 = vignetted pixel, 128 = pixel in
        overscan region, 1024 = small blemish, 2048 = more than 30% of
        background pixels rejected by sigma-clipping in the data reduction.

    Returns
    -------
    spliced_spectrum_table (``astropy.Table`` object):
        Astropy Table containing the spliced spectrum. Only returned if
        ``output_file`` is ``None``.
    """
    # Read the data
    sections = read_spectrum(dataset, prefix)

    unique_sections, overlap_pair_sections, overlap_trio_sections = \
        find_overlap_all(sections)

    # Merge the overlapping spectral sections
    merged_pairs = [
        merge_overlap(overlap_pair_sections[k], acceptable_dq_flags)
        for k in range(len(overlap_pair_sections))
    ]

    if len(overlap_trio_sections) > 0:
        merged_trios = [
            merge_overlap(overlap_trio_sections[k], acceptable_dq_flags)
            for k in range(len(overlap_trio_sections))
        ]
    else:
        merged_trios = []

    # By now we have two lists: unique_sections and merged_sections. The next
    # step is to concatenate everything in the correct order.

    # Finally splice the unique and merged sections
    wavelength, flux, uncertainty, dq = splice(unique_sections, merged_pairs,
                                               merged_trios)

    # Instantiate the spectrum dictionary
    spectrum_dict = \
        {'WAVELENGTH': wavelength, 'FLUX': flux, 'ERROR': uncertainty, 'DQ': dq}
    spliced_spectrum_table = Table(spectrum_dict)

    # This feature has not been tested yet! Use carefully!
    if update_fits is True:
        fits.append(prefix + '%s_x1d.fits' % dataset,
                    data=spliced_spectrum_table)

    # Return or output the result
    if output_file is None:
        return spliced_spectrum_table
    else:
        spliced_spectrum_table.write(output_file, format='ascii')
