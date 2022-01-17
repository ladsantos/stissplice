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


__all__ = ["read_spectrum", "find_overlap", "merge_overlap", "splice",
           "splice_pipeline"]


# Read the Echelle spectrum based on dataset name and a prefix for file location
def read_spectrum(filename, prefix):
    """
    This is a fairly straightforward function to read the spectrum from a `x1d`
    FITS file.

    Parameters
    ----------
    filename (`str`):
        Name of the `*_x1d.fits` file containing the spectrum.

    prefix (`str`):
        Path to the `*_x1d.fits` file containing the spectrum.

    Returns
    -------
    spectrum (`list`):
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
    n_orders = len(wavelength)

    # We index the spectral regions as a `dict` in order to avoid confusion with
    # too many numerical indexes
    spectrum = [{'wavelength': wavelength[i],
                 'flux': flux[i],
                 'uncertainty': uncertainty[i]} for i in range(n_orders)]
    return spectrum


# Identify overlapping regions in each order
def find_overlap(order_pair):
    """

    Parameters
    ----------
    order_pair

    Returns
    -------
    sections
    """
    order_0, order_1 = order_pair

    # Identify the wavelength values in the borders of each order
    borders_0 = \
        np.array([min(order_0['wavelength']), max(order_0['wavelength'])])
    borders_1 = \
        np.array([min(order_1['wavelength']), max(order_1['wavelength'])])

    # Identify the indexes where the orders overlap
    i0 = tools.nearest_index(order_0['wavelength'], borders_1[1])
    i1 = tools.nearest_index(order_1['wavelength'], borders_0[0])
    overlap_index_0 = np.arange(0, i0, 1)
    overlap_index_1 = np.arange(i1, 1024, 1)

    # Break down the order pair into four sections: two are unique spectral
    # sections, and the other two are the overlapping spectral sections
    overlap_0 = {'wavelength': order_0['wavelength'][overlap_index_0],
                 'flux': order_0['flux'][overlap_index_0],
                 'uncertainty': order_0['uncertainty'][overlap_index_0]}
    overlap_1 = {'wavelength': order_1['wavelength'][overlap_index_1],
                 'flux': order_1['flux'][overlap_index_1],
                 'uncertainty': order_1['uncertainty'][overlap_index_1]}
    unique_0 = {'wavelength': order_0['wavelength'][i0:],
                'flux': order_0['flux'][i0:],
                'uncertainty': order_0['uncertainty'][i0:]}
    unique_1 = {'wavelength': order_1['wavelength'][:i1],
                'flux': order_1['flux'][:i1],
                'uncertainty': order_1['uncertainty'][:i1]}
    sections = [unique_0, unique_1, overlap_0, overlap_1]

    return sections


# Merge overlapping sections
def merge_overlap(overlap_0, overlap_1, inconsistency_sigma=3, outlier_sigma=5,
                  correct_inconsistent_fluxes=True, 
                  correct_outlier_fluxes=True):
    """
    
    Parameters
    ----------
    overlap_0
    overlap_1
    inconsistency_sigma
    outlier_sigma
    correct_inconsistent_fluxes
    correct_outlier_fluxes

    Returns
    -------

    """
    # First we need to determine which spectrum has a lower SNR
    avg_snr = np.array([np.mean(overlap_0['flux'] / overlap_0['uncertainty']),
                        np.mean(overlap_1['flux'] / overlap_1['uncertainty'])])

    # We interpolate the higher-SNR spectrum to the wavelength bins of the lower
    # SNR spectrum. This is to avoid degrading the spectrum that is already 
    # the worst
    if avg_snr[0] > avg_snr[1]:
        overlap_interp = overlap_0
        overlap_ref = overlap_1
    else:
        overlap_ref = overlap_0
        overlap_interp = overlap_1
    f_interp = np.interp(overlap_ref['wavelength'],
                         overlap_interp['wavelength'], overlap_interp['flux'])
    err_interp = np.interp(overlap_ref['wavelength'],
                           overlap_interp['wavelength'],
                           overlap_interp['uncertainty'])

    # Check if the fluxes are consistent with each other
    delta_f = f_interp - overlap_ref['flux']
    err_delta = (err_interp ** 2 + overlap_ref['uncertainty'] ** 2) ** 0.5
    inconsistent = np.where(abs(delta_f) > inconsistency_sigma * err_delta)[0]

    # Merge the spectra, the inconsistent pixels will be corrected in the next
    # step
    wl_merge = np.copy(overlap_ref['wavelength'])
    f_merge = (f_interp + overlap_ref['flux']) / 2
    err_merge = ((err_interp ** 2 + overlap_ref['uncertainty'] ** 2) ** 0.5) / 2

    # For the pixels where we have inconsistencies, we use the ones with 
    # the highest SNR
    snr_interp = f_interp / err_interp
    snr_ref = overlap_ref['flux'] / overlap_ref['uncertainty']
    if correct_inconsistent_fluxes is True:
        for i in inconsistent:
            if snr_interp[i] > snr_ref[i]:
                f_merge[i] = f_interp[i]
                err_merge[i] = err_interp[i]
            else:
                f_merge[i] = overlap_ref['flux'][i]
                err_merge[i] = overlap_ref['uncertainty'][i]
    else:
        pass

    # There could still be pixels with discrepant fluxes, if the higher SNR
    # one was an outlier. So we check and correct for them. The next code lines
    # are not really neat and nice to read, so bear with me for now.
    surround_mean = np.array([(f_merge[i] + f_merge[i + 2]) / 2
                              for i in range(len(f_merge) - 2)])
    sm_err = np.array([((err_merge[i] ** 2 + err_merge[i + 2] ** 2) ** 0.5) / 2
                       for i in range(len(f_merge) - 2)])
    surround_diff = abs(f_merge[1:-1] - surround_mean)
    s_diff_err = ((err_merge[1:-1] ** 2 + sm_err ** 2) ** 0.5) / 2
    outlier_indexes = \
        np.where(surround_diff > outlier_sigma * s_diff_err)[0] + 1
    # For the outliers, we use the flux from the lower SNR spectrum instead
    if correct_outlier_fluxes is True:
        for i in outlier_indexes:
            if snr_interp[i] > snr_ref[i]:
                f_merge[i] = overlap_ref['flux'][i]
                err_merge[i] = overlap_ref['uncertainty'][i]
            else:
                f_merge[i] = f_interp[i]
                err_merge[i] = err_interp[i]
    else:
        pass

    overlap_merged = {'wavelength': wl_merge, 'flux': f_merge,
                      'uncertainty': err_merge}

    return overlap_merged


# Splice the spectra
def splice(unique_spectra_list, merged_overlap_list):
    """

    Parameters
    ----------
    unique_spectra_list
    merged_overlap_list

    Returns
    -------

    """
    n_overlap = len(merged_overlap_list)
    all_spectra = []
    for i in range(n_overlap):
        all_spectra.append(unique_spectra_list[i])
        all_spectra.append(merged_overlap_list[i])
    all_spectra.append(unique_spectra_list[-1])
    spliced_wavelength = \
        np.concatenate([spectrum['wavelength'] for spectrum in all_spectra])
    spliced_flux = \
        np.concatenate([spectrum['flux'] for spectrum in all_spectra])
    spliced_uncertainty = \
        np.concatenate([spectrum['uncertainty'] for spectrum in all_spectra])

    return spliced_wavelength, spliced_flux, spliced_uncertainty


# The splice pipeline does everything
def splice_pipeline(dataset, prefix='./', update_fits=False, output_file=None,
                    inconsistency_sigma=3, outlier_sigma=5,
                    correct_inconsistent_fluxes=True,
                    correct_outlier_fluxes=True):
    """

    Parameters
    ----------
    dataset
    prefix
    update_fits
    output_file
    inconsistency_sigma
    outlier_sigma
    correct_inconsistent_fluxes
    correct_outlier_fluxes

    Returns
    -------

    """
    # Read the data
    sections = read_spectrum(dataset, prefix)
    n_sections = len(sections)

    # Separate spectral sections into pairs
    pairs = [[sections[i], sections[i + 1]] for i in range(n_sections - 1)]
    n_pairs = len(pairs)

    # Identify unique and overlapping sections pair by pair
    unique_sections = []
    overlap_pairs = []
    for i in range(n_pairs):
        splices = find_overlap(pairs[i])
        pairs[i][0] = splices[0]
        pairs[i][1] = splices[1]
        unique_sections.append(pairs[i][0])
        if i == n_pairs - 1:
            unique_sections.append(pairs[i][1])
        else:
            pairs[i + 1][0] = splices[1]
        overlap_pairs.append([splices[2], splices[3]])

    # Merge the overlapping spectral sections
    merged_sections = [
        merge_overlap(overlap_pairs[i][0], overlap_pairs[i][1],
                      inconsistency_sigma, outlier_sigma,
                      correct_inconsistent_fluxes, correct_outlier_fluxes)
        for i in range(len(overlap_pairs))
    ]

    # By now we have two lists: unique_sections and merged_sections. The next
    # step is to concatenate everything in the correct order. Since the spectra
    # are listed in reverse order in the `x1d` file, we unreverse them here
    unique_sections.reverse()
    merged_sections.reverse()

    # Finally splice the unique and merged sections
    wavelength, flux, uncertainty = splice(unique_sections, merged_sections)
    spectrum_dict = \
        {'WAVELENGTH': wavelength, 'FLUX': flux, 'ERROR': uncertainty}
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
