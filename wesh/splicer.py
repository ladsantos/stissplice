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
    filename (``str``):
        Name of the ``*_x1d.fits`` file containing the spectrum.

    prefix (``str``):
        Path to the ``*_x1d.fits`` file containing the spectrum.

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

    # We index the spectral regions as a `dict` in order to avoid confusion with
    # too many numerical indexes
    spectrum = [{'wavelength': wavelength[i],
                 'flux': flux[i],
                 'uncertainty': uncertainty[i],
                 'data_quality': data_quality[i]} for i in range(n_orders)]
    return spectrum


# Identify overlapping regions in each order
def find_overlap(order_pair):
    """
    Find and return the overlapping sections of the Echelle spectrum.

    Parameters
    ----------
    order_pair (Sequence):
        The pair of spectral sections containing two ``dict`` objects, one for
        each order. Can be either an array, list, or sequence.

    Returns
    -------
    sections (``list``):
        List of sections in the following order: unique sections for the first
        and second spectral order followed by the overlapping sections for the
        first and second order. "First" and "second" here do not refer to the
        actual order numbers of the spectrum, but rather in the pair of
        neighboring orders.
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
                 'uncertainty': order_0['uncertainty'][overlap_index_0],
                 'data_quality': order_0['data_quality'][overlap_index_0]
                 }
    overlap_1 = {'wavelength': order_1['wavelength'][overlap_index_1],
                 'flux': order_1['flux'][overlap_index_1],
                 'uncertainty': order_1['uncertainty'][overlap_index_1],
                 'data_quality': order_1['data_quality'][overlap_index_1]
                 }
    unique_0 = {'wavelength': order_0['wavelength'][i0:],
                'flux': order_0['flux'][i0:],
                'uncertainty': order_0['uncertainty'][i0:],
                'data_quality': order_0['data_quality'][i0:]
                }
    unique_1 = {'wavelength': order_1['wavelength'][:i1],
                'flux': order_1['flux'][:i1],
                'uncertainty': order_1['uncertainty'][:i1],
                'data_quality': order_1['data_quality'][:i1]
                }
    sections = [unique_0, unique_1, overlap_0, overlap_1]

    return sections


# Merge overlapping sections
def merge_overlap(overlap_0, overlap_1, inconsistency_sigma=3, outlier_sigma=5,
                  correct_inconsistent_fluxes=True, 
                  correct_outlier_fluxes=True,
                  acceptable_dq_flags=(0, 64, 128, 1024, 2048)):
    """
    Merges overlapping spectral regions. The basic workflow of this function
    is to interpolate the two sections into a common wavelength table and
    calculate the weighted mean flux for each wavelength bin. If the fluxes are
    inconsistent between each other, the code can use the flux with higher SNR
    instead of the mean. If there are still outlier fluxes (compared to
    neighboring pixels), the code uses the flux from the lower SNR section
    instead. Co-added (merged) pixels will have their DQ flag set to `32768` if
    they are the result of combining good pixels (according to the list of
    acceptable flags). Their DQ flag will be set to `65536` if the combined
    pixels do not have acceptable DQ flag.
    
    Parameters
    ----------
    overlap_0 (``dict``):
        Dictionary containing the overlapping spectrum of the first order in a
        pair of neighboring orders.

    overlap_1 (``dict``):
        Dictionary containing the overlapping spectrum of the second order in a
        pair of neighboring orders.

    inconsistency_sigma (``float`` or ``int``, optional):
        Threshold standard deviation to determine if two fluxes in the same
        wavelength bin are inconsistent with each other. Default is ``3``.

    outlier_sigma (``float`` or ``int``, optional):
        Threshold standard deviation to determine if the flux in a given pixel
        is inconsistent with the average between the neighboring pixels. Default
        is ``5``.

    correct_inconsistent_fluxes (``bool``, optional):
        Parameter that decides whether to correct or not correct inconsistent
        fluxes. Default is ``True``.

    correct_outlier_fluxes (``bool``, optional):
        Parameter that decides whether to correct or not correct outlier fluxes.
        Default is ``True``.

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
    # First we need to determine which spectrum has a lower SNR
    avg_snr = np.array([np.mean(overlap_0['flux'] / overlap_0['uncertainty']),
                        np.mean(overlap_1['flux'] / overlap_1['uncertainty'])])

    # We interpolate the higher-SNR spectrum to the wavelength bins of the lower
    # SNR spectrum. This is to avoid degrading the spectrum that already has
    # the worst quality
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
    # step. We will take the weighted averaged, with weights equal to the
    # uncertainties squared. We multiply everything by a large number to avoid
    # numerical issues.
    scale = 1E10
    weights_interp = 1 / (err_interp * scale) ** 2
    weights_ref = 1 / (overlap_ref['uncertainty'] * scale) ** 2

    # We create a new data-quality array filled with 32768, which is what we
    # establish as the flag for co-added pixels
    dq_merge = np.ones_like(f_interp, dtype=int) * 32768

    # Here we deal with the data-quality flags. We only accept flags that are
    # listed in `acceptable_dq_flags`. Let's initialize the dq flag arrays
    dq_ref = overlap_ref['data_quality']
    # The interpolated dq flag array is a bit more involved. First we
    # interpolate the original array to the new wavelength grid, and then we
    # round all the values to the nearest integer
    dq_interp = np.rint(np.interp(overlap_ref['wavelength'],
                        overlap_interp['wavelength'],
                        overlap_interp['data_quality']))
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
        dq_weights_ref[np.where(dq_ref == adq)[0]] = 1
        dq_weights_interp[np.where(dq_interp == adq)[0]] = 1

    # Now we need to verify if we are setting the dq weighting to zero in both
    # the reference and the interpolated dqs. If this is the case, we will
    # set their weights to one and then flag these pixels
    sum_dq_weights = np.copy(dq_weights_ref + dq_weights_interp)
    dq_weights_ref[sum_dq_weights < 1] = 1
    dq_weights_interp[sum_dq_weights < 1] = 1
    dq_merge[sum_dq_weights < 1] = 65536

    # And then we multiply the original weights by the dq weights
    weights_interp *= dq_weights_interp
    weights_ref *= dq_weights_ref

    # This following array will be important later
    sum_weights = weights_interp + weights_ref

    # Finally co-add the overlaps
    wl_merge = np.copy(overlap_ref['wavelength'])
    f_merge = (f_interp * weights_interp +
               overlap_ref['flux'] * weights_ref) / sum_weights
    err_merge = ((err_interp ** 2 * weights_interp ** 2 +
                  overlap_ref['uncertainty'] ** 2 * weights_ref ** 2) ** 0.5) \
        / sum_weights

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

    # There could still be pixels with discrepant fluxes if the higher SNR
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
                      'uncertainty': err_merge, 'data_quality': dq_merge}

    return overlap_merged


# Splice the spectra
def splice(unique_spectra_list, merged_overlap_list):
    """
    Concatenate the unique and the (merged) overlapping spectra.

    Parameters
    ----------
    unique_spectra_list (``list``):
        List of unique spectra.

    merged_overlap_list (``list``):
        List of merged overlapping spectra.

    Returns
    -------
    spliced_wavelength (``numpy.ndarray``):
        Array containing the wavelengths in the entire spectrum.

    spliced_flux (``numpy.ndarray``):
        Array containing the fluxes in the entire spectrum.

    spliced_uncertainty (``numpy.ndarray``):
        Array containing the flux uncertainties in the entire spectrum.
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
    spliced_data_quality = \
        np.concatenate([spectrum['data_quality'] for spectrum in all_spectra])

    return spliced_wavelength, spliced_flux, spliced_uncertainty, \
        spliced_data_quality


# The splice pipeline does everything
def splice_pipeline(dataset, prefix='./', update_fits=False, output_file=None,
                    inconsistency_sigma=3, outlier_sigma=5,
                    correct_inconsistent_fluxes=True,
                    correct_outlier_fluxes=True,
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

    inconsistency_sigma (``float`` or ``int``, optional):
        Threshold standard deviation to determine if two fluxes in the same
        wavelength bin are inconsistent with each other. Default is ``3``.

    outlier_sigma (``float`` or ``int``, optional):
        Threshold standard deviation to determine if the flux in a given pixel
        is inconsistent with the average between the neighboring pixels. Default
        is ``5``.

    correct_inconsistent_fluxes (``bool``, optional):
        Parameter that decides whether to correct or not correct inconsistent
        fluxes. Default is ``True``.

    correct_outlier_fluxes (``bool``, optional):
        Parameter that decides whether to correct or not correct outlier fluxes.
        Default is ``True``.

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
    n_sections = len(sections)

    # Separate spectral sections into pairs
    pairs = [[sections[j], sections[j + 1]] for j in range(n_sections - 1)]
    n_pairs = len(pairs)

    # Identify unique and overlapping sections pair by pair
    unique_sections = []
    overlap_pairs = []
    for i in range(n_pairs):
        splices = find_overlap(pairs[i])
        pairs[i][0] = splices[0]
        pairs[i][1] = splices[1]
        unique_sections.append(splices[0])
        if i == n_pairs - 1:
            unique_sections.append(splices[1])
        else:
            pairs[i + 1][0] = splices[1]
        overlap_pairs.append([splices[2], splices[3]])

    # Merge the overlapping spectral sections
    merged_sections = [
        merge_overlap(overlap_pairs[k][0], overlap_pairs[k][1],
                      inconsistency_sigma, outlier_sigma,
                      correct_inconsistent_fluxes, correct_outlier_fluxes,
                      acceptable_dq_flags)
        for k in range(len(overlap_pairs))
    ]

    # By now we have two lists: unique_sections and merged_sections. The next
    # step is to concatenate everything in the correct order. Since the spectra
    # are listed in reverse order in the `x1d` file, we un-reverse them here
    unique_sections.reverse()
    merged_sections.reverse()

    # Finally splice the unique and merged sections
    wavelength, flux, uncertainty, dq = splice(unique_sections, merged_sections)
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
