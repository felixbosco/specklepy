import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, IRAFStarFinder

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logging


def find_sources(image, noise_threshold, fwhm, starfinder='DAO', background_subtraction=True, writeto=None, verbose=True):
    """Finding sources in an image with a StarFinder routine.

    Long description...

    Args:
        image (np.ndarray):
        noise_threshold (float):
        fwhm (float):
        starfinder (str, optional): Choose whether the 'DAO' or 'IRAF'
            StarFinder implementations from photutils shall be used. Default
            is 'DAO'.
        background_subtraction (bool, optional): Default is True.
        writeto (str, optional): If provided as a str, the list of identified
            sources  is saved to this file.
        verbose (bool, optional): Set to False, if reducing the terminal output.
            Default is True.

    Returns:
        sources (astropy.table.Table): Table of identified sources, None if no
            sources are detected.
    """

    # Input parameters
    if not isinstance(image, np.ndarray):
        if isinstance(image, str):
            logging.info("The argument image '{}' is interpreted as file name.".format(image))
            filename = image
            image = fits.getdata(filename)
        else:
            raise SpecklepyTypeError('find_sources()', argname='image', argtype=type(image), expected='np.ndarray or str')
            # raise TypeError("The function sourceextraction.find_stars received \
            #                     argument image of type {}, but needs to be \
            #                     either np.ndarray or str type!".format(type(image)))
    else:
        filename = 'current cube'

    # Prepare noise statistics
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    logging.info("Noise statistics for {}:\n\tMean = {:.3}\n\tMedian = {:.3}\n\tStdDev = {:.3}".format(filename, mean, median, std))

    # Instantiate starfinder object
    if not isinstance(starfinder, str):
        raise TypeError("The function sourceextraction.find_stars received \
                            argument starfinder of type {}, but needs to be \
                            str type!".format(type(image)))
    if 'dao' in starfinder.lower():
        starfinder = DAOStarFinder(fwhm=fwhm, threshold=noise_threshold*std)
    elif 'iraf' in starfinder.lower():
        starfinder = IRAFStarFinder(fwhm=fwhm, threshold=noise_threshold*std)
    else:
        raise ValueError("The function sourceextraction.find_stars received \
                            argument starfinder of {}, but needs to be \
                            either 'DAO' or 'IRAF'!".format(type(image)))

    # Find stars
    # daofind = DAOStarFinder(fwhm=starfinder_fwhm, threshold=noise_threshold*std)
    if background_subtraction:
        logging.info("Subtracting background...")
        image -= median
    logging.info("Finding sources...")
    sources = starfinder(image)

    # Reformatting sources table
    sources.sort('flux', reverse=True)
    sources.rename_column('xcentroid', 'x')
    sources.rename_column('ycentroid', 'y')
    sources.keep_columns(['x', 'y', 'flux'])
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    if verbose:
        logging.info("Found {} sources:\n{}".format(len(sources), sources))
    else:
        logging.info("Found {} sources".format(len(sources)))

    # Save sources table to file, if writeto is provided
    if writeto is not None:
        logging.info("Writing list of sources to file {}".format(writeto))
        sources.write(writeto, format='ascii', overwrite=True)

    return sources
