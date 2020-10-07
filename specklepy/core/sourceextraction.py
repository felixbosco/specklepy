from IPython import embed
import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from photutils import DAOStarFinder, IRAFStarFinder

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger


def extract_sources(image, noise_threshold, fwhm, star_finder='DAO', image_var=None, background_subtraction=True,
                    write_to=None, debug=True):
    """Extract sources from an image with a StarFinder routine.

    Long description...

    Args:
        image (np.ndarray or str):
            Image array or the name of a file containing the image array.
        noise_threshold (float):
            Multiple of the uncertainty/ standard deviation of the image.
        fwhm (float):
            Expected full width at half maximum (FWHM) of the sources in units of pixels.
        star_finder (str, optional):
            Choose whether the 'DAO' or 'IRAF' StarFinder implementations from photutils shall be used. Default is
            'DAO'.
        image_var (float or str):
            Variance of the image used for the StarFinder threshold (=noise_threshold * sqrt(image_var)). If not
            provided, the code extracts this value from sigma clipped stats. If provided as str-type, the code tries to
            use this as a key to the FITS file HDU list.
        background_subtraction (bool, optional):
            Let the StarFinder consider the background subtraction. Set False for ignoring background flux. Default is
            `True`.
        write_to (str, optional):
            If provided as a str, the list of identified sources is saved to this file.
        debug (bool, optional):
            Show debugging information. Default is `False`.

    Returns:
        sources (astropy.table.Table): Table of identified sources, None if no
            sources are detected.
    """

    # Set logger level
    if debug:
        logger.setLevel('DEBUG')

    # Input parameters
    if isinstance(image, np.ndarray):
        filename = 'current cube'
    elif isinstance(image, str):
        logger.info("The argument image '{}' is interpreted as file name.".format(image))
        filename = image
        image = fits.getdata(filename)
        image = image.squeeze()
    else:
        raise SpecklepyTypeError('extract_sources()', argname='image', argtype=type(image),
                                 expected='np.ndarray or str')

    # Prepare noise statistics
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    logger.info(f"Noise statistics for {filename}:\n\tMean = {mean:.3}\n\tMedian = {median:.3}\n\tStdDev = {std:.3}")

    # Set detection threshold
    if image_var is None:
        threshold = noise_threshold * std
    else:
        if isinstance(image_var, str):
            # Try to load variance extension from file
            image_var = fits.getdata(filename, image_var)
            image_var = np.mean(image_var)
        threshold = noise_threshold * np.sqrt(image_var)

    # Set sky background
    if background_subtraction:
        logger.info(f"Considering mean sky background of {mean}")
        sky = mean
    else:
        sky = 0.0

    # Instantiate StarFinder object
    if not isinstance(star_finder, str):
        raise SpecklepyTypeError('extract_sources', argname='starfinder', argtype=type(star_finder), expected='str')
    if 'dao' in star_finder.lower():
        star_finder = DAOStarFinder(fwhm=fwhm, threshold=threshold, sky=sky)
    elif 'iraf' in star_finder.lower():
        star_finder = IRAFStarFinder(fwhm=fwhm, threshold=threshold, sky=sky)
    else:
        raise SpecklepyValueError('extract_sources', argname='star_finder', argvalue=star_finder,
                                  expected="'DAO' or 'IRAF")

    # Find stars
    logger.info("Extracting sources...")
    sources = star_finder(image)

    # Reformatting sources table
    sources.sort('flux', reverse=True)
    sources.rename_column('xcentroid', 'x')
    sources.rename_column('ycentroid', 'y')
    sources.keep_columns(['x', 'y', 'flux'])

    # Add terminal output
    logger.info(f"Extracted {len(sources)} sources")
    logger.debug(sources)

    # Save sources table to file, if requested
    if write_to is not None:
        logger.info("Writing list of sources to file {}".format(write_to))
        sources.write(write_to, format='ascii.fixed_width', overwrite=True)

    return sources
