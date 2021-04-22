import numpy as np

from specklepy.logging import logger
from specklepy.core.sourceextraction import SourceExtractor


def extract_sources(image, noise_threshold, fwhm, algorithm='DAO', image_var=None, background_subtraction=True,
                    show=False, select=False, matching_threshold=30, collapse=False, write_to=None, cast_dtype=None,
                    debug=False):
    """Extract sources from an image with a StarFinder routine.

    Long description...

    Args:
        image (np.ndarray or str):
            Image array or the name of a file containing the image array.
        noise_threshold (float):
            Multiple of the uncertainty/ standard deviation of the image.
        fwhm (float):
            Expected full width at half maximum (FWHM) of the sources in units of pixels.
        algorithm (str, optional):
            Choose whether the `'DAO'` or `'IRAF'` StarFinder algorithms from `photutils` shall be used. Also the
            `specklepy` algorithms `'PeakFinder'` and `'ManualFinder'` are available. Default is 'DAO'.
        image_var (float or str):
            Variance of the image used for the StarFinder threshold (=noise_threshold * sqrt(image_var)). If not
            provided, the code extracts this value from sigma clipped stats. If provided as str-type, the code tries to
            use this as a key to the FITS file HDU list.
        background_subtraction (bool, optional):
            Let the StarFinder consider the background subtraction. Set False for ignoring background flux. Default is
            `True`.
        show (bool, optional):
            Create a plot of the identified sources on the image.
        select (bool or dict, optional):
            Create a plot of the identified sources on the image, with the option of selecting apertures for future
            purposes. These are selected by clicking on the image.
        matching_threshold (int, optional):
            Radial threshold for accepting a source as matched or, when manually selecting peaks, the radius of the
            search box.
        collapse (bool, optional):
            Collapse the a data cube along the third axis prior to analysis.
        write_to (str, optional):
            If provided as a str, the list of identified sources is saved to this file.
        cast_dtype (str, optional):
            Str-representation of a data type to cast the image to, prior to the extraction process.
        debug (bool, optional):
            Show debugging information. Default is `False`.

    Returns:
        sources (astropy.table.Table):
            Table of identified sources, `None` if no sources are detected.
        selected (list, optional):
            List of selected sources, returned only if `select!=False`.
    """

    # Set logger level
    if debug:
        logger.setLevel('DEBUG')

    # Initialize the extractor
    extractor = SourceExtractor(algorithm=algorithm, fwhm=fwhm, sigma=noise_threshold)
    extractor.initialize_image(image, extension=None, var=image_var, dtype=cast_dtype, collapse=collapse)
    extractor.initialize_star_finder()

    # Reset parameters if requested
    if not background_subtraction:
        extractor.image.sky_bkg = 0.0
        logger.debug(f"Resetting sky background to {extractor.image.sky_bkg}")

    # Find sources
    sources = extractor.find_sources(uncertainties=True)

    # Save sources table to file, if requested
    if write_to is not None:
        extractor.write_to(write_to)

    if show:
        extractor.show()

    if select:
        # Interpret select type
        message = None
        save_to = None
        if isinstance(select, dict):
            message = select.get('message', 'Please select sources!')
            save_to = select.get('save_to', None)
        elif isinstance(select, str):
            message = 'Please select sources!'
            save_to = select

        # Select sources
        selected = extractor.select(save_to=save_to, message=message, radius=matching_threshold)
        logger.info(f"Selected {len(selected)} sources")
        logger.debug(f"\n{selected}")
        return sources, selected

    return sources
