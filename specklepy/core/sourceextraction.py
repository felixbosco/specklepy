from IPython import embed
import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from photutils import DAOStarFinder, IRAFStarFinder
from photutils import CircularAperture

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.logging import logger
from specklepy.plotting.plot import StarFinderPlot
from specklepy.utils import save_eval


def extract_sources(image, noise_threshold, fwhm, star_finder='DAO', image_var=None, background_subtraction=True,
                    show=False, select=False, write_to=None, cast_dtype=None, debug=False):
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
        show (bool, optional):
            Create a plot of the identified sources on the image.
        select (bool, optional):
            Create a plot of the identified sources on the image, with the option of selecting apertures for future
            purposes. These are selected by clicking on the image.
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
            List of selected sources, returned only if `select==True`.
    """

    # Set logger level
    if debug:
        logger.setLevel('DEBUG')

    # Input parameters
    if isinstance(image, np.ndarray):
        filename = 'current cube'
    elif isinstance(image, str):
        logger.info(f"The argument image '{image!r}' is interpreted as file name.")
        filename = image
        image = fits.getdata(filename)
        logger.debug(f"Data type of file input is {image.dtype}")
        image = image.squeeze()
    else:
        raise SpecklepyTypeError('extract_sources()', argname='image', argtype=type(image),
                                 expected='np.ndarray or str')

    # Cast the image data type if requested
    if cast_dtype is not None:
        image = image.astype(save_eval(cast_dtype))

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
    logger.debug(f"Extraction parameters:\n\tFWHM = {fwhm}\n\tThreshold = {threshold}\n\tSky = {sky}")
    sources = star_finder(image)

    # Reformatting sources table
    sources.sort('flux', reverse=True)
    sources.rename_column('xcentroid', 'x')
    sources.rename_column('ycentroid', 'y')
    sources.keep_columns(['x', 'y', 'flux'])

    # Add terminal output
    logger.info(f"Extracted {len(sources)} sources")
    logger.debug(f"\n{sources}")

    # Save sources table to file, if requested
    if write_to is not None:
        logger.info(f"Writing list of sources to file {write_to!r}")
        sources.write(write_to, format='ascii.fixed_width', overwrite=True)

    if show:
        plot = StarFinderPlot(image_data=image)
        positions = np.transpose((sources['x'], sources['y']))
        plot.add_apertures(positions=positions, radius=fwhm/2)
        plot.show()

    if select:
        plot = StarFinderPlot(image_data=image)
        positions = np.transpose((sources['x'], sources['y']))
        plot.add_apertures(positions=positions, radius=fwhm/2)
        selected = plot.select_apertures(marker_size=100)
        return sources, selected

    return sources


class SourceExtractor(object):

    def __init__(self, fwhm, algorithm='DAO', sigma=5.0):

        # Store input
        self.fwhm = fwhm
        self.sigma = sigma
        self.algorithm = algorithm

        # Initialize attributes
        self.star_finder = None
        self.image = None
        self.sources = None

    def __call__(self, source, extension=None, dtype=None):
        self.initialize_image(source, extension=extension, dtype=dtype)
        self.initialize_star_finder()
        return self.find_sources()

    @property
    def threshold(self):
        return self.sigma * self.image.stddev

    def initialize_image(self, source, extension=None, dtype=None):
        if isinstance(source, str):
            self.image = StarFinderImage.from_file(source, extension=extension)
        else:
            self.image = StarFinderImage(source)

        # Cast data type
        if dtype is not None:
            self.image.data = self.image.data.astype(save_eval(dtype))

        # Initialize statistics
        self.image.sigma_clipped_statistics()

    def initialize_star_finder(self):

        # Build parameter dictionary
        params = {'fwhm': self.fwhm,'threshold': self.threshold, 'sky': self.image.sky_bkg}

        # Type and value check on algorithm
        if not isinstance(self.algorithm, str):
            raise SpecklepyTypeError('extract_sources', argname='starfinder', argtype=type(self.algorithm),
                                     expected='str')
        if 'dao' in self.algorithm.lower():
            self.star_finder = DAOStarFinder(**params)
        elif 'iraf' in self.algorithm.lower():
            self.star_finder = IRAFStarFinder(**params)
        else:
            raise SpecklepyValueError('extract_sources', argname='star_finder', argvalue=self.algorithm,
                                      expected="'DAO' or 'IRAF")

        return self.star_finder

    def find_sources(self):
        # Find stars
        logger.info("Extracting sources...")
        logger.debug(f"Extraction parameters:"
                     f"\n\tFWHM = {self.fwhm}"
                     f"\n\tThreshold = {self.threshold}"
                     f"\n\tSky = {self.image.sky_bkg}")
        sources = self.star_finder(self.image.data)

        # Reformatting sources table
        sources.sort('flux', reverse=True)
        sources.rename_column('xcentroid', 'x')
        sources.rename_column('ycentroid', 'y')
        sources.keep_columns(['x', 'y', 'flux'])

        # Add terminal output
        logger.info(f"Extracted {len(sources)} sources")
        logger.debug(f"\n{sources}")

        # Store source table
        self.sources = sources
        return sources


class StarFinderImage(object):

    def __init__(self, image, filename=None):
        self.data = image
        self.filename = filename
        self._stddev = None
        self._sky_bkg = None

    @property
    def stddev(self):
        if self._stddev is None:
            self.sigma_clipped_statistics()
        return self._stddev

    @property
    def sky_bkg(self):
        if self._sky_bkg is None:
            self.sigma_clipped_statistics()
        return self._sky_bkg

    @classmethod
    def from_file(cls, filename, extension=None):
        logger.info(f"Read FITS image from file {filename!r} [{str(extension)}]")
        image = fits.getdata(filename, extension)
        logger.debug(f"Data type of file input is {image.dtype}")
        image = image.squeeze()
        return cls(image=image, filename=filename)

    def sigma_clipped_statistics(self, sigma=3.0):
        mean, median, std = sigma_clipped_stats(data=self.data, sigma=sigma)
        logger.info(f"Noise statistics for {self.filename!r}:"
                    f"\n\tMean = {mean:.3}\n\tMedian = {median:.3}\n\tStdDev = {std:.3}")
        self._sky_bkg = mean
        self._stddev = std
        return mean, std
