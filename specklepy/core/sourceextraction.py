from IPython import embed
import numpy as np
import os

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from photutils import DAOStarFinder, IRAFStarFinder

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

    # Initialize the extractor
    extractor = SourceExtractor(algorithm=star_finder, fwhm=fwhm, sigma=noise_threshold)
    extractor.initialize_image(image, extension=None, dtype=cast_dtype)
    extractor.initialize_star_finder()

    # Reset parameters if requested
    if image_var is not None:
        if isinstance(image_var, str):
            logger.warning(f"'image_var' was provided as str-type. Interpreting as the name of a FITS extension "
                           f"containing the image variance.")
            logger.info(f"Reading data from {extractor.image.filename!r} [{image_var}]")
            image_var = fits.getdata(extractor.image.filename, image_var)
            image_var = np.sqrt(np.mean(image_var))
        extractor.stddev = np.sqrt(image_var)
        logger.debug(f"Resetting image standard deviation to {extractor.stddev}")
    if not background_subtraction:
        extractor.sky_bkg = 0.0
        logger.debug(f"Resetting sky background to {extractor.sky_bkg}")

    # Find sources
    sources = extractor.find_sources()

    # Save sources table to file, if requested
    if write_to is not None:
        extractor.write_to(write_to)

    if show:
        extractor.show()

    if select:
        save_selected_to = select if isinstance(select, str) else None
        selected = extractor.select(save_selected_to)
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

    @property
    def sky_bkg(self):
        return self.image.sky_bkg

    @sky_bkg.setter
    def sky_bkg(self, value):
        self.image.sky_bkg = value

    @property
    def stddev(self):
        return self.image.stddev

    @stddev.setter
    def stddev(self, value):
        self.image.std_dev = value

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

    def show(self):
        plot = StarFinderPlot(image_data=self.image.data)
        if self.sources is not None:
            positions = np.transpose((self.sources['x'], self.sources['y']))
            plot.add_apertures(positions=positions, radius=self.fwhm / 2)
        plot.show()

    def select(self, save_to=None):
        if self.sources is None:
            logger.warning("SourceExtractor does not store identified sources yet! Finding sources now...")
            self.find_sources()

        # Create plot
        plot = StarFinderPlot(image_data=self.image.data)
        positions = np.transpose((self.sources['x'], self.sources['y']))
        plot.add_apertures(positions=positions, radius=self.fwhm / 2)

        # Graphically select apertures
        selected = plot.select_apertures(marker_size=100)
        selected = self.cross_match(selected)

        # Save results
        if save_to is not None:
            selected.write(save_to, format='ascii.fixed_width', overwrite=True)
        return selected

    def cross_match(self, positions):
        indexes = []
        for pos in positions:
            distances = self.compute_distances(pos)
            indexes.append(np.argmin(distances))
        return self.sources[indexes]

    def write_to(self, filename):
        """Save source table to a file"""
        if filename == 'default':
            filename = self.default_file_name(self.image.filename)
        if self.sources is not None:
            logger.info(f"Writing list of sources to file {filename!r}")
            self.sources.write(filename, format='ascii.fixed_width', overwrite=True)

    @staticmethod
    def default_file_name(filename):
        return 'sources_' + os.path.basename(filename).replace('.fits', '.dat')

    def compute_distances(self, position):
        dx = position[0] - self.sources['x'].data
        dy = position[1] - self.sources['y'].data
        return np.sqrt(np.add(np.square(dx), np.square(dy)))


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

    @stddev.setter
    def stddev(self, value):
        self._stddev = value

    @property
    def sky_bkg(self):
        if self._sky_bkg is None:
            self.sigma_clipped_statistics()
        return self._sky_bkg

    @sky_bkg.setter
    def sky_bkg(self, value):
        self._sky_bkg = value

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
