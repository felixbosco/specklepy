import numpy as np
import os

from scipy import signal

from astropy.stats import sigma_clipped_stats
from astropy.table import Table

from photutils import DAOStarFinder, IRAFStarFinder, find_peaks

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.fits import get_data
from specklepy.logging import logger
from specklepy.plotting.plot import StarFinderPlot
from specklepy.utils import save_eval
from specklepy.utils.array import peak_index
from specklepy.utils.box import Box
from specklepy.utils.moment import moment_2d


def extract_sources(image, noise_threshold, fwhm, star_finder='DAO', image_var=None, background_subtraction=True,
                    show=False, select=False, collapse=False, write_to=None, cast_dtype=None, debug=False):
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
        select (bool or dict, optional):
            Create a plot of the identified sources on the image, with the option of selecting apertures for future
            purposes. These are selected by clicking on the image.
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
    extractor = SourceExtractor(algorithm=star_finder, fwhm=fwhm, sigma=noise_threshold)
    extractor.initialize_image(image, extension=None, dtype=cast_dtype, collapse=collapse)
    extractor.initialize_star_finder()

    # Reset parameters if requested
    if image_var is not None:
        if isinstance(image_var, str):
            logger.info(f"Reading data from {extractor.image.filename!r} [{image_var}]")
            image_var = get_data(extractor.image.filename, extension=image_var)
        extractor.stddev = np.sqrt(np.mean(image_var))
        logger.debug(f"Resetting image standard deviation to {extractor.stddev}")
    if not background_subtraction:
        extractor.sky_bkg = 0.0
        logger.debug(f"Resetting sky background to {extractor.sky_bkg}")

    # Find sources
    sources = extractor.find_sources(uncertainties=True, image_var=image_var)

    # Save sources table to file, if requested
    if write_to is not None:
        extractor.write_to(write_to)

    if show:
        extractor.show()

    if select:
        message = None
        save_to = None
        if isinstance(select, dict):
            message = select.get('message', 'Please select sources!')
            save_to = select.get('save_to', None)
        selected = extractor.select(save_to=save_to, message=message)
        logger.info(f"Selected {len(selected)} sources")
        logger.debug(f"\n{selected}")
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

    def __call__(self, source, extension=None, dtype=None, collapse=False):
        self.initialize_image(source, extension=extension, dtype=dtype, collapse=collapse)
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

    def initialize_image(self, source, extension=None, dtype=None, collapse=False):
        if isinstance(source, str):
            self.image = StarFinderImage.from_file(source, extension=extension)
        else:
            self.image = StarFinderImage(source)

        # Cast data type
        if dtype is not None:
            self.image.data = self.image.data.astype(save_eval(dtype))

        # Collapse data cube
        if collapse:
            self.image.collapse_image()

        # Initialize statistics
        self.image.sigma_clipped_statistics()

    def initialize_star_finder(self):

        # Build parameter dictionary
        params = {'fwhm': self.fwhm, 'threshold': self.threshold}  # , 'sky': self.image.sky_bkg}

        # Type and value check on algorithm
        if not isinstance(self.algorithm, str):
            raise SpecklepyTypeError('extract_sources', argname='starfinder', argtype=type(self.algorithm),
                                     expected='str')
        if 'dao' in self.algorithm.lower():
            self.star_finder = DAOStarFinder(**params)
        elif 'iraf' in self.algorithm.lower():
            self.star_finder = IRAFStarFinder(**params)
        elif 'peak' in self.algorithm.lower():
            self.star_finder = PeakFinder(**params)
        else:
            raise SpecklepyValueError('extract_sources', argname='star_finder', argvalue=self.algorithm,
                                      expected="'DAO', 'IRAF', or 'PeakFinder'")
        logger.debug(f"StarFinder algorithm is initialized as: {self.star_finder}")

        return self.star_finder

    def find_sources(self, uncertainties=False, image_var=None):
        # Find stars
        logger.info("Extracting sources...")
        logger.debug(f"Extraction parameters:"
                     f"\n\tFWHM = {self.fwhm}"
                     f"\n\tThreshold = {self.threshold}"
                     f"\n\tSky = {self.image.sky_bkg}")
        sources = self.star_finder(np.ma.masked_less(self.image.data - self.image.sky_bkg, 0).filled(0))

        # Reformatting sources table
        sources.sort('flux', reverse=True)
        sources.rename_column('xcentroid', 'x')
        sources.rename_column('ycentroid', 'y')
        sources.keep_columns(['x', 'y', 'flux'])

        # Re-measure positions and flux to obtain uncertainties
        if uncertainties:
            sources = self.estimate_uncertainties(sources=sources, image_var=image_var)

        # Add terminal output
        logger.info(f"Extracted {len(sources)} sources")
        logger.debug(f"\n{sources}")

        # Store source table
        self.sources = sources
        return sources

    def estimate_uncertainties(self, sources, image_var):
        new_sources = Table(names=['x', 'dx', 'y', 'dy', 'flux', 'dflux'])
        radius = round(self.fwhm / 2.35 * 1.5)  # 1.5 times the standard deviation
        logger.info(f"Measuring uncertainties over {(radius * 2 + 1)}x{radius * 2 + 1} pixels...")

        # Iterate through sources
        for source in sources:
            pos = int(round(source['y'])), int(round(source['x']))
            box = Box.centered_at(pos[1], pos[0], radius=radius)
            box.crop_to_shape(shape=self.image.data.shape)
            aperture = box(self.image.data)
            try:
                var = box(image_var)
            except TypeError:
                var = None

            # Compute moments and uncertainties, and add to new table
            try:
                f, df, x, dx, y, dy = moment_2d(aperture, var=var)
                new_sources.add_row([x + box.x_min, dx, y + box.y_min, dy, f, df])
            except ValueError:
                logger.warning(f"Unsuccessful to measure uncertainties for source in box {box}")
                new_sources.add_row([source['x'], 0, source['y'], 0, source['flux'], 0])

        return new_sources

    def show(self):
        plot = StarFinderPlot(image_data=self.image.data)
        if self.sources is not None:
            positions = np.transpose((self.sources['x'], self.sources['y']))
            plot.add_apertures(positions=positions, radius=self.fwhm / 2)
        plot.show()

    def select(self, save_to=None, message=None, cross_match=True, radius=None):
        if self.sources is None:
            logger.warning("SourceExtractor does not store identified sources yet! Finding sources now...")
            self.find_sources()

        # Create plot
        plot = StarFinderPlot(image_data=self.image.data)
        positions = np.transpose((self.sources['x'], self.sources['y']))
        plot.add_apertures(positions=positions, radius=self.fwhm / 2)

        # Prompt message
        if message is not None:
            logger.info(message)

        # Graphically select apertures
        positions = plot.select_apertures(marker_size=100)
        if cross_match:
            selected = self.cross_match(positions)
        else:
            selected = self.find_closest_peaks(positions, radius=radius)

        # Save results
        if save_to is not None and len(selected) > 0:
            selected.write(save_to, format='ascii.fixed_width', overwrite=True)

        return selected

    def cross_match(self, positions):
        indexes = []
        for pos in positions:
            distances = self.compute_distances(pos)
            indexes.append(np.argmin(distances))
        return self.sources[indexes]

    def find_closest_peaks(self, positions, radius):
        # Initialize table of identified peaks
        peaks_table = Table(names=['x', 'dx', 'y', 'dy', 'flux', 'dflux'])

        # Iterate through guesses and identify peaks within a radius from the guess
        for pos in positions:
            peak = self.find_closest_peak(image=self.image.data, guess=pos, radius=radius)
            peaks_table.add_row([peak[1], None, peak[0], None, None, None])

        # Estimate uncertainties
        peaks_table = self.estimate_uncertainties(peaks_table, image_var=np.square(self.image.var))

        return peaks_table

    @staticmethod
    def find_closest_peak(image, guess, radius):
        box = Box.centered_at(guess[1], guess[0], radius=radius)
        aperture = box(image)
        peak = peak_index(aperture)
        return peak[0] + box.y_min, peak[1] + box.x_min

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

    def cross_match_table(self, table):

        # Initialize array of summed distances
        minimum_distances_total = np.empty((len(self.sources)))

        # Reference pattern: Stellar positions relative to first star
        table_dx = table['x'] - table['x'][0]
        table_dy = table['y'] - table['y'][0]

        # Iterate through sources to identify reference star from `table`
        for j in range(len(self.sources)):
            # Stellar positions relative to source `j`
            sources_dx = self.sources['x'] - self.sources['x'][j]
            sources_dy = self.sources['y'] - self.sources['y'][j]

            # Compute distances
            for i, (xi, yi) in enumerate(zip(table_dx, table_dy)):
                ri = np.sqrt((xi - sources_dx) ** 2 + (yi - sources_dy) ** 2)
                j_min = np.argmin(ri)
                minimum_distances_total[j] += ri[j_min]

        # Estimate distance-minimizing counter part to reference star
        best_match = np.argmin(minimum_distances_total)

        # Estimate shift
        shift = table[0]['y'] - self.sources[best_match]['y'], table[0]['x'] - self.sources[best_match]['x']

        return shift


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
        # image = fits.getdata(filename, extension)
        image = get_data(filename, extension=extension)
        logger.debug(f"Data type of file input is {image.dtype}")
        image = image.squeeze()
        return cls(image=image, filename=filename)

    def collapse_image(self):
        self.data = np.sum(self.data, axis=0)

    def sigma_clipped_statistics(self, sigma=3.0):
        mean, median, std = sigma_clipped_stats(data=self.data, sigma=sigma)
        logger.info(f"Noise statistics for {self.filename!r}:"
                    f"\n\tMean = {mean:.3}\n\tMedian = {median:.3}\n\tStdDev = {std:.3}")
        self._sky_bkg = mean
        self._stddev = std
        return mean, std


class PeakFinder(object):

    def __init__(self, fwhm, threshold):
        self.fwhm = fwhm
        self.threshold = threshold

    def __call__(self, image, *args, **kwargs):
        # peaks = find_peaks_2d(array=image, height=self.threshold, width=self.fwhm)#, threshold=5, distance=1)
        # sources = Table(names=['xcentroid', 'ycentroid', 'flux'])
        # for peak in peaks:
        #     sources.add_row([peak[0], peak[1], image[peak[1], peak[0]]])
        sources = find_peaks(data=image, threshold=self.threshold, box_size=1.5*self.fwhm)
        sources.rename_column('x_peak', 'xcentroid')
        sources.rename_column('y_peak', 'ycentroid')
        sources.rename_column('peak_value', 'flux')
        return sources

    def __repr__(self):
        return f"PeakFinder(fwhm: {self.fwhm:.2f}, threshold: {self.threshold:.2f})"


def find_peaks_1d(array, axis=None, **kwargs):
    if axis == 1:
        array = array.transpose()
    peaks, props = signal.find_peaks(array.flatten(), **kwargs)
    peaks = np.unravel_index(peaks, array.shape)
    if axis == 1:
        peaks = peaks[1], peaks[0]

    return peaks


def find_peaks_2d(array, **kwargs):
    peaks_x = find_peaks_1d(array=array, axis=0, **kwargs)
    peaks_y = find_peaks_1d(array=array, axis=1, **kwargs)

    matches = []
    for pxx, pxy in zip(peaks_x[0], peaks_x[1]):
        good_x = pxx == peaks_y[0]
        good_y = pxy == peaks_y[1]
        matching = good_x & good_y
        if np.any(matching):
            matches.append([peaks_y[0][matching], peaks_y[1][matching]])

    matches = np.array(matches).transpose()[0]

    return matches
