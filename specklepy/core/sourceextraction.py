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
from specklepy.utils.point import Point
from specklepy.utils.vector import Vector


class SourceExtractor(object):

    def __init__(self, fwhm, algorithm='DAO', sigma=5.0):

        # Store input
        self.fwhm = fwhm
        self.sigma = sigma
        self.algorithm = algorithm

        # Initialize attributes
        self.image = None
        self.sources = None

    @property
    def threshold(self):
        return self.sigma * np.mean(self.image.stddev)

    def initialize_image(self, source, extension=None, var=None, dtype=None, collapse=False):
        if isinstance(source, str):
            self.image = StarFinderImage.from_file(source, extension=extension, var=var, dtype=save_eval(dtype))
        else:
            self.image = StarFinderImage(source, var=var)

        # Collapse data cube
        if collapse:
            self.image.collapse_image()

    def initialize_star_finder(self):

        # Build parameter dictionary
        params = {'fwhm': self.fwhm, 'threshold': self.threshold}  # , 'sky': self.image.sky_bkg}

        # Type and value check on algorithm
        if not isinstance(self.algorithm, str):
            raise SpecklepyTypeError('extract_sources', argname='starfinder', argtype=type(self.algorithm),
                                     expected='str')
        if 'dao' in self.algorithm.lower():
            star_finder = DAOStarFinder(**params)
        elif 'iraf' in self.algorithm.lower():
            star_finder = IRAFStarFinder(**params)
        elif 'peak' in self.algorithm.lower():
            star_finder = PeakFinder(**params)
        elif 'manu' in self.algorithm.lower():
            star_finder = ManualFinder(**params)
        else:
            raise SpecklepyValueError('extract_sources', argname='star_finder', argvalue=self.algorithm,
                                      expected="'DAO', 'IRAF', or 'PeakFinder'")
        logger.debug(f"StarFinder algorithm is initialized as: {star_finder}")

        return star_finder

    def find_sources(self, uncertainties=False):

        # Initialize StarFinder algorithm
        star_finder = self.initialize_star_finder()

        # Find stars
        logger.info("Extracting sources...")
        logger.debug(f"Extraction parameters:"
                     f"\n\tFWHM = {self.fwhm}"
                     f"\n\tThreshold = {self.threshold}"
                     f"\n\tSky = {self.image.sky_bkg}")
        sources = star_finder(np.ma.masked_less(self.image.data - self.image.sky_bkg, 0).filled(0))

        # Reformatting sources table
        sources.sort('flux', reverse=True)
        sources = self.rename_columns(sources)
        sources.keep_columns(['x', 'y', 'flux'])

        # Re-measure positions and flux to obtain uncertainties
        if uncertainties:
            sources = self.estimate_uncertainties(sources=sources)

        # Add terminal output
        logger.info(f"Extracted {len(sources)} sources")
        logger.debug(f"\n{sources}")

        # Store source table
        self.sources = sources
        return sources

    @staticmethod
    def rename_columns(table):
        aliases = {'xcentroid': 'x', 'ycentroid': 'y', 'x_peak': 'x', 'y_peak': 'y', 'peak_value': 'flux'}
        for key, alias in aliases.items():
            if key in table.colnames:
                table.rename_column(key, alias)
        return table

    def estimate_uncertainties(self, sources, weighting='gaussian'):
        """Estimate the astrometric and photometric uncertainties of the sources in the provided table.

        This function uses a two-dimensional moment with propagated uncertainties to re-measure the flux and position
        values of the sources as the zeroth and first moment within an aperture, centered on the provided position of
        the sourcse.

        Args:
            sources (astropy.table.Table):
                Table of sources.
            weighting (str, optional):
                Weighting scheme. Can be 'gaussian' (default), 'uniform', or a weighting kernel.

        Returns:
            sources (astropy.table.Table):
                Table of sources with updated positions and flux values, plus corresponding uncertainties
        """

        # Initialize output table
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
                var = box(np.square(self.image.stddev))
            except TypeError:
                var = None
            except IndexError:
                var = np.full(aperture.shape, np.square(self.image.stddev))

            # Compute moments and uncertainties, and add to new table
            try:
                f, df, x, dx, y, dy = moment_2d(aperture, var=var, weights=weighting, fwhm=self.fwhm)
                if box.x_min is not None:
                    x += box.x_min
                if box.y_min is not None:
                    y += box.y_min
                new_sources.add_row([x, dx, y, dy, f, df])
            except ValueError:
                logger.warning(f"Unsuccessful to measure uncertainties for source in box {box}")
                new_sources.add_row([source['x'], 0, source['y'], 0, source['flux'], 0])

        return new_sources

    def show(self):
        plot = StarFinderPlot(image_data=self.image.data)
        if self.sources is not None:
            positions = np.transpose((self.sources['y'], self.sources['x']))
            plot.add_apertures(positions=positions, radius=self.fwhm / 2)
        plot.show()

    def select(self, save_to=None, message=None, radius=None):
        if self.sources is None:
            logger.warning("SourceExtractor does not store identified sources yet! Finding sources now...")
            self.find_sources()

        # Create plot
        plot = StarFinderPlot(image_data=self.image.data)
        positions = np.transpose((self.sources['y'], self.sources['x']))
        plot.add_apertures(positions=positions, radius=self.fwhm / 2)

        # Prompt message
        if message is not None:
            logger.info(message)

        # Graphically select apertures
        positions = plot.select_apertures(marker_size=100)
        selected = self.cross_match(positions, radius_threshold=radius)

        # Save results
        if save_to is not None and len(selected) > 0:
            selected.write(save_to, format='ascii.fixed_width', overwrite=True)

        return selected

    def cross_match(self, guesses, radius_threshold=None):
        indexes = []
        for guess in guesses:
            distances = self.compute_distances(guess)
            try:
                if np.min(distances) > radius_threshold:
                    # Reject source to far away from known sources
                    continue
            except TypeError:
                pass
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
        dx = position[1] - self.sources['x'].data
        dy = position[0] - self.sources['y'].data
        return np.sqrt(np.add(np.square(dx), np.square(dy)))

    def cross_match_table(self, table, threshold=5, penalty=100):
        """Cross match the `sources` table with a reference `table` in order to measure the shift between two images.

        The function iterates through the table of identified `sources`, counts how many stars have a counterpart if
        the current star is the brightest of the reference stars. The best match then serves as a measure to estimate
        the shift between the two sets of data.s

        Args:
            table (astropy.table.Table):
                Table of source to which the table of sources shall be cross-matched.
            threshold (int, optional):
                Threshold in offset pixels, serving as a measure whether a source has a counterpart or not.
            penalty:
                Penalty that individual combinations receive for one source not having a counterpart.

        Returns:
            shift (tuple):
                Shift between the sources in the `sources` table and the comparison `table`, after identification of the
                corresponding sources.
        """

        # Initialize array of summed distances
        minimum_distances_total = np.empty((len(self.sources)))

        # Reference pattern: Stellar positions relative to brightest star
        table.sort('flux', reverse=True)
        reference = Vector(table[0]['y'], table[0]['x'])
        table_dx = table['x'] - reference[1]
        table_dy = table['y'] - reference[0]

        # Iterate through sources to identify reference star from `table`
        for j in range(len(self.sources)):
            # Stellar positions relative to source `j`
            sources_dx = self.sources['x'] - self.sources['x'][j]
            sources_dy = self.sources['y'] - self.sources['y'][j]

            # Compute distances
            for i, (dxi, dyi) in enumerate(zip(table_dx, table_dy)):
                ri = np.sqrt(np.square(dxi - sources_dx) + np.square(dyi - sources_dy))
                j_min = np.argmin(ri)
                if ri[j_min] < threshold:
                    minimum_distances_total[j] += ri[j_min]
                else:
                    minimum_distances_total[j] += penalty

        # Estimate distance-minimizing counter part to reference star
        best_match = np.argmin(minimum_distances_total)
        best_match = Vector(self.sources[best_match]['y'], self.sources[best_match]['x'])

        # Estimate shift
        shift = (reference - best_match).astype(tuple)

        return shift


class StarFinderImage(object):

    def __init__(self, image, filename=None, var=None):
        self.data = image
        self.filename = filename
        self._stddev = np.sqrt(var) if var is not None else None
        self._sky_bkg = None

    @property
    def stddev(self):
        if self._stddev is None:
            _, self._stddev = self.sigma_clipped_statistics()
        return self._stddev

    @stddev.setter
    def stddev(self, value):
        self._stddev = value

    @property
    def var(self):
        return np.square(self.stddev)

    @property
    def sky_bkg(self):
        if self._sky_bkg is None:
            self._sky_bkg, _ = self.sigma_clipped_statistics()
        return self._sky_bkg

    @sky_bkg.setter
    def sky_bkg(self, value):
        self._sky_bkg = value

    @classmethod
    def from_file(cls, filename, extension=None, var=None, dtype=None):
        logger.info(f"Read FITS image from file {filename!r} [{str(extension)}]")
        image = get_data(filename, extension=extension, dtype=dtype, squeeze=True)
        logger.debug(f"Data type of file input is {image.dtype}")
        if isinstance(var, str):
            var = get_data(filename, extension=var, dtype=dtype, squeeze=True)
        return cls(image=image, filename=filename, var=var)

    def collapse_image(self):
        self.data = np.sum(self.data, axis=0)

    def sigma_clipped_statistics(self, sigma=3.0):
        mean, median, std = sigma_clipped_stats(data=self.data, sigma=sigma)
        logger.info(f"Noise statistics for {self.filename!r}:"
                    f"\n\tMean = {mean:.3}\n\tMedian = {median:.3}\n\tStdDev = {std:.3}")
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


class ManualFinder(object):

    def __init__(self, radius=None, **kwargs):

        self.radius = radius
        if self.radius is None:
            self.radius = kwargs.get('fwhm', 20)

    def __call__(self, image):

        # Initialize table of identified peaks
        peaks_table = Table(names=['x', 'dx', 'y', 'dy', 'flux', 'dflux'])

        # Create plot and graphically select apertures
        plot = StarFinderPlot(image_data=image)
        logger.info("Select guesses for iterative source selection!")
        guesses = plot.select_apertures(marker_size=100)

        # Iterate through guesses and identify peaks within a radius from the guess
        for guess in guesses:
            peak = self.find_closest_peak(image=image, guess=guess, radius=self.radius)
            peaks_table.add_row([peak[1], None, peak[0], None, image[peak], None])

        return peaks_table

    @staticmethod
    def find_closest_peak(image, guess, radius, threshold=0.1):
        guess = Point(*guess, order='yx')
        peak = None

        done = False
        while not done:
            box = Box.centered_at(x0=guess.x, y0=guess.y, radius=radius)
            box.crop_to_shape(image.shape)
            aperture = box(image)
            peak = peak_index(aperture)
            peak = Point(*peak, order='yx')

            # Shift into full frame reference
            if box.y_min is not None:
                peak.y += box.y_min
            if box.x_min is not None:
                peak.x += box.x_min

            # Test for convergence
            if peak.distance_to(guess) < threshold:
                done = True
            else:
                guess = peak
        return peak.y, peak.x
