from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder

from holopy.logging import logging


class SourceExtraction(object):

    def __init__(self):
        pass


    def __call__(self, *args, **kwargs):
        self.find_sources(*args, **kwargs)


    def find_sources(self, image, starfinder_fwhm, noise_threshold, background_subtraction=True, verbose=True):
        # Prepare noise statistics
        if isinstance(image, str):
            logging.info("The argument image '{}' is interpreted as file name.".format(image))
            image = fits.getdata(image)
        mean, median, std = sigma_clipped_stats(image, sigma=3.0)
        logging.info("Noise statistics for {}:\n\tMean = {:.3}\n\tMedian = {:.3}\n\tStdDev = {:.3}".format(image, mean, median, std))

        # Find stars
        daofind = DAOStarFinder(fwhm=starfinder_fwhm, threshold=noise_threshold*std)
        if background_subtraction:
            logging.info("Subtracting background...")
            image -= median
        logging.info("Finding sources...")
        self.sources = daofind(image)
        self.sources.sort('flux', reverse=True)
        self.sources.rename_column('xcentroid', 'x')
        self.sources.rename_column('ycentroid', 'y')
        self.sources.keep_columns(['x', 'y', 'flux'])
        for col in self.sources.colnames:
            self.sources[col].info.format = '%.8g'  # for consistent table output
        if verbose:
            logging.info("Found {} sources:\n{}".format(len(self.sources), self.sources))


    def writeto(self, file):
        logging.info("Writing list of sources to file {}".format(file))
        self.sources.write(file, format='ascii', overwrite=True)
