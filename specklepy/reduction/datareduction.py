import os

from specklepy.logging import logger


class DataReduction(object):

    def __init__(self, **kwargs):
        self.paths = kwargs.get('PATHS')
        self.flat_fielding = kwargs.get('FLAT')
        self.sky_subtraction = kwargs.get('SKY')
        self.options = kwargs.get('OPTIONS')

    def run(self):
        self.run_post_correlation()
        self.run_dark_correction()
        self.initialize_directories()
        self.run_flat_fielding()
        self.run_linearization()
        self.run_sky_subtraction()

        # Close reduction
        logger.info("Reduction finished!")

    def initialize_directories(self):
        if not os.path.isdir(self.paths.get('outDir')):
            logger.debug(f"Making directory {self.paths.get('outDir')}")
            os.makedirs(self.paths.get('outDir'))
        if not os.path.isdir(self.paths.get('tmpDir')):
            logger.debug(f"Making directory {self.paths.get('tmpDir')}")
            os.makedirs(self.paths.get('tmpDir'))

    def run_post_correlation(self):
        pass

    def run_dark_correction(self):
        pass

    def run_flat_fielding(self):
        pass

    def run_linearization(self):
        pass

    def run_sky_subtraction(self):
        pass
