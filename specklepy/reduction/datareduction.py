from IPython import embed
import os

from specklepy.io.config import read
from specklepy.io.filearchive import ReductionFileArchive
from specklepy.logging import logger
from specklepy.reduction import dark


class DataReduction(object):

    config_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '../config/reduction.cfg'))

    def __init__(self, paths=None, dark=None, flat=None, sky=None, options=None, **kwargs):
        self.paths = paths
        self.dark = dark
        self.flat = flat
        self.sky = sky
        self.options = options
        for attr in ['paths', 'dark', 'flat', 'sky', 'options']:
            if getattr(self, attr) is None:
                logger.debug(f"Transforming parameters from section {attr.upper()!r} in the config file to {attr!r}...")
                setattr(self, attr, kwargs.get(attr.upper()))

        # Initialize file archive
        self.files = ReductionFileArchive(file_list=self.paths.get('fileList'), in_dir=self.paths.get('filePath'),
                                          out_dir=self.paths.get('outDir'))

    @classmethod
    def from_file(cls, file_name):
        logger.info(f"Configuring data reduction from config file {cls.config_file!r}")
        config = read(par_file=cls.config_file)
        logger.info(f"Updating data reduction configuration from parameter file {file_name!r}")
        params = read(par_file=file_name)
        for section in params:
            config[section].update(**params[section])
        return cls(**config)

    def run(self):
        self.initialize_directories()
        # self.run_post_correlation()
        self.run_dark_correction()
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

    def set_up(self):
        pass

    def run_post_correlation(self):
        pass

    def run_dark_correction(self):
        darks = self.files.filter({'OBSTYPE': 'DARK'})
        master_dark = dark.MasterDark(file_list=darks, file_path=self.files.in_dir,
                                      file_name=self.dark.get('masterDarkFile'), out_dir=self.paths.get('tmpDir'))
        master_dark.combine(number_frames=self.dark.get('numberFrames'))
        master_dark.write()

    def run_flat_fielding(self):
        pass

    def run_linearization(self):
        pass

    def run_sky_subtraction(self):
        pass
