import numpy as np
import os

from astropy.io import fits

from specklepy.logging import logger


class MasterDark(object):

    def __init__(self, file_list, file_name='MasterFlat.fits', file_path=None, out_dir=None, new=True):
        self.files = file_list
        self.file_name = file_name
        self.file_path = file_path
        self.out_dir = out_dir

        # Initialize maps
        self.means = None
        self.image = None
        self.var = None
        self.mask = None

    def combine(self):
        logger.info("Combining master dark frame...")

        # Iterate through files
        for file in self.files:
            logger.info(f"Reading dark frames from file {file!r}...")
            path = os.path.join(self.file_path, file)
            with fits.open(path) as hdu_list:
                cube = hdu_list[0].data
                self.means.append(np.mean(cube, axis=0))
                var = cube.var(axis=0)
                self.combine_var(var)
                self.combine_mask(var == 0)

        self.image = np.mean(np.array(self.means), axis=0)

    def combine_var(self, new_var):
        if self.var is None:
            self.var = new_var
        else:
            self.var = np.add(self.var, new_var)

    def combine_mask(self, new_mask):
        if self.mask is None:
            self.mask = new_mask
        else:
            np.logical_or(self.mask, new_mask)

    def write(self):
        pass
