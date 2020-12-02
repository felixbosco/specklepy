import numpy as np
import os

from astropy.io import fits

from specklepy.logging import logger


class MasterDark(object):

    extensions = {'variance': 'VAR', 'mask': 'MASK'}

    def __init__(self, file_list, file_name='MasterDark.fits', file_path=None, out_dir=None, new=True):
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
        self.means = []

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

    def write(self, overwrite=True):

        # Build primary HDU
        header = fits.Header()
        for index, file in enumerate(self.files):
            header.set(f"HIERARCH SPECKLEPY SOURCE FILE{index:04} NAME", os.path.basename(file))
        primary = fits.PrimaryHDU(data=self.image, header=header)

        # Build HDU list
        hdu_list = fits.HDUList([primary])

        # Build variance HDU
        if self.var is not None:
            var_hdu = fits.ImageHDU(data=self.var, name=self.extensions.get('variance'))
            hdu_list.append(var_hdu)

        # Build mask HDU
        if self.mask is not None:
            mask_hdu = fits.ImageHDU(data=self.mask, name=self.extensions.get('mask'))
            hdu_list.append(mask_hdu)

        # Write HDU list to file
        logger.info(f"Writing master dark frame to file {self.file_name!r}")
        hdu_list.writeto(self.file_name, overwrite=overwrite)
