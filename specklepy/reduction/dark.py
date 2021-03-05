import numpy as np
import os

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from specklepy.logging import logger
from specklepy.utils.time import default_time_stamp


class MasterDark(object):

    extensions = {'variance': 'VAR', 'mask': 'MASK'}

    def __init__(self, file_list, file_name='MasterDark.fits', file_path=None, out_dir=None, setup=None, new=True):
        self.files = file_list
        self.file_name = self.insert_setup_to_file_name(file_name=file_name, setup=setup)
        self.file_path = file_path if file_path is not None else ''
        self.out_dir = out_dir if out_dir is not None else ''

        # Initialize maps
        self.means = None
        self.image = None
        self.var = None
        self.mask = None

    @classmethod
    def from_file(cls, file_path):
        # Create object from path information
        out_dir, file_name = os.path.split(file_path)
        obj = cls(file_list=None, file_name=file_name, out_dir=out_dir, setup=None)

        # Load data from file
        obj.image = fits.getdata(obj.path)
        try:
            obj.var = fits.getdata(obj.path, obj.extensions.get('variance'))
        except KeyError:
            logger.debug(f"Loading MasterDark from file {obj.path!r} without {obj.extensions.get('variance')!r} "
                         f"extension")
        try:
            obj.mask = fits.getdata(obj.path, obj.extensions.get('mask')).astype(bool)
        except KeyError:
            logger.debug(f"Loading MasterDark from file {obj.path!r} without {obj.extensions.get('mask')!r} "
                         f"extension")

        return obj

    @property
    def path(self):
        return os.path.join(self.out_dir, self.file_name)

    @staticmethod
    def insert_setup_to_file_name(file_name, setup=None):
        if setup is None:
            return file_name
        else:
            base, ext = os.path.splitext(file_name)
            return f"{base}_{setup}{ext}"

    def combine(self, number_frames=None):
        logger.info("Combining master dark frame...")
        if number_frames is not None:
            logger.debug(f"Using only the first {number_frames} frames of each cube")
        self.means = []

        # Iterate through files
        for file in self.files:
            logger.info(f"Reading dark frames from file {file!r}...")
            path = os.path.join(self.file_path, file)
            with fits.open(path) as hdu_list:
                cube = hdu_list[0].data
                logger.info("Computing sigma-slipped statistics of data cube...")
                mean, _, std = sigma_clipped_stats(data=cube[:number_frames], axis=0)

                self.means.append(mean)
                self.combine_var(np.square(std))
                self.combine_mask(std == 0)

        self.image = np.mean(np.array(self.means), axis=0)
        self.var /= np.square(len(self.files))  # For correctly propagating the uncertainties

    def combine_var(self, new_var):
        if self.var is None:
            self.var = new_var
        else:
            self.var = np.add(self.var, new_var)

    def combine_mask(self, new_mask):
        if self.mask is None:
            self.mask = new_mask
        else:
            self.mask = np.logical_or(self.mask, new_mask)

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
            mask_hdu = fits.ImageHDU(data=self.mask.astype(np.int16), name=self.extensions.get('mask'))
            hdu_list.append(mask_hdu)

        # Write HDU list to file
        logger.info(f"Writing master dark frame to file {self.path!r}")
        hdu_list.writeto(self.path, overwrite=overwrite)

    def subtract(self, file_path, extension=None):
        """Subtract the master dark from a file containing image data.

        The master dark is subtracted from the image or each frame in a data cube. Then uncertainties are propagated.

        Arguments:
            file_path (str):
                Path to the file, containing image data.
            extension (str, optional):
                Classifier for the image data extension.
        """

        logger.info(f"Subtracting master dark {self.file_name} from file at {file_path!r}")

        # Load image data
        data = fits.getdata(file_path, extension)

        # Subtract
        if data.ndim == 2:
            data = data - self.image
        elif data.ndim == 3:
            for f, frame in enumerate(data):
                data[f] = frame - self.image

        # Propagate variances
        try:
            var = fits.getdata(file_path, self.extensions.get('variance'))
            has_var_hdu = True
            var += self.var
        except KeyError:
            has_var_hdu = False
            var = self.var

        # Propagate mask
        try:
            mask = fits.getdata(file_path, self.extensions.get('mask')).astype(bool)
            has_mask_hdu = True
            mask = np.logical_or(mask, self.mask)
        except KeyError:
            has_mask_hdu = False
            mask = self.mask

        # Store data to cube
        with fits.open(file_path, mode='update') as hdu_list:
            # Update header
            hdu_list[0].header.set('HIERARCH SPECKLEPY REDUCTION DARKCORR', default_time_stamp())

            # Image data
            hdu_list[0].data = data

            # Variance data
            if has_var_hdu:
                hdu_list[self.extensions.get('variance')] = var
            else:
                var_hdu = fits.ImageHDU(data=var, name=self.extensions.get('variance'))
                hdu_list.append(var_hdu)

            # Mask data
            if has_mask_hdu:
                hdu_list[self.extensions.get('mask')] = mask.as
type(np.int16)
            else:
                mask_hdu = fits.ImageHDU(data=mask.astype(np.int16), name=self.extensions.get('mask'))
                hdu_list.append(mask_hdu)

            # Write HDU list to file
            logger.info(f"Updating dark subtraction in file {file_path!r}")
            hdu_list.flush()
