import numpy as np
import os

from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats

from specklepy.logging import logger
from specklepy.reduction.subwindow import SubWindow
from specklepy.utils.time import default_time_stamp


class MasterDark(object):

    extensions = {'variance': 'VAR', 'mask': 'MASK'}

    def __init__(self, file_list, file_name='MasterDark.fits', file_path=None, out_dir=None, setup=None,
                 sub_window=None, new=True):
        self.files = file_list
        self.file_name = self.insert_setup_to_file_name(file_name=file_name, setup=setup)
        self.file_path = file_path if file_path is not None else ''
        self.out_dir = out_dir if out_dir is not None else ''

        # Store sub-window
        if isinstance(sub_window, str):
            self.sub_window = sub_window
        else:
            self.sub_window = np.unique(sub_window)[0]

        # Initialize maps
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
        obj.sub_window = fits.getheader(obj.path)["HIERARCH SPECKLEPY REDUCTION SUBWIN"]

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

    def combine(self, max_number_frames=None, rejection_threshold=10):
        logger.info("Combining master dark frame...")
        # if max_number_frames is not None:
        #     logger.debug(f"Using only the first {max_number_frames} frames of each cube")
        means = []
        vars = []
        number_frames = []

        # Iterate through files
        for file in self.files:
            logger.info(f"Reading DARK frames from file {file!r}...")
            path = os.path.join(self.file_path, file)
            with fits.open(path) as hdu_list:
                data = hdu_list[0].data.squeeze()

                if data.ndim == 2:
                    means.append(data)
                    vars.append(np.zeros(data.shape))
                    # self.combine_mask(np.zeros(data.shape, dtype=bool))
                    number_frames.append(1)

                elif data.ndim == 3:
                    logger.info("Computing statistics of data cube...")
                    mean = np.mean(data, axis=0)
                    std = np.std(data, axis=0)

                    # Identify outliers based on sigma-clipping
                    mean_mask = sigma_clip(mean, sigma=rejection_threshold, masked=True).mask
                    std_mask = sigma_clip(std, sigma=rejection_threshold, masked=True).mask
                    mask = np.logical_or(mean_mask, std_mask)
                    mask_indexes = np.array(np.where(mask)).transpose()

                    # Re-compute the identified pixels
                    logger.info(f"Re-measuring {len(mask_indexes)} outliers...")
                    for mask_index in mask_indexes:
                        # Extract t-series for the masked pixel
                        arr = data[:, mask_index[0], mask_index[1]]

                        # Compute sigma-clipped statistics for this pixel
                        arr_mean, _, arr_std = sigma_clipped_stats(arr, sigma=rejection_threshold)
                        mean[mask_index[0], mask_index[1]] = arr_mean
                        std[mask_index[0], mask_index[1]] = arr_std

                    mean = sigma_clip(mean, sigma=rejection_threshold, masked=True)
                    std = sigma_clip(std, sigma=rejection_threshold, masked=True)

                    # Store results into lists
                    means.append(mean)
                    vars.append(np.square(std))
                    # self.combine_mask(np.logical_or(mean.mask, std.mask))
                    number_frames.append(data.shape[0])

                else:
                    raise ValueError(f"Shape of data {data.shape} is not understood. Data must be either 2 or "
                                     f"3-dimensional!")

        # Combine data and propagate the uncertainties
        self.image = np.average(np.ma.masked_array(means), axis=0, weights=number_frames)
        self.var = np.average(np.ma.masked_array(vars), axis=0, weights=number_frames)

        # Combine mask
        self.mask = np.logical_or(self.image.mask, self.var.mask)
        self.image = self.image.data
        self.var = self.var.data

    # def combine_var(self, new_var):
    #     if self.var is None:
    #         self.var = new_var
    #     else:
    #         self.var = np.add(self.var, new_var)
    #
    # def combine_mask(self, new_mask):
    #     if self.mask is None:
    #         self.mask = new_mask
    #     else:
    #         self.mask = np.logical_or(self.mask, new_mask)

    def write(self, overwrite=True):

        # Build primary HDU
        header = fits.Header()
        for index, file in enumerate(self.files):
            header.set(f"HIERARCH SPECKLEPY SOURCE FILE{index:04} NAME", os.path.basename(file))
        header.set("HIERARCH SPECKLEPY REDUCTION SUBWIN", self.sub_window)
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

    def subtract(self, file_path, extension=None, sub_window=None):
        """Subtract the master dark from a file containing image data.

        The master dark is subtracted from the image or each frame in a data cube. Then uncertainties are propagated.

        Arguments:
            file_path (str):
                Path to the file, containing image data.
            extension (str, optional):
                Classifier for the image data extension.
        """

        logger.info(f"Subtracting master dark {self.file_name!r} from file at {file_path!r}")

        # Construct sub-window
        sub_window = SubWindow.from_str(sub_window, full=self.sub_window)

        # Construct good pixel mask
        if self.mask is None:
            gpm = np.ones(sub_window(self.image).shape, dtype=bool)
        else:
            gpm = ~self.mask

        # Load image data
        data = fits.getdata(file_path, extension)

        # Subtract
        if data.ndim == 2:
            data = np.subtract(data, sub_window(self.image), where=gpm)
        elif data.ndim == 3:
            for f, frame in enumerate(data):
                data[f] = np.subtract(frame, sub_window(self.image), where=gpm)

        # Propagate variances
        try:
            var = fits.getdata(file_path, self.extensions.get('variance'))
            has_var_hdu = True
            var = np.add(var, sub_window(self.var), where=gpm)
        except KeyError:
            has_var_hdu = False
            var = sub_window(self.var)

        # Propagate mask
        try:
            mask = fits.getdata(file_path, self.extensions.get('mask')).astype(bool)
            has_mask_hdu = True
            mask = np.logical_or(mask, sub_window(self.mask))
        except KeyError:
            has_mask_hdu = False
            mask = sub_window(self.mask)

        # Store data to cube
        with fits.open(file_path, mode='update') as hdu_list:
            # Update header
            hdu_list[0].header.set('HIERARCH SPECKLEPY REDUCTION DARKCORR', default_time_stamp())

            # Image data
            hdu_list[0].data = data

            # Variance data
            if has_var_hdu:
                hdu_list[self.extensions.get('variance')].data = var
            else:
                var_hdu = fits.ImageHDU(data=var, name=self.extensions.get('variance'))
                hdu_list.append(var_hdu)

            # Mask data
            if has_mask_hdu:
                hdu_list[self.extensions.get('mask')].data = mask.astype(np.int16)
            else:
                mask_hdu = fits.ImageHDU(data=mask.astype(np.int16), name=self.extensions.get('mask'))
                hdu_list.append(mask_hdu)

            # Write HDU list to file
            logger.info(f"Updating dark subtraction in file {file_path!r}")
            hdu_list.flush()
