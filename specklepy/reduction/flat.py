from datetime import datetime
from IPython import embed
import numpy as np
import os
from tqdm import trange

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip, sigma_clipped_stats

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.masterfile import MasterFile
from specklepy.logging import logger
from specklepy.reduction.subwindow import SubWindow
from specklepy.utils.combine import combine_masks
from specklepy.utils.time import default_time_stamp


class MasterFlat(object):

    fill_value = 0  # was np.nan
    mask_threshold = 100

    extensions = {'variance': 'VAR', 'mask': 'MASK'}

    def __init__(self, file_list, file_name='MasterFlat.fits', file_path=None, out_dir=None, filter=None,
                 sub_window=None, new=True):
        """

        Arguments:
            sub_window (str, optional):
                String describing the full window in the same format as the sub-window strings for providing relative
                coordinates.
        """

        # Store input parameters
        if isinstance(file_list, (list, np.ndarray)) or file_list is None:
            self.files = file_list
        elif isinstance(file_list, Table):
            is_flat_file = file_list['OBSTYPE'] == 'FLAT'
            self.files = file_list['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_list', type(file_list), 'astropy.table.Table')

        if isinstance(file_name, str):
            self.file_name = self.insert_filter_to_file_name(file_name, filter=filter)
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_name', type(file_name), 'str')

        if isinstance(file_path, str) or file_path is None:
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_path', type(file_path), 'str')
        self.out_dir = out_dir

        # Store sub-window
        if isinstance(sub_window, str):
            self.sub_window = sub_window
        else:
            self.sub_window = np.unique(sub_window)[0]

        # Create an output file
        # self.master_file = MasterFile(self.file_name, files=self.files, in_dir=file_path, out_dir=out_dir,
        #                               initialize=new)

        # Initialize maps
        self.means = None
        self.image = None
        self.var = None
        self.mask = None

        # Initialize parameters
        self.normalized = False

    @classmethod
    def from_file(cls, file_path):
        # Create object from path information
        out_dir, file_name = os.path.split(file_path)
        obj = cls(file_list=None, file_name=file_name, out_dir=out_dir, filter=None)

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
    def insert_filter_to_file_name(file_name, filter=None):
        if filter is None:
            return file_name
        else:
            base, ext = os.path.splitext(file_name)
            return f"{base}_{filter}{ext}"

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

    def combine(self, method='clip'):
        """Combine the frames of the stored files to a master flat.

        Args:
            method (str, optional):
                Method for the frame combination. Can be either 'median' for a
                conventional median combination without propagation of
                uncertainties (since the median does not allow for this) or
                'clip' for sigma clipping and a subsequent estimate of the
                variance of the cube, followed by a mean combination.
        """

        # Type check
        if not isinstance(method, str):
            raise SpecklepyTypeError('MasterFlat.combine()', argname='method', argtype=type(method), expected='str')

        # Read image frames from file
        logger.info("Combining the following files to a master flat:")
        self.means = []

        # flats = None
        for index, file in enumerate(self.files):
            logger.info(f"{index:4}: {file!r}")
            data = fits.getdata(os.path.join(self.file_path, file))
            # data -= np.min(data)  # Avoid negative values that screw up the normalization
            logger.debug(f"{file} (dtype: {data.dtype}, shape: {data.shape})")

            if data.ndim == 2:
                self.means.append(data)
                self.combine_var(np.zeros(data.shape))
                self.combine_mask(np.zeros(data.shape, dtype=bool))

            elif data.ndim == 3:
                logger.info("Computing sigma-slipped statistics of data cube...")
                mean, _, std = sigma_clipped_stats(data=data, axis=0)
                self.means.append(mean)
                self.combine_var(np.square(std))
                self.combine_mask(np.zeros(data.shape, dtype=bool))

        #     # Create a master flat
        #     if flats is None:
        #         flats = data
        #     else:
        #         flats = np.append(flats, data, axis=0)
        # logger.debug(f"'flats' cube has shape: {flats.shape}")

        # Collapse master flat along axis 0
        logger.info(f"Combining flats with {method!r} method...")
        if method == 'median':
            self.image = np.median(np.array(self.means), axis=0)
            # master_flat_var = None
        elif method == 'clip':
            flats = sigma_clip(np.array(self.means), axis=0, masked=True)
            self.image = np.mean(flats, axis=0)
            # master_flat_var = np.var(flats, axis=0)
        else:
            raise SpecklepyValueError('MasterFlat.combine()', argname='method', argvalue=method,
                                      expected="'clip' or 'median'")

        # Normalize the master flat
        logger.info(f"Normalizing master flat in {method!r} mode...")
        if method == 'median':
            norm = np.median(self.image)
            self.image = np.divide(self.image, norm)
            self.var = None
            self.normalized = True
        elif method == 'clip':
            clipped = sigma_clip(self.image)
            norm = np.mean(clipped)
            norm_var = np.var(clipped)
            master_flat_normed = np.divide(self.image, norm)
            self.var = np.divide(self.var, np.square(norm)) + \
                                     np.divide(np.square(self.image), np.power(norm, 4)) * norm_var
            self.image = master_flat_normed
            self.normalized = True

        # Replace masked values by fill_value
        self.mask = combine_masks(sigma_clip(self.image, masked=True, sigma=self.mask_threshold), self.image, self.var)
        logger.debug(f"Replacing masked values by {self.fill_value}...")
        if self.image is not None:
            self.image.mask = self.mask
            self.image, mask_im = self.fill_masked(self.image)
            self.mask = np.logical_or(self.mask, mask_im)
        if self.var is not None:
            self.var.mask = self.mask
            self.var, mask_var = self.fill_masked(self.var)
            self.mask = np.logical_or(self.mask, mask_var)

        # # Store variance in extension
        # self.master_file.data = master_flat_normed
        # if master_flat_normed_var is not None:
        #     self.master_file.new_extension('VAR', data=master_flat_normed_var)
        # self.master_file.new_extension('MASK', data=mask.astype(np.int16))

    def fill_masked(self, masked_array):
        try:
            return masked_array.filled(fill_value=self.fill_value), masked_array.mask
        except AttributeError:
            return masked_array, np.zeros(masked_array.shape, dtype=bool)

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
        logger.info(f"Writing master flat frame to file {self.path!r}")
        hdu_list.writeto(self.path, overwrite=overwrite)

    def run_correction(self, file_list, file_path=None, sub_windows=None):
        """Executes the flat field correction on files in a list.

        Args:
            file_list (list):
                List of files to receive the flat field correction.
            file_path (str, optional):
                Path to the files in file_list.
            sub_windows (list, optional):
                List of sub-window strings, where the entries indicate, how the individual exposure is positioned
                within the master flat field. Should be provided as str-types.
        """

        # Input parameters
        if isinstance(file_list, (list, np.ndarray)):
            pass
        elif isinstance(file_list, str):
            file_list = [file_list]
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_list', type(file_list), 'list')

        # Load master flat field
        gpm = (self.mask == 0)

        # Initialize dummy list if sub_windows is not provided
        if sub_windows is None:
            sub_windows = [None] * len(file_list)
        elif isinstance(sub_windows, str):
            sub_windows = [sub_windows]

        # Iterate through product files
        for file, sub_window_str in zip(file_list, sub_windows):
            logger.info(f"Applying flat field correction on file {file!r}")

            # Add path to file name
            if file_path:
                file = os.path.join(file_path, file)

            # Construct sub-window
            sub_window = SubWindow.from_str(sub_window_str, full=self.sub_window)

            # Open the product files and update
            with fits.open(file, mode='update') as hdu_list:
                extension = 0
                cube = hdu_list[extension].data.astype(float)

                # Check sizes of frames and sub-window
                if sub_window.shape is not None and \
                        (cube.shape[-2] > sub_window.shape[-2] or cube.shape[-1] > sub_window.shape[-1]):
                    logger.warning(f"Unable to apply flat field correction to file {file!r}. Reason may be that "
                                   f"the sub-window covered by the master flat is smaller than the image.")
                else:

                    # Expand 2D image to a one frame cube
                    if cube.ndim == 2:
                        hdu_list[extension].data = np.divide(cube, sub_window(self.image), where=sub_window(gpm))
                    else:
                        # Normalize image/ cube frames by master flat
                        for f in trange(cube.shape[0], desc='Updating frame'):
                            frame = cube[f]
                            try:
                                hdu_list[extension].data[f] = np.divide(frame, sub_window(self.image),
                                                                        where=sub_window(gpm))
                            except ValueError as e:
                                logger.error(e)
                                embed()
                            hdu_list.flush()

                    # Propagate variance
                    if 'VAR' in hdu_list:
                        cube_var = hdu_list['VAR'].data

                        if self.var is None:
                            # The cube stores its variance map already
                            pass
                        else:
                            # Combine variance maps of the cube and the master flat
                            # TODO double check the variance propagation formula
                            # out_var = np.multiply(np.square(np.divide(1, np.square(master_flat))), master_flat_var)
                            out_var = np.add(sub_window(self.var), cube_var)
                            hdu_list['VAR'].data = out_var
                    else:
                        # TODO: Implement a measure for the cube variance
                        if self.var is None:
                            # No variance map is available
                            pass
                        else:
                            # Use the variance map of the master flat
                            hdu_list.append(fits.ImageHDU(name='VAR', data=sub_window(self.var)))

                    # Store the mask
                    if 'MASK' in hdu_list:
                        hdu_list['MASK'].data = np.logical_or(hdu_list['MASK'].data.astype(bool),
                                                              sub_window(gpm)).astype(np.int16)
                    else:
                        hdu_list.append(fits.ImageHDU(data=sub_window(self.mask), name='MASK'))

                    # Update FITS header and store updates
                    hdu_list[0].header.set('HIERARCH SPECKLEPY REDUCTION FLATCORR', default_time_stamp())
                    hdu_list.flush()
