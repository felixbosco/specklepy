from datetime import datetime
from IPython import embed
import numpy as np
import os
from tqdm import trange

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.masterfile import MasterFile
from specklepy.logging import logger
from specklepy.reduction.subwindow import SubWindow


class MasterFlat(object):

    def __init__(self, file_list, file_name='MasterFlat.fits', file_path=None, out_dir=None, new=True):

        # Store input parameters
        if isinstance(file_list, (list, np.ndarray)):
            self.files = file_list
        elif isinstance(file_list, Table):
            is_flat_file = file_list['OBSTYPE'] == 'FLAT'
            self.files = file_list['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_list', type(file_list), 'astropy.table.Table')

        if isinstance(file_name, str):
            self.file_name = file_name
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_name', type(file_name), 'str')

        if isinstance(file_path, str) or file_path is None:
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_path', type(file_path), 'str')

        # Create an output file
        self.master_file = MasterFile(self.file_name, files=self.files, in_dir=file_path, out_dir=out_dir,
                                      initialize=new)

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
        logger.info("Combining the following file list to a master flat:")
        flats = None
        for index, file in enumerate(self.files):
            logger.info(f"{index:4}: {file}")
            data = fits.getdata(os.path.join(self.file_path, file))
            data -= np.min(data)  # Avoid negative values that screw up the normalization
            logger.debug(f"{file} (dtype: {data.dtype}, shape: {data.shape})")

            # Create a master flat
            if flats is None:
                flats = data
            else:
                flats = np.append(flats, data, axis=0)
        logger.debug(f"'flats' cube has shape: {flats.shape}")

        # Collapse master flat along axis 0
        logger.info(f"Combining flats with {method!r} method...")
        if method == 'median':
            master_flat = np.median(flats, axis=0)
            master_flat_var = None
        elif method == 'clip':
            flats = sigma_clip(flats, axis=0, masked=True)
            master_flat = np.mean(flats, axis=0)
            master_flat_var = np.var(flats, axis=0)
        else:
            raise SpecklepyValueError('MasterFlat.combine()', argname='method', argvalue=method,
                                      expected="'clip' or 'median'")
        del flats

        # Normalize the master flat
        logger.info(f"Normalizing master flat in {method!r} mode...")
        if method == 'median':
            norm = np.median(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = None
        elif method == 'clip':
            clipped = sigma_clip(master_flat)
            norm = np.mean(clipped)
            norm_var = np.var(clipped)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = np.divide(master_flat_var, np.square(norm)) + \
                                     np.divide(np.square(master_flat), np.power(norm, 4)) * norm_var
        else:
            master_flat_normed = None
            master_flat_normed_var = None

        # Replace masked values by NaNs
        logger.debug(f"Replacing masked values by Nan...")
        if master_flat_normed is not None and hasattr(master_flat_normed, 'mask'):
            if np.sum(master_flat_normed.mask) == 0:
                master_flat_normed = master_flat_normed.data
            else:
                master_flat_normed = np.where(master_flat_normed.mask, np.nan, master_flat_normed)
                master_flat_normed = master_flat_normed.data
        if master_flat_normed_var is not None and hasattr(master_flat_normed_var, 'mask'):
            if np.sum(master_flat_normed_var.mask) == 0:
                master_flat_normed_var = master_flat_normed_var.data
            else:
                master_flat_normed_var = np.where(master_flat_normed_var.mask, np.nan, master_flat_normed_var)
                master_flat_normed_var = master_flat_normed_var.data

        # Store variance in extension
        self.master_file.data = master_flat_normed
        if master_flat_normed_var is not None:
            self.master_file.new_extension('VAR', data=master_flat_normed_var)

    def run_correction(self, file_list, file_path=None, sub_windows=None, full_window=None):
        """Executes the flat field correction on files in a list.

        Args:
            file_list (list):
                List of files to receive the flat field correction.
            file_path (str, optional):
                Path to the files in file_list.
            sub_windows (list, optional):
                List of sub-window strings, where the entries indicate, how the individual exposure is positioned
                within the master flat field. Should be provided as str-types.
            full_window (str, optional):
                String describing the full window in the same format as the sub-window strings for providing relative
                coordinates.
        """

        # Input parameters
        if isinstance(file_list, (list, np.ndarray)):
            pass
        elif isinstance(file_list, str):
            file_list = [file_list]
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_list', type(file_list), 'list')

        # Load master flat field
        master_flat = self.master_file.data
        if self.master_file.has_extension('VAR'):
            master_flat_var = self.master_file['VAR']
            propagate_uncertainties = True
        else:
            propagate_uncertainties = False

        # Initialize dummy list if sub_windows is not provided
        if sub_windows is None:
            sub_windows = [None] * len(file_list)
        elif isinstance(sub_windows, str):
            sub_windows = [sub_windows]

        # Iterate through product files
        for file, sub_window_str in zip(file_list, sub_windows):
            logger.info(f"Applying flat field correction on file {file}")

            # Add path to file name
            if file_path:
                file = os.path.join(file_path, file)

            # Construct sub-window
            sub_window = SubWindow.from_str(sub_window_str, full=full_window)

            # Open the product files and update
            with fits.open(file, mode='update') as hdu_list:
                extension = 0
                cube = hdu_list[extension].data.astype(float)

                # Expand 2D image to a one frame cube
                if cube.ndim == 2:
                    hdu_list[extension].data = np.divide(cube, sub_window(master_flat))
                else:
                    # Normalize image/ cube frames by master flat
                    for f in trange(cube.shape[0], desc='Updating frame'):
                        frame = cube[f]
                        try:
                            hdu_list[extension].data[f] = np.divide(frame, sub_window(master_flat))
                        except ValueError as e:
                            logger.error(e)
                            embed()
                        hdu_list.flush()

                # Propagate variance if the master flat has this information
                if self.master_file.has_extension('VAR'):
                    master_flat_var = self.master_file['VAR']

                    # Create VAR HDU or propagate the variances
                    if 'VAR' in hdu_list:
                        cube_var = hdu_list['VAR'].data
                        # TODO double check the variance propagation formula
                        # hdu_list['VAR'].data = np.multiply(np.square(np.divide(1, np.square(master_flat))),
                        # master_flat_var)
                        hdu_list['VAR'].data += master_flat_var
                    else:
                        var_hdu = fits.ImageHDU(data=master_flat_var, name='VAR', header=None)
                        hdu_list.append(var_hdu)

                hdu_list[0].header.set('FLATCORR', str(datetime.now()))
                hdu_list.flush()
