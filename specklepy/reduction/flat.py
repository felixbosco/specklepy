import numpy as np
import os
from tqdm import trange

from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.masterfile import MasterFile
from specklepy.io.reductionfile import ReductionFile
from specklepy.logging import logger


class MasterFlat(object):

    def __init__(self, file_list, file_name='MasterFlat.fits', file_path=None, out_dir=None):

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

        if isinstance(file_path, str):
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_path', type(file_path), 'str')

        # Create an output file
        self.master_file = MasterFile(self.file_name, files=self.files, out_dir=out_dir)

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
        for index, file in enumerate(self.files):
            logger.info("{:4}: {}".format(index, file))
            data = fits.getdata(os.path.join(self.file_path, file))

            # Create a master flat
            if index is 0:
                flats = data
            else:
                np.append(flats, data, axis=0)

        # Collapse master flat along axis 0
        if method is 'median':
            master_flat = np.median(flats, axis=0)
            master_flat_var = None
        elif method is 'clip':
            flats = sigma_clip(flats, axis=0, masked=True)
            master_flat = np.mean(flats, axis=0)
            master_flat_var = np.var(flats, axis=0)
        else:
            raise SpecklepyValueError('MasterFlat.combine()', argname='method', argvalue=method,
                                      expected="'clip' or 'median'")
        del flats

        # Normalize the master flat
        logger.info(f"Normalizing master flat in {method} mode...")
        if method is 'median':
            norm = np.median(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = None
        elif method is 'clip':
            norm = np.mean(master_flat)
            norm_var = np.var(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = np.divide(master_flat_var, np.square(norm)) + \
                                     np.divide(np.square(master_flat), np.power(norm, 4)) * norm_var
        else:
            master_flat_normed = None
            master_flat_normed_var = None

        # # Store master flat to file
        # if not hasattr(self, 'masterfile'):
        #     self.master_file = MasterFile(self.file_name, files=self.files, shape=master_flat_normed.shape,
        #                                   header_card_prefix='HIERARCH SPECKLEPY')

        # Replace masked values by NaNs
        master_flat_normed = np.where(master_flat_normed.mask, np.nan, master_flat_normed)
        self.master_file.data = master_flat_normed

        # Store variance in extension
        if 'master_flat_normed_var' in locals():
            # Replace masked values by NaNs
            master_flat_normed_var = np.where(master_flat_normed_var.mask, np.nan, master_flat_normed_var)
            self.master_file.new_extension('VAR', data=master_flat_normed_var)

    def run_correction(self, file_list, file_path=None):
        """Executes the flat field correction on files in a list.

        Args:
            file_list (list):
                List of files to receive the flat field correction.
            file_path (str, optional):
                Path to the files in file_list.
        """

        # Input parameters
        if isinstance(file_list, (list, np.ndarray)):
            pass
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_list', type(file_list), 'list')

        # Load master flat field
        master_flat = self.master_file.data
        if self.master_file.has_extension('VAR'):
            master_flat_var = self.master_file['VAR']
            propagate_uncertainties = True
        else:
            propagate_uncertainties = False

        # Iterate through product files
        for file in file_list:
            logger.info(f"Applying flat field correction on file {file}")

            # Add path to file name
            if file_path:
                file = os.path.join(file_path, file)

            # Open the product files and update
            with fits.open(file, mode='update') as hdu_list:
                extension = 0
                cube = hdu_list[extension].data

                # Expand 2D image to a one frame cube
                if cube.ndim == 2:
                    cube = np.expand_dims(cube, 0)

                # Normalize image/ cube frames by master flat
                for f in trange(cube.shape[0], desc='Updating frame'):
                    frame = hdu_list[extension].data[f]
                    hdu_list[extension].data[f] = np.divide(frame, master_flat)
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
