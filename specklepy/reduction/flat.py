import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.masterfile import MasterFile
from specklepy.io.reductionfile import ReductionFile
from specklepy.logging import logger


class MasterFlat(object):

    def __init__(self, file_list, filename='MasterFlat.fits', file_path=None):
        # Store input parameters
        if isinstance(file_list, (list, np.ndarray)):
            self.files = file_list
        elif isinstance(file_list, Table):
            is_flat_file = file_list['OBSTYPE'] == 'FLAT'
            self.files = file_list['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'filelist', type(file_list), 'astropy.table.Table')

        if isinstance(filename, str):
            self.filename = filename
        else:
            raise SpecklepyTypeError('MasterFlat', 'filename', type(filename), 'str')

        if isinstance(file_path, str):
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_path', type(file_path), 'str')

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
        elif method is 'clip':
            flats = sigma_clip(flats, axis=0, masked=True)
            master_flat = np.mean(flats, axis=0)
            master_flat_var = np.var(flats, axis=0)
        else:
            raise SpecklepyValueError('MasterFlat.combine()', argname='method', argtype=method,
                                      expected="'clip' or 'median'")
        del flats

        # Normalize the master flat
        logger.info(f"Normalizing master flat in {method} mode...")
        if method is 'median':
            norm = np.median(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
        elif method is 'clip':
            norm = np.mean(master_flat)
            norm_var = np.var(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = np.divide(master_flat_var, np.square(norm)) + \
                                     np.divide(np.square(master_flat), np.power(norm, 4)) * norm_var

        # Store master flat to file
        if not hasattr(self, 'masterfile'):
            self.masterfile = MasterFile(self.filename, files=self.files, shape=master_flat_normed.shape,
                                         header_card_prefix='HIERARCH SPECKLEPY')
        # Replace masked values by NaNs
        master_flat_normed = np.where(master_flat_normed.mask, np.nan, master_flat_normed)
        self.masterfile.data = master_flat_normed

        # Store variance in extension
        if 'master_flat_normed_var' in locals():
            # Replace masked values by NaNs
            master_flat_normed_var = np.where(master_flat_normed_var.mask, np.nan, master_flat_normed_var)
            self.masterfile.new_extension('VAR', data=master_flat_normed_var)

    def run_correction(self, file_list, filter=None, prefix=None, save_dir=None):
        """Executes the flat field correction on files in a list.

        Args:
            file_list (list or astropy.Table):
                List of files to receive the flat field correction.
            filter (dict, optional):
                Dictionary passed to file filtering function within
                FileManager.
            prefix (str, optional):
                File prefix for the output files.
            save_dir (str, optional):
                Directory, in which the files will be stored.
        """

        # Input parameters
        if isinstance(file_list, (list, np.ndarray)):
            pass
        # elif isinstance(file_list, Table):
        #     is_flat_file = file_list['OBSTYPE'] == 'FLAT'
        #     files = file_list['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'filelist', type(file_list), 'astropy.table.Table')

        master_flat = self.masterfile.data

        flatfield_corrected_files = {}
        for file in file_list:
            logger.info(f"Applying flat field correction on file {file}")

            # Read data and header information
            image, header = fits.getdata(os.path.join(self.file_path, file), header=True)

            # Apply flat field correction
            flatfielded_data = np.divide(image, master_flat)

            corrected_file = ReductionFile(file=os.path.join(self.file_path, file),
                                           data=flatfielded_data,
                                           prefix=prefix,
                                           path=self.file_path,
                                           reduction='FLATCORR')

            # Propagate variance if the master flat has this information
            if self.masterfile.has_extension('VAR'):
                master_flat_var = self.masterfile['VAR']
                # TODO double check the variance propagation formula
                image_var = np.multiply(np.square(np.divide(image, np.square(master_flat))), master_flat_var)
                corrected_file.new_extension(name='VAR', data=image_var)
            flatfield_corrected_files[file] = corrected_file.filename

        return flatfield_corrected_files
