import os
import numpy as np
from datetime import datetime
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clip

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.masterfile import MasterFile
from specklepy.logging import logging
from specklepy.utils.plot import imshow



class MasterFlat(object):

    def __init__(self, filelist, filename='MasterFlat.fits', file_path=None):
        # Store input parameters
        if isinstance(filelist, (list, np.ndarray)):
            self.files = filelist
        elif isinstance(filelist, Table):
            is_flat_file = filelist['OBSTYPE'] == 'FLAT'
            self.files = filelist['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'filelist', type(filelist), 'astropy.table.Table')

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

        if isinstance(method, str):
            if method not in ['clip', 'median']:
                raise SpecklepyValueError('MasterFlat.combine()', argname='method', argtype=method, expected="'clip' or 'median'")
        else:
            raise SpecklepyTypeError('MasterFlat.combine()', argname='method', argtype=type(method), expected='str')

        logging.info("Combining the following filelist to a master flat:")

        # Read image frames from file
        for index, file in enumerate(self.files):
            logging.info("{:4}: {}".format(index, file))
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
        del flats

        # Normalize the master flat
        logging.info(f"Normalizing master flat in {method} mode...")
        if method is 'median':
            norm = np.median(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
        elif method is 'clip':
            norm = np.mean(master_flat)
            norm_var = np.var(master_flat)
            master_flat_normed = np.divide(master_flat, norm)
            master_flat_normed_var = np.divide(master_flat_var, np.square(norm)) + np.divide(np.square(master_flat), np.power(norm, 4)) * norm_var


        # Store master flat to file
        if not hasattr(self, 'masterfile'):
            self.masterfile = MasterFile(self.filename, files=self.files, shape=master_flat_normed.shape, header_prefix='HIERARCH SPECKLEPY')
        master_flat_normed = np.where(master_flat_normed.mask, np.nan, master_flat_normed) # Replace masked values by NaNs
        self.masterfile.data = master_flat_normed

        # Store variance in extension
        if 'master_flat_normed_var' in locals():
            master_flat_normed_var = np.where(master_flat_normed_var.mask, np.nan, master_flat_normed_var) # Replace masked values by NaNs
            self.masterfile.new_extension('VAR', data=master_flat_normed_var)



    def run_correction(self, filelist, filter=None, prefix=None, savedir=None):
        """Executes the flat field correction of a filelist.

        Args:
            filelist (list or astropy.Table):
                List of files to receive the flatfield correction.
            filter (dict, optional):
                Dictionary passed to file filtering function within
                FileManager.
            prefix (str, optional):
                File prefix for the output files.
            savedir (str, optional):
                Directory, in which the files will be stored.
        """

        # Input parameters
        if isinstance(filelist, (list, np.ndarray)):
            pass
        # elif isinstance(filelist, Table):
        #     is_flat_file = filelist['OBSTYPE'] == 'FLAT'
        #     files = filelist['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'filelist', type(filelist), 'astropy.table.Table')

        master_flat = self.masterfile.data
        try:
            master_flat_var = self.masterfile['VAR']
        except KeyError:
            # masterfile carries no variance information
            pass

        flatfield_corrected_files = {}
        for file in filelist:
            logging.info(f"Applying flat field correction on file {file}")
            # Create output file name
            corrected_file = prefix + file
            # if savedir is not None:
            #     corrected_file = os.path.join(savedir, corrected_file)
            # else:
            #     corrected_file = os.path.join(self.file_path, corrected_file)
            flatfield_corrected_files[file] = corrected_file

            # Read data and header information
            image, header = fits.getdata(os.path.join(self.file_path, file), header=True)

            # Apply flat field correction
            if 'master_flat_var' in locals():
                # TODO double check the variance propagation formula
                image_var = np.multiply(np.square(np.divide(image, np.square(master_flat))), master_flat_var)
            image = np.divide(image, master_flat)

            # Save corrected file to update in the FileManager
            header.set('PIPELINE', 'SPECKLEPY')
            header.set('FLATCORR', str(datetime.now()))
            primary_hdu = fits.PrimaryHDU(data=image, header=header)
            hdulist = fits.HDUList(hdus=[primary_hdu])
            if 'image_var' in locals():
                var_hdu = fits.ImageHDU(data=image_var, name='VAR')
                hdulist.append(var_hdu)
            if savedir is not None:
                corrected_file = os.path.join(savedir, corrected_file)
            else:
                corrected_file = os.path.join(self.file_path, corrected_file)
            hdulist.writeto(corrected_file)

        return flatfield_corrected_files

