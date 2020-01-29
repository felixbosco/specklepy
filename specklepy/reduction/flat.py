import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

from specklepy.exceptions import SpecklepyTypeError
from specklepy.io.masterfile import MasterFile
from specklepy.logging import logging



class MasterFlat(object):

    def __init__(self, files, filename='MasterFlat.fits', file_path=None):
        # Store input parameters
        if isinstance(files, list):
            self.files = files
        elif isinstance(files, Table):
            is_flat_file = files['OBSTYPE'] == 'FLAT'
            self.files = files['FILE'][is_flat_file]
        else:
            raise SpecklepyTypeError('MasterFlat', 'files', type(files), 'astropy.table.Table')

        if isinstance(filename, str):
            self.filename = filename
        else:
            raise SpecklepyTypeError('MasterFlat', 'filename', type(filename), 'str')

        if isinstance(file_path, str):
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('MasterFlat', 'file_path', type(file_path), 'str')


    def make_master_flat(self):
        logging.info("Combining the following files to a master flat:")

        # Read image frames from file
        for index, file in enumerate(self.files):
            logging.info("{:4}: {}".format(index, file))
            data = fits.getdata(os.path.join(self.file_path, file))

            # Create a master flat
            if index is 0:
                master_flat = data
            else:
                np.append(master_flat, data, axis=0)

        # Collapse master flat
        master_flat = np.median(master_flat, axis=0)

        # Normalize the master flat
        logging.info("Normalizing master flat {}".format(self.filename))
        master_flat /= np.median(master_flat)

        # Store master flat to file
        if not hasattr(self, 'outfile'):
            self.outfile = MasterFile(self.filename, files=self.files, shape=master_flat.shape, header_prefix='HIERARCH SPECKLEPY')
        self.outfile.data = master_flat

