import os
import numpy as np
from astropy.io import fits

from specklepy.io.outfile import Outfile



def make_master_flat(params):

    is_flat_file = params.fileList['OBSTYPE'] == 'FLAT'
    flat_files = params.fileList['FILE'][is_flat_file]

    for index, file in enumerate(flat_files):
        # Read data from file and median combine if cube
        data = fits.getdata(os.path.join(params.filePath, file))
        if data.ndim == 3:
            data = np.median(data, axis=0)

        # Create a master flat
        if index is 0:
            master_flat = Outfile(params.masterFlatFile, shape=(1024, 1024), hprefix='HIERARCH SPECKLEPY')

        master_flat.data += data

    # Normalize the master flat
    logging.info("Normalizing master flat {}".format(params.masterFlatFile))
    master_flat.data /= np.median(master_flat.data)

