import numpy as np
from astropy.io import fits
from datetime import datetime

from holopy.logging import logging


class Outfile(object):

    def __init__(self, file_list, filename=None):
        if filename is None:
            self.filename = self._make_time_stamp() + ".fits"
        else:
            self.filename = filename

        # Initialize HDU
        hdu = fits.PrimaryHDU()

        # Initialize header
        hdr_input = fits.getheader(file_list[0])
        hdu.header.set('NAXIS', 2)
        hdu.header.set('NAXIS1', hdr_input['NAXIS1'])
        hdu.header.set('NAXIS2', hdr_input['NAXIS2'])

        # Initialize data
        hdu.data = np.zeros( (hdr_input['NAXIS1'], hdr_input['NAXIS2']) )

        # Write to files
        hdulist = fits.HDUList([hdu])
        logging.info("Saving result to: {}".format(self.filename))
        hdulist.writeto(self.filename, overwrite=True)
        logging.info("Saving successful!")


    def _make_time_stamp(self):
        """
        The helper function _make_time_stamp() returns a string:
        'YYYYMMDD_HHMMSS'.
        """
        tmp = str(datetime.now())
        tmp = tmp.split('.')[0]
        tmp = tmp.replace(' ', '_')
        tmp = tmp.replace('-', '')
        tmp = tmp.replace(':', '')
        return tmp
