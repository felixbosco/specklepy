import numpy as np
from os import path
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logging


class Outfile(object):

    def __init__(self, files, filename=None, cards={}, header_prefix="HIERARCH specklepy "):
        if filename is None:
            self.filename = self._make_time_stamp() + ".fits"
        else:
            self.filename = filename

        # Initialize HDU
        hdu = fits.PrimaryHDU()

        # Initialize header
        hdr_input = fits.getheader(files[0])
        # hdu.header.set('NAXIS', 2)
        # hdu.header.set('NAXIS1', hdr_input['NAXIS1'])
        # hdu.header.set('NAXIS2', hdr_input['NAXIS2'])

        # Add cards from cards dictionary to header
        for key in cards:
            hdu.header.set(header_prefix + key, cards[key])

        # Add list of files to header
        for index, file in enumerate(files):
            hdu.header.set(header_prefix + "FILE {} NAME".format(index), path.basename(file))
            hdu.header.set(header_prefix + "FILE {} FRAMENUMBER".format(index), fits.getheader(file)['NAXIS3'])

        # Initialize data
        hdu.data = np.zeros( (hdr_input['NAXIS1'], hdr_input['NAXIS2']) )
        hdu.header.set('DATE', str(datetime.now()))

        # Write to files
        logging.info("Saving results to: {}".format(self.filename))
        hdu.writeto(self.filename, overwrite=True)
        logging.info("File initialized!")


    @property
    def data(self):
        return fits.getdata(self.filename)

    @data.setter
    def data(self, data):
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[0].data = data
            hdulist[0].header.set('DATE', str(datetime.now()))
            hdulist.flush()
        logging.info("Updating data in {}".format(self.filename))


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
