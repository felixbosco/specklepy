import numpy as np
import os
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logging


class PSFfile(object):

    def __init__(self, inFile, outDir, frame_shape, cards={}, header_prefix="HIERARCH specklepy "):

        # Create PSF directory, if not existing yet
        # outDir += 'psf/'
        if not os.path.exists(outDir):
            logging.info('Creating PSF directory {}'.format(outDir))
            os.makedirs(outDir)

        # Adapt filename to form the name of the outfile
        _, outfile = os.path.split(inFile)
        outfile = outfile.replace('.fits', '_psfs.fits')
        self.filename = outDir + outfile

        # Type assertion
        if not isinstance(frame_shape, tuple):
            raise TypeError("frame_shape argument must have type tuple but was given as {}.".format(type(frame_shape)))

        # Initialize HDU
        hdu = fits.PrimaryHDU()

        # Initialize header
        hdr_input = fits.getheader(inFile)

        # Add cards from cards dictionary to header
        for key in cards:
            hdu.header.set(header_prefix + key, cards[key])

        # Add name of parent file to header
        hdu.header.set(header_prefix + "FILE NAME", os.path.basename(inFile))

        # Add parameter of reference stars to header

        # Initialize data
        hdu.data = np.zeros( (hdr_input['NAXIS3'], frame_shape[0], frame_shape[1]) )
        hdu.header.set('DATE', str(datetime.now()))

        # Write to files
        logging.info("Initializing PSF file {}".format(self.filename))
        hdu.writeto(self.filename, overwrite=True)


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


    def __getitem__(self, index):
        return fits.getdata(self.filename)[index]


    def __setitem__(self, index, data):
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[0].data[index] = data
            hdulist[0].header.set('UPDATED', str(datetime.now()))
            hdulist.flush()


    def update_frame(self, frame_index, data):
        self[frame_index] = data
