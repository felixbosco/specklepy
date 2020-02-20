import numpy as np
import os
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logger
from specklepy.io.outfile import Outfile



class PSFfile(Outfile):

    def __init__(self, inFile, outDir, frame_shape, cards=None, header_card_prefix=None):

        # Create PSF directory, if not existing yet
        if not os.path.exists(outDir):
            logger.info('Creating PSF directory {}'.format(outDir))
            os.makedirs(outDir)

        # Adapt filename to form the name of the outfile
        _, outfile = os.path.split(inFile)
        outfile = outfile.replace('.fits', '_psfs.fits')
        # self.filename = outDir + outfile

        # Type assertion
        if not isinstance(frame_shape, tuple):
            raise TypeError("frame_shape argument must have type tuple but was given as {}.".format(type(frame_shape)))

        if cards is None:
            cards = {}
        elif not isinstance(cards, dict):
            raise TypeError("PSFfile received cards argument of type {}, but needs to be dict type!".format(type(cards)))

        if header_card_prefix is None:
            header_card_prefix = ""
        elif not isinstance(header_card_prefix, str):
            raise TypeError("PSFfile received header_card_prefix argument of type {}, but needs to be str type!".format(type(header_card_prefix)))

        # Add name of parent file to header
        cards["FILE NAME"] = os.path.basename(inFile)

        hdr_input = fits.getheader(inFile)
        shape = (hdr_input['NAXIS3'], frame_shape[0], frame_shape[1])

        super().__init__(filename=outDir + outfile, shape=shape, cards=cards, header_card_prefix=header_card_prefix)
