import numpy as np
from os import path
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logging
from specklepy.io.outfile import Outfile


class RECfile(Outfile):

    def __init__(self, files, filename=None, cards=None, header_prefix="HIERARCH SPECKLEPY"):

        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards["FILE {} NAME".format(index)] = path.basename(file)
            cards["FILE {} FRAMENUMBER".format(index)] = fits.getheader(file)['NAXIS3']

        hdr_input = fits.getheader(files[0])
        shape = (hdr_input['NAXIS1'], hdr_input['NAXIS2'])

        super().__init__(filename=filename, shape=None, extensions=None, cards=cards, timestamp=False, hprefix=header_prefix)
