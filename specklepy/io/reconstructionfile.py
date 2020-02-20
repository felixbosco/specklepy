import numpy as np
from os import path
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logger
from specklepy.io.outfile import Outfile


class ReconstructionFile(Outfile):

    def __init__(self, filename, files, cards=None, header_card_prefix="HIERARCH SPECKLEPY"):

        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards["FILE {}".format(index)] = path.basename(file)
            cards["FILE {} FRAMES".format(index)] = fits.getheader(file)['NAXIS3']

        hdr_input = fits.getheader(files[0])
        shape = (hdr_input['NAXIS1'], hdr_input['NAXIS2'])

        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False, header_card_prefix=header_card_prefix)
