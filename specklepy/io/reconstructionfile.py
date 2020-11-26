import os
from astropy.io import fits

from specklepy.io.outfile import Outfile


class ReconstructionFile(Outfile):

    def __init__(self, filename, files, shape=None, cards=None, in_dir=None, header_card_prefix="HIERARCH SPECKLEPY"):

        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards[f"FILE{index:04} NAME"] = os.path.basename(file)
            if in_dir:
                file = os.path.join(in_dir, file)
            cards[f"FILE{index:04} FRAMES"] = fits.getheader(file)['NAXIS3']

        # Derive shape from FITS header
        if shape is None:
            # Read header information
            if in_dir:
                hdr_input = fits.getheader(os.path.join(in_dir, files[0]))
            else:
                hdr_input = fits.getheader(files[0])

            # Derive shape from header entries
            shape = self.extract_frame_shape(hdr_input)

        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False,
                         header_card_prefix=header_card_prefix)
