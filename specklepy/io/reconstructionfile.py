import os
from astropy.io import fits

from specklepy.io.outfile import Outfile


class ReconstructionFile(Outfile):

    def __init__(self, filename, files, shape=None, cards=None, in_dir=None, header_card_prefix="HIERARCH SPECKLEPY"):

        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards[f"SOURCE FILE{index:04} NAME"] = os.path.basename(file)
            if in_dir:
                file = os.path.join(in_dir, file)
            cards[f"SOURCE FILE{index:04} FRAMES"] = self.extract_frame_number(fits.getheader(file))

        # Derive frame shape from FITS header
        if shape is None:
            example_path = os.path.join(in_dir, files[0]) if in_dir is not None else files[0]
            shape = self.extract_frame_shape(fits.getheader(example_path))

        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False,
                         header_card_prefix=header_card_prefix)
