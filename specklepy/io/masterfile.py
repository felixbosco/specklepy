import os

from astropy.io import fits

from specklepy.io.outfile import Outfile


class MasterFile(Outfile):

    def __init__(self, filename, files, shape=None, cards=None, in_dir=None, out_dir=None,
                 header_card_prefix="HIERARCH SPECKLEPY", initialize=True):

        # Apply fall back value
        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards[f"SOURCE FILE{index:04} NAME"] = os.path.basename(file)

        # Derive frame shape from FITS header
        if shape is None:
            example_path = os.path.join(in_dir, files[0]) if in_dir is not None else files[0]
            shape = self.extract_frame_shape(fits.getheader(example_path))

        # Initialize as parent class instance
        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False, path=out_dir,
                         header_card_prefix=header_card_prefix, initialize=initialize)
