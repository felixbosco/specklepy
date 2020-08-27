import os

from astropy.io import fits

from specklepy.io.outfile import Outfile


class MasterFile(Outfile):

    def __init__(self, filename, files, shape=None, cards=None, in_dir=None, out_dir=None,
                 header_card_prefix="HIERARCH SPECKLEPY"):

        # Apply fall back value
        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards["FILE {}".format(index)] = os.path.basename(file)

        # Derive shape from FITS header
        if shape is None:
            # Read header information
            if in_dir:
                hdr_input = fits.getheader(os.path.join(in_dir, files[0]))
            else:
                hdr_input = fits.getheader(files[0])

            # Derive shape from header entries
            if 'NAXIS3' in hdr_input:
                shape = (hdr_input['NAXIS2'], hdr_input['NAXIS3'])
            else:
                shape = (hdr_input['NAXIS1'], hdr_input['NAXIS2'])

        # Initialize as parent class instance
        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False, path=out_dir,
                         header_card_prefix=header_card_prefix)
