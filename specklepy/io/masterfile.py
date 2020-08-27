from os import path
from astropy.io import fits

from specklepy.io.outfile import Outfile


class MasterFile(Outfile):

    def __init__(self, filename, files, shape=None, cards=None, out_dir=None, header_card_prefix="HIERARCH SPECKLEPY"):

        if cards is None:
            cards = {}

        # Add list of files to header
        for index, file in enumerate(files):
            cards["FILE {}".format(index)] = path.basename(file)

        if shape is None:
            hdr_input = fits.getheader(files[0])
            try:
                shape = (hdr_input['NAXIS2'], hdr_input['NAXIS3'])
            except KeyError:
                shape = (hdr_input['NAXIS1'], hdr_input['NAXIS2'])

        super().__init__(filename=filename, shape=shape, extensions=None, cards=cards, timestamp=False, path=out_dir,
                         header_card_prefix=header_card_prefix)
