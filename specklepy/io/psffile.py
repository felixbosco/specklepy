import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logger
from specklepy.io.outfile import Outfile


class PSFFile(Outfile):

    """Outfile with some default parameters and init behaviour for PSF files."""

    def __init__(self, in_file, out_dir, frame_shape, in_dir=None, cards=None, header_card_prefix=None):
        """Create a PSFFile instance.

        Args:
            in_file (str):
                Name of the parent file.
            out_dir (str):
                Name of the directory that the file will be stored in.
            frame_shape (tuple):
                Shape of the PSF frames, which is the box size.
            in_dir (str, optional):
                Path to the input file.
            cards (dict, optional):
                Dictionary of header cards.
            header_card_prefix (str, optional):
        """

        # Create PSF directory, if not existing yet
        if not os.path.exists(out_dir):
            logger.info(f"Creating PSF directory {out_dir}")
            os.makedirs(out_dir)

        # Adapt filename to form the name of the out_file
        # _, out_file = os.path.split(in_file)
        # out_file = out_file.replace('.fits', '_psfs.fits')
        out_file = 'psf_' + os.path.basename(in_file)
        # self.filename = outDir + out_file

        # Type assertion
        if not isinstance(frame_shape, tuple):
            raise SpecklepyTypeError('PSFFile', 'frame_shape', type(frame_shape), 'tuple')

        if cards is None:
            cards = {}
        elif not isinstance(cards, dict):
            raise SpecklepyTypeError('PSFFile', 'cards', type(cards), 'dict')

        if header_card_prefix is None:
            header_card_prefix = ""
        elif not isinstance(header_card_prefix, str):
            raise SpecklepyTypeError('PSFFile', 'header_card_prefix', type(header_card_prefix), 'str')

        # Add name of parent file to header
        cards["FILE NAME"] = os.path.basename(in_file)

        # Derive data shape
        if in_dir is not None:
            hdr_input = fits.getheader(os.path.join(in_dir, in_file))
        else:
            hdr_input = fits.getheader(in_file)
        shape = (hdr_input['NAXIS3'], frame_shape[0], frame_shape[1])

        super().__init__(filename=out_file, path=out_dir, shape=shape, cards=cards,
                         header_card_prefix=header_card_prefix)
