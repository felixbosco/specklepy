import os
from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logger
from specklepy.io.outfile import Outfile


class PSFfile(Outfile):

    """Outfile with some default parameters and init behaviour for PSF files."""

    def __init__(self, inFile, outDir, frame_shape, cards=None, header_card_prefix=None):
        """Create a PSFfile instance.

        Args:
            inFile (str):
                Name of the parent file.
            outDir (str):
                Name of the directory that the file will be stored in.
            frame_shape (tuple):
                Shape of the PSF frames, which is the box size.
            cards (dict):
                Dictionary of header cards.
            header_card_prefix (str, optional):
        """

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
            raise SpecklepyTypeError('PSFfile', 'frame_shape', type(frame_shape), 'tuple')

        if cards is None:
            cards = {}
        elif not isinstance(cards, dict):
            raise SpecklepyTypeError('PSFfile', 'cards', type(cards), 'dict')

        if header_card_prefix is None:
            header_card_prefix = ""
        elif not isinstance(header_card_prefix, str):
            raise SpecklepyTypeError('PSFfile', 'header_card_prefix', type(header_card_prefix), 'str')

        # Add name of parent file to header
        cards["FILE NAME"] = os.path.basename(inFile)

        hdr_input = fits.getheader(inFile)
        shape = (hdr_input['NAXIS3'], frame_shape[0], frame_shape[1])

        super().__init__(filename=os.path.join(outDir, outfile), shape=shape, cards=cards,
                         header_card_prefix=header_card_prefix)
