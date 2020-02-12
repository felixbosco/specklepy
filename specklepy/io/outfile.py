import numpy as np
from os import path
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError



class Outfile(object):

    def __init__(self, filename, shape=None, extensions=None, cards=None, timestamp=False, hprefix=None):
        """Instantiate a generic outfile.

        Args:
            filename (str):
            shape (tuple, dtype=int, optional):
            cards (dict, optional):
            timestamp (bool, optional):
                Set to True to automatically add a time stamp to the file name.
                Default is False.
            hprefix (str, optional):
                Prefix of header cards. Default is None.
        """

        # Input parameters
        if filename is None:
            raise RuntimeError("Outfile did not receive a filename!")
        else:
            self.filename = filename

        if shape is None or isinstance(shape, tuple):
            self.shape = shape
        elif isinstance(shape, list):
            self.shape = shape[0]
            self.extshapes = shape[1:]
            if len(self.extshapes) != len(extensions):
                raise ValueError('Outfile received shape as list type, but the number of shapes does not match to the number of extensions!')
        else:
            raise SpecklepyTypeError('Outfile', 'shape', type(shape), 'tuple')

        if isinstance(extensions, str):
            self.extensions = [extensions]
        elif extensions is None or isinstance(extensions, list):
            self.extensions = extensions
        else:
            raise SpecklepyTypeError('Outfile', 'extensions', type(extensions), 'list')

        if cards is None:
            self.cards = {}
        elif isinstance(cards, dict):
            self.cards = cards
        else:
            raise SpecklepyTypeError('Outfile', 'cards', type(cards), 'dict')

        if not isinstance(timestamp, bool):
            raise SpecklepyTypeError('Outfile', 'timestamp', type(timestamp), 'bool')
        else:
            if timestamp:
                self.filename = filename.replace('.fits', '_{}.fits'.format(self.time_stamp()))

        if hprefix is None:
            self.hprefix = ""
        elif isinstance(hprefix, str):
            # Assert that there are gaps between prefix and card keywords
            if hprefix[-1] != ' ':
                hprefix = hprefix + ' '
            self.hprefix = hprefix
        else:
            raise SpecklepyTypeError('Outfile', 'hprefix', type(shape), 'str')

        # Initialize primary HDU
        hdu = fits.PrimaryHDU()
        for key in self.cards:
            hdu.header.set(self.hprefix + key, self.cards[key])
        if self.shape is not None:
            hdu.data = np.zeros(self.shape)
        hdu.header.set('DATE', str(datetime.now()))

        # Create a HDU list with the primary and append extensions
        hdulist = fits.HDUList([hdu])

        if self.extensions is not None:
            for index, extension in enumerate(self.extensions):
                exthdu = fits.ImageHDU(name=extension)
                if hasattr(self, 'extshapes'):
                    exthdu.data = np.zeros(self.extshapes[index +1])
                elif self.shape is not None:
                    exthdu.data = np.zeros(self.shape)
                hdulist.append(exthdu)

        # Write HDU list to file
        logging.info("Initializing file {}".format(self.filename))
        hdulist.writeto(self.filename, overwrite=True)


    def time_stamp(self):
        """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
        return datetime.now().strftime('%Y%m%d_%H%M%S')


    @property
    def data(self):
        return fits.getdata(self.filename)


    @data.setter
    def data(self, data):
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[0].data = data
            hdulist[0].header.set('UPDATED', str(datetime.now()))
            hdulist.flush()
        logging.info("Updating data in {}".format(self.filename))


    def __getitem__(self, extension):
        return fits.getdata(self.filename, extension)#[index]


    def __setitem__(self, extension, data):
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[extension].data = data
            hdulist[extension].header.set('UPDATED', str(datetime.now()))
            hdulist.flush()


    def update_frame(self, frame_index, data):
        self.data[frame_index] = data


    def new_extension(self, name, data=None, header=None, index=None):
        """Create a new fits extension.

        All arguments are passed to the new HDU in the HDUlist of instances
        file, i.e. 'self.filename'.

        Args:
            name (str):
                Name of the new extension/ HDU.
            data (array, optional):
                Data array, passed to the new HDU.
            header (fits.header, optional):
                Fits header, passed to the new HDU.
            index (int, optional):
                Index of the extension, before which the new hdu is inserted.
                None translates into appending the new HDU at the end of the
                HDUlist.
        """

        with fits.open(self.filename, mode='update') as hdulist:
            hdulist.info()
            hdu = fits.ImageHDU(data=data, name=name, header=header)
            if index is None:
                hdulist.append(hdu=hdu)
            else:
                hdulist.insert(index=index, hdu=hdu)
            hdulist.info()
            hdulist.flush()
