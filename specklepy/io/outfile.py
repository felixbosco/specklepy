import numpy as np
from os import path
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError


class Outfile(object):

    def __init__(self, filename, path=None, data=None, shape=None, extensions=None, header=None, cards=None, timestamp=False, header_card_prefix=None, verbose=True):
        """Instantiate a generic outfile.

        Args:
            filename (str):
                Name of the file that is always referred to in the background.
            path (str, optional):
                Target path under which the new file will be stored.
            data (np.ndarray, optional):
                If provided, the file will be initialized with these data and
                shape argument is ignored. Default is None.
            shape (tuple, dtype=int, optional):
                If provided, the file will be initialized with an empty array
                of this shape. This argument will be ignored if data argument
                is provided. Default is None.
            extensions (list of str, optional):
                List of name(s) of extension(s) that will be initialized right
                away. Default is None.
            header (fits.header, optional):
                Header that will be used to initialize the new Primary HDU.
                Default is None.
            cards (dict, optional):
                Dictionary of cards that will be added to the fits header.
                Default is None.
            timestamp (bool, optional):
                Set to True to automatically add a time stamp to the file name.
                Default is False.
            header_card_prefix (str, optional):
                Prefix of header cards. Default is None.
            verbose (bool, optional):
                Set to False to suppress logging output. Default is True.
        """

        # Input parameters
        if filename is None:
            raise RuntimeError("Outfile did not receive a filename!")
        else:
            self.filename = filename

        if path is None or isinstance(path, str):
            self.path = path
        else:
            raise SpecklepyTypeError('Outfile', argname='path', argtype=type(path), expected='str')

        if data is None or isinstance(data, np.ndarray):
            pass
        else:
            raise SpecklepyTypeError('Outfile', 'data', type(data), 'np.ndarray')

        if shape is None or isinstance(shape, tuple):
            self.shape = shape
        elif isinstance(shape, list):
            self.shape = shape[0]
            self.extshapes = shape[1:]
            if len(self.extshapes) != len(extensions):
                raise ValueError('Outfile received shape as list type, but the number of shapes does not match to the number of extensions!')
        else:
            raise SpecklepyTypeError('Outfile', 'shape', type(shape), 'tuple')

        if shape is not None and data is not None:
            logging.info("Outfile instance ignores shape input as data is provided for initialization!")

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

        if header is None or isinstance(header, fits.header):
            pass
        else:
            raise SpecklepyTypeError('Outfile', 'header', type(header), 'fits.header')

        if not isinstance(timestamp, bool):
            raise SpecklepyTypeError('Outfile', 'timestamp', type(timestamp), 'bool')
        else:
            if timestamp:
                self.filename = filename.replace('.fits', '_{}.fits'.format(self.time_stamp()))

        if header_card_prefix is None:
            self.header_card_prefix = ""
        elif isinstance(header_card_prefix, str):
            # Assert that there are gaps between prefix and card keywords
            if header_card_prefix[-1] != ' ':
                header_card_prefix = header_card_prefix + ' '
            self.header_card_prefix = header_card_prefix
        else:
            raise SpecklepyTypeError('Outfile', 'header_card_prefix', type(header_card_prefix), 'str')

        if isinstance(verbose, bool):
            self.verbose = verbose
        else:
            raise SpecklepyTypeError('Outfile', 'verbose', type(verbose), 'bool')


        # Initialize primary HDU
        hdu = fits.PrimaryHDU(header=header)
        for key in self.cards:
            hdu.header.set(self.header_card_prefix + key, self.cards[key])
        if self.shape is not None and data is None:
            hdu.data = np.zeros(self.shape)
        if data is not None:
            hdu.data = data
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
        if self.verbose:
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
        if self.verbose:
            logging.info("Updating data in {}".format(self.filename))


    def __getitem__(self, extension):
        return fits.getdata(self.filename, extension)#[index]


    def __setitem__(self, extension, data):
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[extension].data = data
            hdulist[extension].header.set('UPDATED', str(datetime.now()))
            hdulist.flush()


    def update_frame(self, frame_index, data, extension=None):
        if extension is None:
            extension = 0
        with fits.open(self.filename, mode='update') as hdulist:
            hdulist[extension].data[frame_index] = data
            hdulist[extension].header.set('UPDATED', str(datetime.now()))
            hdulist.flush()
        # self.data[frame_index] = data


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
            hdu = fits.ImageHDU(data=data, name=name, header=header)
            if index is None:
                hdulist.append(hdu=hdu)
            else:
                hdulist.insert(index=index, hdu=hdu)
            hdulist.flush()
