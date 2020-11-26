import numpy as np
import os
from astropy.io import fits
from datetime import datetime

from specklepy.logging import logger
from specklepy.exceptions import SpecklepyTypeError
from specklepy.utils.time import default_time_stamp


class Outfile(object):

    def __init__(self, filename, path=None, data=None, shape=None, extensions=None, header=None, cards=None,
                 initialize=True, timestamp=False, header_card_prefix=None, verbose=True):
        """Instantiate a generic outfile.

        Args:
            filename (str):
                Name of the file that is always referred to in the background.
            path (str, optional):
                Target path under which the new file will be stored.
            data (np.ndarray, optional):
                If provided, the file will be initialized with these data and shape argument is ignored. Default is
                None.
            shape (tuple, dtype=int, optional):
                If provided, the file will be initialized with an empty array of this shape. This argument will be
                ignored if data argument is provided. Default is None.
            extensions (list of dict, optional):
                List of dictionaries of names (and data or shapes) for each extension that shall be initialized right
                away. Could look like this:
                [{'name': 'VAR', 'data': array},
                 {'name': 'SUB', 'shape': (1024, 1024)}]
                Default is None.
            header (fits.header, optional):
                Header that will be used to initialize the new Primary HDU. Default is None.
            cards (dict, optional):
                Dictionary of cards that will be added to the fits header. Default is None.
            initialize (bool, optional):
                Set to `False` to suppress creating a new file.
            timestamp (bool, optional):
                Set to True to automatically add a time stamp to the file name. Default is False.
            header_card_prefix (str, optional):
                Prefix of header cards. Default is None.
            verbose (bool, optional):
                Set to False to suppress logger output. Default is True.
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
        else:
            raise SpecklepyTypeError('Outfile', 'shape', type(shape), 'tuple')

        if shape is not None and data is not None:
            logger.info("Outfile instance ignores shape input as data is provided for initialization!")

        if isinstance(extensions, dict):
            self.extensions = [extensions]
        elif extensions is None or isinstance(extensions, list):
            self.extensions = extensions
        else:
            raise SpecklepyTypeError('Outfile', 'extensions', type(extensions), 'list of dict')

        if cards is None:
            self.cards = {}
        elif isinstance(cards, dict):
            self.cards = cards
        else:
            raise SpecklepyTypeError('Outfile', 'cards', type(cards), 'dict')

        if header is None or isinstance(header, fits.header.Header):
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

        if initialize:
            self.initialize_file()

    def initialize_file(self, header=None, data=None):

        # Initialize primary HDU
        hdu = fits.PrimaryHDU(header=header)
        for key in self.cards:
            hdu.header.set(self.header_card_prefix + key, self.cards[key])
        if self.shape is not None and data is None:
            hdu.data = np.zeros(self.shape)
        if data is not None:
            hdu.data = data
        if 'DATE' not in hdu.header.keys():
            hdu.header.set('DATE', str(datetime.now()))

        # Create a HDU list with the primary and write to file
        hdu_list = fits.HDUList([hdu])
        if self.verbose:
            logger.info("Initializing file {}".format(self.file_path))
        hdu_list.writeto(self.file_path, overwrite=True)

        if self.extensions is not None:
            for extension in self.extensions:
                if 'name' not in extension:
                    raise ValueError(f"Outfile initialization with extensions failed, because extension did not "
                                     f"contain a name!")
                else:
                    name = extension['name']
                if 'data' in extension:
                    data = extension['data']
                elif 'shape' in extension:
                    data = np.empty(extension['shape'])
                else:
                    data = None
                if 'header' in extension:
                    header = extension['header']
                else:
                    header = None
                self.new_extension(name=name, data=data, header=header)

    @classmethod
    def from_file(cls, filename, path=None):
        return cls(filename=filename, path=path, initialize=False)

    @staticmethod
    def time_stamp():
        """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
        return default_time_stamp()

    @property
    def file_path(self):
        if self.path is None:
            return self.filename
        else:
            return os.path.join(self.path, self.filename)

    @property
    def data(self):
        return fits.getdata(self.file_path)

    @data.setter
    def data(self, data):
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[0].data = data
            hdu_list[0].header.set('UPDATED', str(datetime.now()))
            hdu_list.flush()
        if self.verbose:
            logger.info(f"Updating data in {self.file_path!r}")

    def __getitem__(self, extension):
        return fits.getdata(self.file_path, extension)  # [index]

    def __setitem__(self, extension, data):
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[extension].data = data
            hdu_list[extension].header.set('UPDATED', str(datetime.now()))
            hdu_list.flush()

    def update_frame(self, frame_index, data, extension=None):
        if extension is None:
            extension = 0
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[extension].data[frame_index] = data
            hdu_list[extension].header.set('UPDATED', str(datetime.now()))
            hdu_list.flush()
        # self.data[frame_index] = data

    def new_extension(self, name, data=None, header=None, index=None):
        """Create a new fits extension.

        All arguments are passed to the new HDU in the HDUlist of instances file, i.e. 'self.filename'.

        Args:
            name (str):
                Name of the new extension/ HDU.
            data (array, optional):
                Data array, passed to the new HDU.
            header (fits.header, optional):
                Fits header, passed to the new HDU.
            index (int, optional):
                Index of the extension, before which the new hdu is inserted. None translates into appending the new
                HDU at the end of the HDUlist.
        """

        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu = fits.ImageHDU(data=data, name=name, header=header)
            if index is None:
                hdu_list.append(hdu=hdu)
            else:
                hdu_list.insert(index=index, hdu=hdu)
            hdu_list.flush()

    def update_extension(self, ext_name, data):
        if self.verbose:
            logger.info(f"Updating data in {self.file_path}[{ext_name}]")
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[ext_name].data = data
            hdu_list.flush()

    def has_extension(self, ext_name):
        with fits.open(self.file_path) as hdu_list:
            return ext_name in hdu_list

    @staticmethod
    def extract_frame_shape(header):
        number_axes = header.get('NAXIS')
        if number_axes == 2:
            shape = (header.get('NAXIS1'), header.get('NAXIS2'))
        elif number_axes == 3:
            shape = (header.get('NAXIS2'), header.get('NAXIS3'))
        else:
            raise NotImplementedError(f"Frame shape extraction is not defined for FITS cubes with {number_axes!r} "
                                      f"axes!")
        return shape
