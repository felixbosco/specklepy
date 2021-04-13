import numpy as np
import os

from astropy.io import fits

from specklepy.io.fits import get_data, get_frame_number
from specklepy.logging import logger
from specklepy.utils.array import frame_number, frame_shape
from specklepy.utils.time import default_time_stamp


class FileStream(object):

    variance_extension = 'VAR'
    mask_extension = 'MASK'
    header_card_prefix = 'HIERARCH SPECKLEPY'

    def __init__(self, file_name, path=None):
        self.file_name = file_name
        self.path = path

    @property
    def file_path(self):
        if self.path is None:
            return self.file_name
        else:
            return os.path.join(self.path, self.file_name)

    @property
    def frame_number(self):
        try:
            data = get_data(self.file_path)
        except FileNotFoundError as e:
            raise e
        return frame_number(data)

    @property
    def frame_shape(self):
        try:
            data = get_data(self.file_path)
        except FileNotFoundError as e:
            raise e
        return frame_shape(data)

    def exists(self):
        return os.path.exists(self.file_path)

    def initialize(self, data=None, shape=None, header=None, cards=None, overwrite=True):

        # Initialize data
        if data is None and shape is not None:
            data = np.empty(shape)

        # Initialize primary HDU
        hdu = fits.PrimaryHDU(data=data, header=header)

        # Add cards to header
        if isinstance(cards, dict):
            for (card, value) in cards.items():
                hdu.header.set(f"{self.header_card_prefix} {card}", value=value)
        hdu.header.set('DATE', default_time_stamp())

        # Store new HDU to file
        logger.info(f"Initializing file {self.file_path!r}")
        hdu.writeto(self.file_path, overwrite=overwrite)

    def get_data(self, extension=None, squeeze=False, dtype=None):
        return get_data(self.file_path, extension=extension, squeeze=squeeze, dtype=dtype)

    def set_data(self, data, extension=None, dtype=None):

        # Cast data to requested data type
        if dtype is not None:
            data = data.astype(dtype=dtype)

        # Update data in file
        with fits.open(name=self.file_path, mode='update') as hdu_list:
            try:
                hdu_list[extension].data = data
            except TypeError as e:
                hdu_list[0].data = data
            except KeyError:
                self.new_extension(name=extension, data=data)
            hdu_list.flush()

    def get_frame(self, frame_index, extension=None):
        return get_data(self.file_path, extension=extension)[frame_index]

    def has_extension(self, extension):
        with fits.open(self.file_path) as hdu_list:
            return extension in hdu_list

    def new_extension(self, name, data=None, header=None, index=None, dtype=None):
        """Create a new FITS file extension.

        All arguments are parsed to the new HDU in the HDUList of instances file, i.e. `self.file_name`.

        Args:
            name (str):
                Name of the new extension/ HDU.
            data (array, optional):
                Data array, passed to the new HDU.
            header (fits.header, optional):
                Fits header, passed to the new HDU.
            index (int, optional):
                Index of the extension, before which the new hdu is inserted. None translates into appending the new
                HDU at the end of the HDUList.
            dtype (type, optional):
                Data type, to which the data are casted prior to storing into the new HDU.
        """

        # Cast data to requested data type
        if dtype is not None:
            data = data.astype(dtype=dtype)

        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu = fits.ImageHDU(data=data, name=name, header=header)
            if index is None:
                hdu_list.append(hdu=hdu)
            else:
                hdu_list.insert(index=index, hdu=hdu)
            hdu_list.flush()

    def update_frame(self, frame_index, data, extension=None):
        if extension is None:
            extension = 0
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[extension].data[frame_index] = data
            hdu_list[extension].header.set('UPDATED', default_time_stamp())
            hdu_list.flush()

    def set_header(self, key, value, comment=None, extension=None):
        if extension is None:
            extension = 0
        with fits.open(self.file_path, mode='update') as hdu_list:
            hdu_list[extension].header.set(key, value, comment)
            hdu_list.flush()

    def insert_header_cards(self, cards, extension=None):
        with fits.open(self.file_path, mode='update') as hdu_list:

            # Load HDU
            try:
                hdu = hdu_list[extension]
            except KeyError:
                hdu = hdu_list[0]

            # Insert header cards
            for key, value in cards.items():
                try:
                    hdu.header.set(key, value)
                except ValueError:
                    hdu.header.set(key, value[0], value[1])

            # Store to file
            hdu_list.flush()

    def build_reconstruction_file_header_cards(self, files, path=None, algorithm=None, card_prefix=None,
                                               insert=True, extension=None):

        # Update the card prefix
        if card_prefix is not None:
            self.header_card_prefix = card_prefix

        # Initialize and fill `cards` dictionary
        cards = {}
        if algorithm is not None:
            cards[f"{self.header_card_prefix} ALGORITHM"] = algorithm
        for index, file in enumerate(files):
            cards[f"{self.header_card_prefix} SOURCE FILE{index:04} NAME"] = os.path.basename(file)
            cards[f"{self.header_card_prefix} SOURCE FILE{index:04} FRAMES"] = get_frame_number(file, path=path)

        # Insert into FITS header
        if insert:
            self.insert_header_cards(cards=cards, extension=extension)

        return cards
