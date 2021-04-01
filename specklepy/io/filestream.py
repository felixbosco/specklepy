import numpy as np
import os

from astropy.io import fits

from specklepy.io.fits import get_data
from specklepy.logging import logger
from specklepy.utils.array import frame_number, frame_shape
from specklepy.utils.time import default_time_stamp


class FileStream(object):

    variance_extension = 'VAR'
    mask_extension = 'MASK'
    header_card_prefix = 'HIERARCH SPECKLEPY '

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
                hdu.header.set(self.header_card_prefix + card, value=value)
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
            hdu_list.flush()

    def has_extension(self, extension):
        with fits.open(self.file_path) as hdu_list:
            return extension in hdu_list

    def new_extension(self, name, data=None, header=None, index=None):
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
        """

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
