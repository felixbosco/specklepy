import numpy as np
import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logger


class Sequence(object):

    def __init__(self, science_files, sky_files, file_path=None):

        if isinstance(science_files, list):
            self.science_files = science_files
        else:
            raise SpecklepyTypeError('Sequence', 'science_files', type(science_files), 'list')

        if isinstance(sky_files, list):
            self.sky_files = sky_files
        else:
            raise SpecklepyTypeError('Sequence', 'sky_files', type(sky_files), 'list')

        if file_path is None:
            self.file_path = ''
        elif isinstance(file_path, str):
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('Sequence', 'file_path', type(file_path), 'str')

    @property
    def files(self):
        return self.sky_files + self.science_files

    def make_master_sky(self):

        if len(self.sky_files) == 1:
            data = fits.getdata(os.path.join(self.file_path, self.sky_files[0]))
            if data.ndim == 3:
                std = np.std(data, axis=0)
                data = np.mean(data, axis=0)
            else:
                std = np.zeros(data.shape)
            self.master_sky = data
            self.master_sky_std = std

        else:
            for index, file in enumerate(self.sky_files):
                data = fits.getdata(os.path.join(self.file_path, file))
                if data.ndim == 3:
                    std = np.std(data, axis=0)
                    data = np.mean(data, axis=0)
                else:
                    std = np.zeros(data.shape)

                if index == 0:
                    shape = (len(self.sky_files)) + data.shape
                    master_sky = np.zeros(shape)
                    master_sky_uncertainty = np.zeros(shape)

                master_sky[index] = data
                master_sky_uncertainty[index] = std

            # Collapse cubes
            self.master_sky = np.divide(np.sum(np.multiply(master_sky, master_sky_uncertainty), axis=0),
                                        np.sum(master_sky_uncertainty, axis=0))
            self.master_sky_std = np.sqrt(np.sum(np.square(master_sky_uncertainty), axis=0))

    def subtract_master_sky(self, save_to=None, filename_prefix=None):

        if filename_prefix is None:
            filename_prefix = ''

        if not hasattr(self, 'master_sky'):
            self.make_master_sky()

        for index, file in enumerate(self.science_files):
            logger.info('Subtracting sky from file {:2}/{:2}) {}'.format(index+1, len(self.science_files), file))
            with fits.open(os.path.join(self.file_path, file)) as hdulist:
                hdulist[0].data = np.subtract(hdulist[0].data, self.master_sky)
                hdulist[0].header.set('SKYSUB', self.time_stamp())
                logger.info(f"Saving sky subtracted data to file {os.path.join(save_to, filename_prefix + file)}")
                hdulist[0].writeto(os.path.join(save_to, filename_prefix + file))

    # def time_stamp(self):
    #     """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
    #     return datetime.now().strftime('%Y%m%d_%H%M%S')