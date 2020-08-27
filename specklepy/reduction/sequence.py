from IPython import embed
import numpy as np
import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError
from specklepy.logging import logger
from specklepy.utils import combine


class Sequence(object):

    def __init__(self, science_files, sky_files, sky_time_stamps, science_time_stamps, source=None, file_path=None):

        if isinstance(science_files, (list, np.ndarray)):
            self.science_files = science_files
        else:
            raise SpecklepyTypeError('Sequence', 'science_files', type(science_files), 'list')

        if isinstance(sky_files, (list, np.ndarray)):
            self.sky_files = sky_files
        else:
            raise SpecklepyTypeError('Sequence', 'sky_files', type(sky_files), 'list')

        if isinstance(science_time_stamps, (list, np.ndarray)):
            self.science_time_stamps = science_time_stamps
        else:
            raise SpecklepyTypeError('Sequence', 'science_time_stamps', type(science_time_stamps), 'list')

        if isinstance(sky_time_stamps, (list, np.ndarray)):
            self.sky_time_stamps = sky_time_stamps
        else:
            raise SpecklepyTypeError('Sequence', 'sky_time_stamps', type(sky_time_stamps), 'list')

        if source is None:
            if len(sky_files) > 0:
                self.source = 'sky'
            else:
                self.source = 'science'
        elif isinstance(source, str):
            self.source = source
        else:
            raise SpecklepyTypeError('Sequence', 'source', type(source), 'str')

        if file_path is None:
            self.file_path = './'
        elif isinstance(file_path, str):
            self.file_path = file_path
        else:
            raise SpecklepyTypeError('Sequence', 'file_path', type(file_path), 'str')

    def __str__(self):
        s = str(type(self)) + "\n"
        s += 'Science files:\n'
        for file in self.science_files:
            s += "> {}\n".format(file)
        s += 'Sky files:\n'
        for file in self.sky_files:
            s += "> {}\n".format(file)
        return s

    @property
    def files(self):
        return np.concatenate((self.sky_files, self.science_files))

    @property
    def time_stamps(self):
        return np.concatenate((self.sky_time_stamps, self.science_time_stamps))

    @property
    def n_sky(self):
        return len(self.sky_files)

    @property
    def n_science(self):
        return len(self.science_files)

    def time_stamps_to_offset_seconds(self):

        # Check whether time stamps free from offset (equivalent to float type)
        if isinstance(self.time_stamps[0], float):
            return 0

        # Find first time stamp
        time_stamps = self.time_stamps
        time_stamps.sort()
        start = time_stamps[0]

        # Remove time offset
        self.sky_time_stamps = combine.time_difference(start, self.sky_time_stamps)
        self.science_time_stamps = combine.time_difference(start, self.science_time_stamps)

    def compute_weights(self, time_scale=300):
        # Assert that times are offset seconds
        self.time_stamps_to_offset_seconds()

        # Compute weights from time differences
        weights = np.zeros((self.n_science, self.n_sky))
        for ii in range(self.n_science):
            delta_t = self.science_time_stamps[ii] - self.sky_time_stamps
            weights[ii] = np.exp(-1/2 * np.square(np.divide(delta_t, time_scale)))

        # # Transform time offsets into weights
        # weights = np.square(weights)

        # Normalize the weights
        for ii in range(self.n_science):
            weights[ii] /= np.sum(weights[ii])

        return weights

    # def make_master_sky(self):
    #
    #     if len(self.sky_files) == 1:
    #         data = fits.getdata(os.path.join(self.file_path, self.sky_files[0]))
    #         if data.ndim == 3:
    #             std = np.std(data, axis=0)
    #             data = np.mean(data, axis=0)
    #         else:
    #             std = np.zeros(data.shape)
    #         self.master_sky = data
    #         self.master_sky_std = std
    #
    #     else:
    #         for index, file in enumerate(self.sky_files):
    #             data = fits.getdata(os.path.join(self.file_path, file))
    #             if data.ndim == 3:
    #                 std = np.std(data, axis=0)
    #                 data = np.mean(data, axis=0)
    #             else:
    #                 std = np.zeros(data.shape)
    #
    #             if index == 0:
    #                 shape = (len(self.sky_files)) + data.shape
    #                 master_sky = np.zeros(shape)
    #                 master_sky_uncertainty = np.zeros(shape)
    #
    #             master_sky[index] = data
    #             master_sky_uncertainty[index] = std
    #
    #         # Collapse cubes
    #         self.master_sky = np.divide(np.sum(np.multiply(master_sky, master_sky_uncertainty), axis=0),
    #                                     np.sum(master_sky_uncertainty, axis=0))
    #         self.master_sky_std = np.sqrt(np.sum(np.square(master_sky_uncertainty), axis=0))
    #
    # def subtract_master_sky(self, save_to=None, filename_prefix=None):
    #
    #     if filename_prefix is None:
    #         filename_prefix = ''
    #
    #     if not hasattr(self, 'master_sky'):
    #         self.make_master_sky()
    #
    #     for index, file in enumerate(self.science_files):
    #         logger.info('Subtracting sky from file {:2}/{:2}) {}'.format(index+1, len(self.science_files), file))
    #         with fits.open(os.path.join(self.file_path, file)) as hdulist:
    #             hdulist[0].data = np.subtract(hdulist[0].data, self.master_sky)
    #             hdulist[0].header.set('SKYSUB', self.time_stamp())
    #             logger.info(f"Saving sky subtracted data to file {os.path.join(save_to, filename_prefix + file)}")
    #             hdulist[0].writeto(os.path.join(save_to, filename_prefix + file))

    # def time_stamp(self):
    #     """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
    #     return datetime.now().strftime('%Y%m%d_%H%M%S')
