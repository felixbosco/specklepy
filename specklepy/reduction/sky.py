import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from datetime import datetime

from specklepy.io.outfile import Outfile
from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError
from specklepy.utils.plot import imshow



def identify_sequences(file_list, file_path='', ignore_time_stamps=False):

    """Identify observational sequences for sky subtraction.

    Typically, observations are organized in sequences such as 'object-sky-object' for optimizing the required
    telescope time. This function now identifies (blocks of) sky files and afterwards allocates science or object files
    to these blocks by choosing estimating the minimum time delay to one of the sky files.

    Args:
        file_list (astropy.table.Table):
        shape (tuple, dtype=int, optional):
        cards (dict, optional):
        timestamp (bool, optional):
            Set to True to automatically add a time stamp to the file name.
            Default is False.
        hprefix (str, optional):
            Prefix of header cards. Default is None.
    """

    # Check input parameters
    if not isinstance(file_list, Table):
        raise SpecklepyTypeError('identify_sequences', 'file_list', type(file_list), 'astropy.table.Table')

    # Iterate over observational setups and assign files to sequences
    sequences = []
    observation_setups = np.unique(file_list['Setup'])
    for setup in observation_setups:

        # Check for setup
        is_in_setup = file_list['Setup'] == setup

        # Check for science or sky type
        is_science_file = file_list['OBSTYPE'] == 'SCIENCE'
        is_science_file &= is_in_setup
        science_files = list(file_list['FILE'][is_science_file])
        is_sky_file = file_list['OBSTYPE'] == 'SKY'
        is_sky_file &= is_in_setup
        sky_files = list(file_list['FILE'][is_sky_file])

        # Check whether there are sky files
        if len(sky_files) is 0:
            raise RuntimeError("There are no sky files in the list!")

        # Assigning sky files to science files
        if ignore_time_stamps:
            sequences.append(Sequence(science_files=science_files, sky_files=sky_files, file_path=file_path))
        else:
            pass

    logging.info("Identified {} sequence(s)...".format(len(sequences)))
    return sequences



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


    def __call__(self, outpath, filename_prefix=None):
        pass



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
            self.master_sky = np.divide(np.sum(np.multiply(master_sky, master_sky_uncertainty), axis=0), np.sum(master_sky_uncertainty, axis=0))
            self.master_sky_std = np.sqrt(np.sum(np.square(master_sky_uncertainty), axis=0))


    def subtract_master_sky(self, saveto=None, filename_prefix=None):

        if filename_prefix is None:
            filename_prefix = ''

        if not hasattr(self, 'master_sky'):
            self.make_master_sky()

        for index, file in enumerate(self.science_files):
            logging.info('Subtracting sky from file {:2}/{:2}) {}'.format(index+1, len(self.science_files), file))
            with fits.open(os.path.join(self.file_path, file)) as hdulist:
                hdulist[0].data = np.subtract(hdulist[0].data, self.master_sky)
                hdulist[0].header.set('SKYSUB', self.time_stamp())
                logging.info('Saving sky subtracted data to file {}'.format(os.path.join(saveto, filename_prefix + file)))
                hdulist[0].writeto(os.path.join(saveto, filename_prefix + file))


    def time_stamp(self):
        """Return a time stamp str of format 'YYYYMMDD_HHMMSS'."""
        return datetime.now().strftime('%Y%m%d_%H%M%S')



def subtract_scalar_background(files, params, prefix=None, debug=False):
    """Estimate and subtract a scalar background."""

    if not isinstance(files, (list, np.ndarray)):
        raise SpecklepyTypeError('substract_scalar_background', argtype=type(files), argname='files', expected='list')
    else:
        if len(files) == 0:
            raise RuntimeError("Sky subtraction received an empty list of files!")

    logging.info("Estimating scalar background and subtract...")
    for file_index, file in enumerate(files):

        image, header = fits.getdata(os.path.join(params.paths.filePath, file), header=True)

        # Update header for and initialize the outfile
        header.set('PIPELINE', 'SPECKLEPY')
        header.set('SKYCORR', str(datetime.now()))
        corrected_file = prefix + file
        corrected_file = os.path.join(params.paths.filePath, corrected_file)
        outfile = Outfile(filename=corrected_file, header=header, shape=image.shape)

        # Estimate scalar background and uncertainties, and subtract
        if image.ndim == 2:
            mean, median, std = sigma_clipped_stats(image, sigma=params.sky.backgroundSigmaClip)
            outfile.data = image - mean
            image_var = np.ones(image.shape) * np.square(std)
            outfile.new_extension(name='VAR', data=image_var)
        elif image.ndim == 3:
            means, medians, stds = sigma_clipped_stats(image, sigma=params.sky.backgroundSigmaClip, axis=(1, 2))
            logging.info(f"Sigma clipped stats:\t{np.mean(means):.2f} +- {np.mean(stds):.2f}")
            outfile.new_extension(name='VAR', data=np.zeros(image.shape))
            tmp_frame = np.ones(image[0].shape)

            for frame_index, frame in enumerate(image):
                print(f"\r\tUpdating frame {frame_index+1:3}...", end='')
                outfile.update_frame(frame_index=frame_index, data=np.subtract(frame, means[frame_index]))
                outfile.update_frame(frame_index=frame_index, data=tmp_frame * np.square(stds[frame_index]), extension='VAR')
            print()
        else:
            raise RuntimeError(f"Images are supposed to have 2 or 3 dimensions but this one has {image.ndim}!")

    logging.info("Scalar background subtraction complete!")



def subtract(params, mode='constant', debug=False):
    logging.info("Subtracting sky from files\n\t{}".format(params.scienceFiles))
    logging.info("Sky files are\n\t{}".format(params.skyFiles))

    if isinstance(params.skyFiles, str):
        skyFiles = [params.skyFiles]
    if len(skyFiles) > 1:
        raise ValueError("sky.subtract does not handle more than one file yet!")

    skyFile = skyFiles[0]
    if mode == 'constant':
        median = np.median(fits.getdata(skyFile))
    elif mode == 'image':
        median = np.median(fits.getdata(skyFile), axis=0)
        imshow(median)

    for scienceFile in params.scienceFiles:
        # Create a copy
        skySubtractedScienceFile = os.path.basename(scienceFile)
        skySubtractedScienceFile = os.path.join(params.tmpDir, "s{}".format(skySubtractedScienceFile))
        logging.info("Subtrcting sky and saving tmp file:\n\t{}".format(skySubtractedScienceFile))
        os.system("cp {} {}".format(scienceFile, skySubtractedScienceFile))

        with fits.open(skySubtractedScienceFile, mode='update') as hdulist:
            for frame_index, frame in enumerate(hdulist[0].data):
                hdulist[0].data[frame_index] -= median
                hdulist.flush()
