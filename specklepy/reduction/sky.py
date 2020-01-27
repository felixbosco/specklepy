import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError
from specklepy.utils.plot import imshow



def identify_sequences(file_list, ignore_time_stamps=False):

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

    if not isinstance(file_list, Table):
        raise SpecklepyTypeError('identify_sequences', 'file_list', type(file_list), 'astropy.table.Table')

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

        # Assigning sky files to science files
        if ignore_time_stamps:
            sequences.append(Sequence(science_files=science_files, sky_files=sky_files))
        else:
            pass

    logging.info("Identified {} sequence(s)...".format(len(sequences)))
    return sequences


class Sequence(object):

    def __init__(self, science_files, sky_files):

        if isinstance(science_files, list):
            self.science_files = science_files
        else:
            raise SpecklepyTypeError('Sequence', 'science_files', type(science_files), 'list')

        if isinstance(sky_files, list):
            self.sky_files = sky_files
        else:
            raise SpecklepyTypeError('Sequence', 'sky_files', type(sky_files), 'list')

    @property
    def files(self):
        return self.sky_files + self.science_files

    def make_master_sky(self):
        self.master_sky = 0

    def subtract_master_sky(self):
        pass




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
