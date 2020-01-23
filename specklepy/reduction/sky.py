import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

from specklepy.logging import logging
from specklepy.exceptions import SpecklepyTypeError
from specklepy.utils.plot import imshow



def identify_sequences(file_list, ordered_sky_subtraction=True):

    if not isinstance(file_list, Table):
        raise SpecklepyTypeError('identify_sequences', 'file_list', type(file_list), 'astropy.table.Table')

    sequences = []
    is_science_file = file_list['OBSTYPE'] == 'SCIENCE'
    science_files = list(file_list['FILE'][is_science_file])
    is_sky_file = file_list['OBSTYPE'] == 'SKY'
    sky_files = list(file_list['FILE'][is_science_file])
    if ordered_sky_subtraction:
        pass
    else:
        sequences.append(Sequence(science_files=science_files, sky_files=sky_files))

    logging.info("Identified {} squence(s)...".format(len(sequences)))
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
