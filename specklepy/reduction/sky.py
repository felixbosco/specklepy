import os
import numpy as np
from astropy.io import fits

from specklepy.logging import logging
from specklepy.utils.plot import imshow



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
