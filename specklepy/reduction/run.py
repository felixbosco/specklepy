from IPython import embed
import matplotlib.pyplot as plt
import numpy as np

from specklepy.exceptions import SpecklepyValueError
from specklepy.io.filemanager import FileManager
from specklepy.logging import logger
from specklepy.reduction import flat, sky

# TODO: Split this function into the parts and sort into the other modules


def full_reduction(params, debug=False):
    """Execute a full reduction following the parameters in the `params` dictionary.

    Args:
        params (dict):
            Dictionary with all the settings for reduction
        debug (bool, optional):
            Show debugging information
    """

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # (0) Read file list table
    logger.info("Reading file list ...")
    in_files = FileManager(params['PATHS']['fileList'])
    logger.info('\n' + str(in_files.table))

    # (1) Initialize reduction files
    # TODO: Implement a data model for the reduction files

    # (2) Flat fielding
    if 'skip' in params['FLAT'] and params['FLAT']['skip']:
        logger.info('Skipping flat fielding as requested from parameter file...')
    else:
        flat_files = in_files.filter({'OBSTYPE': 'FLAT'})
        if len(flat_files) == 0:
            logger.warning("Did not find any flat field observations. No flat field correction will be applied!")
        else:
            logger.info("Starting flat field correction...")
            master_flat = flat.MasterFlat(flat_files, filename=params['FLAT']['masterFlatFile'],
                                          file_path=params['PATHS']['filePath'])
            master_flat.combine()

    # (3) Linearization
    # TODO: Implement linearization

    # (4) Sky subtraction
    if 'skip' in params['SKY'] and params['SKY']['skip']:
        logger.info('Skipping sky background subtraction as requested from parameter file...')
    else:
        logger.info("Starting sky subtraction...")
        sky.subtract_sky_background(**params['SKY'], in_files=in_files, file_path=params['PATHS']['filePath'])

    # Close reduction
    logger.info("Reduction finished...")
