from datetime import datetime
import os

from astropy.io import fits

from specklepy.io.filemanager import FileManager
from specklepy.logging import logger
from specklepy.reduction import flat, sky


def full_reduction(params, debug=False):
    """Execute a full reduction following the parameters in the `params` dictionary.

    TODO: Split this function into the parts and sort into the other modules

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
    if not os.path.isdir(params['PATHS']['outDir']):
        os.makedirs(params['PATHS']['outDir'])
    for file in in_files.filter({'OBSTYPE': 'SCIENCE'}):
        src = os.path.join(params['PATHS']['filePath'], file)
        dest = os.path.join(params['PATHS']['outDir'], 'r'+file)
        logger.info(f"Initializing data product file {dest}")
        os.system(f"cp {src} {dest}")
        with fits.open(dest) as hdu_list:
            hdu_list[0].header.set('PIPELINE', 'Specklepy')
            hdu_list[0].header.set('REDUCED', datetime.now().strftime('%Y-%m-%d_%H:%M:%S'))

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
        try:
            sky.subtract_sky_background(**params['SKY'], in_files=in_files, file_path=params['PATHS']['filePath'])
        except RuntimeError as e:
            raise RuntimeWarning(e)

    # Close reduction
    logger.info("Reduction finished...")
