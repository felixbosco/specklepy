from datetime import datetime
import os

from astropy.io import fits

from specklepy.io.filearchive import FileArchive
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
    in_files = FileArchive(params['PATHS']['fileList'],
                           in_dir=params['PATHS']['filePath'],
                           out_dir=params['PATHS']['outDir'])
    logger.info('\n' + str(in_files.table))

    # (1) Initialize reduction files
    if not os.path.isdir(params['PATHS']['outDir']):
        os.makedirs(params['PATHS']['outDir'])
    product_files = in_files.initialize_product_files()

    # (2) Flat fielding
    if 'skip' in params['FLAT'] and params['FLAT']['skip']:
        logger.info('Skipping flat fielding as requested from parameter file...')
    else:
        flat_files = in_files.get_flats()
        if len(flat_files) == 0:
            logger.warning("Did not find any flat field observations. No flat field correction will be applied!")
        else:
            logger.info("Starting flat field correction...")
            master_flat = flat.MasterFlat(flat_files, file_name=params['FLAT']['masterFlatFile'],
                                          file_path=params['PATHS']['filePath'])
            master_flat.combine()
            master_flat.run_correction(file_list=product_files, file_path=None)

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
