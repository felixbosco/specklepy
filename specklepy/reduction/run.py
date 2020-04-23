from specklepy.io.filemanager import FileManager
from specklepy.logging import logger
from specklepy.reduction import flat, sky

def all(params):
    # (0) Read file list table
    logger.info("Reading file list ...")
    inFiles = FileManager(params.paths.fileList)
    logger.info('\n' + str(inFiles.table))

    # (1) Initialize reduction files
    # TODO: Implement a data model for the reduction files

    # (2) Flat fielding
    if not params.flat.skip:
        flat_files = inFiles.filter({'OBSTYPE': 'FLAT'})
        if len(flat_files) == 0:
            logger.warning("Did not find any flat field observations. No flat field correction will be applied!")
        else:
            logger.info("Starting flat field correction...")
            master_flat = flat.MasterFlat(flat_files, filename=params.flat.masterFlatFile,
                                          file_path=params.paths.filePath)
            master_flat.combine()

    # (3) Linearization
    # TODO: Implement linearization

    # (4) Sky subtraction
    if not params.sky.skip:
        if params.sky.source == 'default':
            sky_files = inFiles.filter({'OBSTYPE': 'SKY'})
            if len(sky_files) == 0:
                logger.warning("Did not find any sky observations. No sky subtraction will be applied!")
            else:
                logger.info("Starting sky subtraction...")
                logger.info(f"Source of sky background measurement: {params.sky.source}")
                for file in sky_files:
                    sky.get_sky_background(file, path=params.paths.filePath)
        elif params.sky.source in ['image', 'frame']:
            # TODO: Implement sky subtraction from image
            pass