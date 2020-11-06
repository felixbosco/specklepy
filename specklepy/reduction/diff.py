import numpy as np
import os

from astropy.io import fits

from specklepy.logging import logger


def differentiate_cube(files, exposure_time_prefix=None, extension=None, debug=False):

    hdu = 0 if extension is None else extension

    # Iterate through files
    for file in files:

        # Make a new copy of the file
        diff_file = 'diff_' + os.path.basename(file)
        logger.info(f"Creating file {diff_file}")
        os.system(f"cp {file} {diff_file}")

        # Load original data and difference
        with fits.open(diff_file, mode='update') as hdu_list:

            # Load original cube
            cube = hdu_list[hdu].data.astype(int)

            # Difference the frames along the time axis
            logger.info("Differencing frames...")
            cube = np.diff(cube, axis=0)
            logger.info(f"New cube has shape {cube.shape}")

            # Update exposure time in header
            if exposure_time_prefix is not None:
                exptime = estimate_frame_exposure_times(hdu_list[hdu].header, exposure_time_prefix)
                hdu_list[hdu].header.set('EXPTIME', exptime)

            # Overwriting data
            logger.info("Storing data to file...")
            hdu_list[hdu].data = cube
            hdu_list.flush()

        # Final terminal output
        logger.info(f"Differencing successful for file {diff_file}")


def estimate_frame_exposure_times(header, common_header_prefix):

    times = []
    for card in header.cards:
        if common_header_prefix in card.keyword:
            times.append(card.value)
    times.pop(0)

    diff_times = np.diff(times)

    # Report statistics
    logger.info(f"Exposure time is: {np.mean(diff_times):.3e} ({np.std(diff_times):.2e})")

    return np.mean(diff_times)
