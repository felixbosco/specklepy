import numpy as np
import os

from astropy.io import fits

from specklepy.logging import logger


def differentiate_cube(files, exposure_time_prefix=None, extension=None, dtype=None, debug=False):

    hdu = 0 if extension is None else extension

    # Iterate through files
    for file in files:

        # Make a new copy of the file
        diff_file = 'diff_' + os.path.basename(file)
        logger.info(f"Creating file {diff_file}")
        os.system(f"cp {file} {diff_file}")

        # Load original data and difference
        with fits.open(diff_file, mode='update') as hdu_list:

            # Load input data cube
            cube = hdu_list[hdu].data
            if dtype is not None:
                cube = cube.astype(eval(dtype))

            # Difference the frames along the time axis
            logger.info("Differencing frames...")
            cube = np.diff(cube, axis=0)
            logger.info(f"New cube has shape {cube.shape}")

            # Update exposure time in header
            if exposure_time_prefix is not None:
                exptime = estimate_frame_exposure_times(hdu_list[hdu].header, exposure_time_prefix)
                hdu_list[hdu].header.set('FEXPTIME', np.around(exptime, 3), 'Frame exposure time (s)')

            # Overwriting data
            logger.info("Storing data to file...")
            hdu_list[hdu].data = cube
            hdu_list.flush()

        # Final terminal output
        logger.info(f"Differencing successful for file {diff_file}")


def estimate_frame_exposure_times(header, common_header_prefix):
    """Estimate a mean exposure time per frame

    Args:
        header (fits.Header):
            FITS header to search for the time stamps.
        common_header_prefix (str):
            Common prefix among the header keywords storing time stamp information.

    Returns:
        mean_exposure_time (float):
            Mean of the exposure times per frame.
    """

    # Differentiate the time stamp values to obtain time deltas
    diff_times = np.diff(extract_time_stamps(header, common_header_prefix))

    # Report statistics
    logger.info(f"Exposure time is: {np.mean(diff_times):.3f} ({np.std(diff_times):.2e})")

    return np.mean(diff_times)


def extract_time_stamps(header, common_header_prefix):
    """Extract the time stamp values from a FITS header.

    Args:
        header (fits.Header):
            FITS header to search for the time stamps.
        common_header_prefix (str):
            Common prefix among the header keywords storing time stamp information.

    Returns:
        time_stamps (list):
            List of time stamp values, typically in units of seconds.
    """

    # Initialize list
    time_stamps = []

    # Iterate through header cards
    for keyword, value in header.items():
        if common_header_prefix in keyword:
            time_stamps.append(value)

    # Remove the zero-th entry
    time_stamps.pop(0)

    return time_stamps
