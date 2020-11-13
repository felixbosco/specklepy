from IPython import embed
import numpy as np
import os
from tqdm import trange

from astropy.io import fits
from astropy.stats import sigma_clip, sigma_clipped_stats

from specklepy.logging import logger


def differentiate_cube(files, exposure_time_prefix=None, extension=None, dtype=None, debug=False):

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # Apply default
    if extension is None:
        extension = 0

    # Iterate through files
    for file in files:

        # Make a new copy of the file
        diff_file = 'diff_' + os.path.basename(file)
        logger.info(f"Creating file {diff_file}")
        os.system(f"cp {file} {diff_file}")

        # Load original data and difference
        with fits.open(diff_file, mode='update') as hdu_list:

            # Load input data cube
            cube = hdu_list[extension].data
            if dtype is not None:
                logger.info(f"Casting data to dtype {dtype!r}")
                cube = cube.astype(eval(dtype))

            # Difference the frames along the time axis
            logger.info("Differencing frames...")
            cube = np.diff(cube, axis=0)
            logger.info(f"New cube has shape {cube.shape}")

            # Update exposure time in header
            if exposure_time_prefix is not None:
                exptime = estimate_frame_exposure_times(hdu_list[extension].header, exposure_time_prefix)
                hdu_list[extension].header.set('FEXPTIME', np.around(exptime, 3), 'Frame exposure time (s)')

            # Overwriting data
            logger.info("Storing data to file...")
            hdu_list[extension].data = cube
            hdu_list.flush()

        # Final terminal output
        logger.info(f"Differencing successful for file {diff_file}")


def differentiate_linear_reg(files, exposure_time_prefix=None, extension=None, dtype=None, debug=False):

    # Set defaults
    slope_mask_sigma = 10
    degree = 1

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # Apply default
    if extension is None:
        extension = 0

    # Iterate through files
    for file in files:

        # Make a new copy of the file
        diff_file = 'diff_' + os.path.basename(file)
        logger.info(f"Creating file {diff_file}")
        os.system(f"cp {file} {diff_file}")

        # Load original data and difference
        with fits.open(diff_file, mode='update') as hdu_list:

            # Load input data cube
            header = hdu_list[extension].header
            cube = hdu_list[extension].data
            if dtype is not None:
                logger.info(f"Casting data to dtype {dtype!r}")
                cube = cube.astype(eval(dtype))

            # Update exposure time in header
            if exposure_time_prefix is not None:
                exptime = estimate_frame_exposure_times(header, exposure_time_prefix)
                hdu_list[extension].header.set('FEXPTIME', np.around(exptime, 3), 'Frame exposure time (s)')

            # Initialize arrays
            time_stamps = extract_time_stamps(header=header, common_header_prefix=exposure_time_prefix)
            time_stamps = np.subtract(time_stamps, time_stamps[0])
            slopes = np.empty(cube[0].shape)
            # slopes_var = np.empty(cube[0].shape)
            intercepts = np.empty(cube[0].shape)
            # intercepts_var = np.empty(cube[0].shape)

            # Linear regression
            logger.info("Applying time-wise linear regression through FITS cube...")
            for row in trange(cube.shape[1]):
                for col in range(cube.shape[2]):
                    flux = cube[:, row, col]
                    coefficients = np.polyfit(time_stamps, flux, deg=degree)
                    # coefficients, cov = np.polyfit(time_stamps, flux, deg=degree, cov=True)
                    slopes[row, col] = coefficients[-2]
                    # slopes_var[row, col] = cov[-2, -2]
                    intercepts[row, col] = coefficients[-1]
                    # intercepts_var[row, col] = cov[-1, -1]
            logger.debug("Linear regression finished successfully")

            # Evaluate statistics on the slopes
            clipped_mean, _, _ = sigma_clipped_stats(slopes)
            logger.info("Subtracting offsets and mean slope...")
            cube = np.subtract(cube, intercepts)
            cube = np.swapaxes(np.subtract(np.swapaxes(cube, 0, 2), time_stamps * clipped_mean), 2, 0)

            # Create mask
            logger.info("Creating pixel mask from sigma-clipping the slopes...")
            clipped = sigma_clip(slopes, sigma=slope_mask_sigma)
            mask = clipped.mask.astype(int)
            mask_hdu = fits.ImageHDU(data=mask, name='MASK')
            logger.debug("Pixel mask created successfully")
            logger.debug("Appending pixel mask HDU..")
            hdu_list.append(mask_hdu)

            # Overwriting data
            logger.info("Storing data to file...")
            hdu_list[extension].data = cube
            logger.debug(f"Updating HDU list in file {diff_file!r}")
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
