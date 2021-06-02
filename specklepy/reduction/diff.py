import numpy as np
import os
import sys

from astropy.io import fits
from astropy.stats import sigma_clip

from specklepy.io import FileStream
from specklepy.logging import logger
from specklepy.utils import default_time_stamp, save_eval


def differentiate_cube(files, delta=1, method='direct', exposure_time_prefix=None, extension=None, dtype=None,
                       debug=False):

    # Set logging level
    if debug:
        logger.setLevel('DEBUG')

    # Apply default
    if extension is None:
        extension = 0

    # Iterate through files
    for file in files:

        # Initialize new cube and its file
        diff_cube = DiffCube(in_file=file, extension=extension)
        diff_cube.initialize_file()

        # Differentiate with the requested method
        if method == 'direct':
            diff_cube.differentiate(delta=delta, exposure_time_prefix=exposure_time_prefix, dtype=dtype)
        elif method == 'linreg':
            # diff_cube.differentiate_linear_reg(exposure_time_prefix=exposure_time_prefix, dtype=dtype, debug=debug)
            raise DeprecationWarning("Cube differentiation in linear regression mode is deprecated!")
        else:
            raise ValueError(f"Differentiation method {method!r} not understood!")

        # Final terminal output
        logger.info(f"Differencing successful for file {diff_cube.file!r}")


class DiffCube(object):

    def __init__(self, in_file, extension=None):

        self.input_file = in_file
        self.file = self.default_output_file(in_file)
        self.extension = extension if extension is not None else 0

    @staticmethod
    def default_output_file(input_file):
        return 'diff_' + os.path.basename(input_file)

    def initialize_file(self):
        logger.info(f"Creating file {self.file!r}")
        os.system(f"cp {self.input_file} {self.file}")

    def differentiate(self, delta=1, exposure_time_prefix=None, extension=None, dtype=None, mask_threshold=None,
                      fill_value=0):

        # Update extension attribute if requested
        if extension is not None:
            self.extension = extension

        # Evaluate the data type
        try:
            dtype = eval(dtype)
        except TypeError:
            pass

        # Load original data and difference
        file = FileStream(self.file)
        cube = file.get_data(extension=self.extension, dtype=dtype)
        logger.info("Differencing frames...")
        cube = np.diff(cube[::delta], axis=0)
        logger.info(f"New cube has shape {cube.shape}")

        # Create mask for outliers in the difference cube
        if mask_threshold is not None:
            bpm = sigma_clip(np.mean(cube, axis=0), sigma=mask_threshold, masked=True).mask
            for f, frame in enumerate(cube):
                cube[f] = np.ma.masked_array(frame, mask=bpm).filled(fill_value=fill_value)
            file.new_extension(name='MASK', data=bpm, dtype=np.int16)

        # Update exposure time in header
        if exposure_time_prefix is not None:
            header = file.get_header(extension=self.extension)
            exptime = self.estimate_frame_exposure_times(header, exposure_time_prefix, delta=delta)
            file.set_header('FEXPTIME', np.around(exptime, 3), comment='Frame exposure time (s)',
                            extension=self.extension)
        file.set_header('FDELTA', delta, 'Index-delta of subsequently subtracted frames', extension=extension)

        # Overwriting data
        logger.info("Storing data to file...")
        file.set_data(cube, extension=self.extension)
        file.set_header('PIPELINE', 'SPECKLEPY', extension=self.extension)
        file.set_header('DIFFDATE', default_time_stamp(), extension=self.extension)

        # with fits.open(self.file, mode='update') as hdu_list:
        #
        #     # Load input data cube
        #     cube = hdu_list[self.extension].data
        #     if dtype is not None:
        #         logger.info(f"Casting data to dtype {dtype!r}")
        #         cube = cube.astype(eval(dtype))
        #
        #     # Difference the frames along the time axis
        #     logger.info("Differencing frames...")
        #     cube = np.diff(cube[::delta], axis=0)
        #     logger.info(f"New cube has shape {cube.shape}")
        #
        #     # Create mask for outliers in the difference cube
        #     if mask_threshold is not None:
        #         bpm = sigma_clip(np.mean(cube, axis=0), sigma=mask_threshold, masked=True).mask
        #         for f, frame in enumerate(cube):
        #             cube[f] = np.ma.masked_array(frame, mask=bpm).filled(fill_value=fill_value)
        #         hdu_list.append(fits.ImageHDU(name='MASK', data=bpm.astype(np.int16)))
        #
        #     # Update exposure time in header
        #     if exposure_time_prefix is not None:
        #         exptime = self.estimate_frame_exposure_times(hdu_list[self.extension].header, exposure_time_prefix,
        #                                                      delta=delta)
        #         hdu_list[self.extension].header.set('FEXPTIME', np.around(exptime, 3), 'Frame exposure time (s)')
        #     hdu_list[self.extension].header.set('FDELTA', delta, 'Index-delta of subsequently subtracted frames')
        #
        #     # Overwriting data
        #     logger.info("Storing data to file...")
        #     hdu_list[self.extension].data = cube
        #     hdu_list[self.extension].header.update(pipeline='SPECKLEPY', diffdate=default_time_stamp())
        #     hdu_list.flush()

    # def differentiate_linear_reg(self, exposure_time_prefix=None, extension=None, dtype=None, debug=False):
    #
    #     raise PendingDeprecationWarning("Cube differentiation in linear regression mode is not working yet and will be "
    #                                     "deprecated soon!")
    #
    #     # Update extension attribute if requested
    #     if extension is not None:
    #         self.extension = extension
    #
    #     # Set defaults
    #     slope_mask_sigma = 10
    #     degree = 1
    #
    #     # Load original data and difference
    #     with fits.open(self.file, mode='update') as hdu_list:
    #
    #         # Load input data cube
    #         header = hdu_list[self.extension].header
    #         cube = hdu_list[self.extension].data
    #         if dtype is not None:
    #             logger.info(f"Casting data to dtype {dtype!r}")
    #             cube = cube.astype(eval(dtype))
    #
    #         # Update exposure time in header
    #         if exposure_time_prefix is not None:
    #             exptime = self.estimate_frame_exposure_times(header, exposure_time_prefix)
    #             hdu_list[self.extension].header.set('FEXPTIME', np.around(exptime, 3), 'Frame exposure time (s)')
    #
    #         # Initialize arrays
    #         time_stamps = self.extract_time_stamps(header=header, common_header_prefix=exposure_time_prefix)
    #         time_stamps = np.subtract(time_stamps, time_stamps[0])
    #         slopes = np.empty(cube[0].shape)
    #         # slopes_var = np.empty(cube[0].shape)
    #         intercepts = np.empty(cube[0].shape)
    #         # intercepts_var = np.empty(cube[0].shape)
    #
    #         # Linear regression
    #         logger.info("Applying time-wise linear regression through FITS cube...")
    #         for row in trange(cube.shape[1]):
    #             for col in range(cube.shape[2]):
    #                 flux = cube[:, row, col]
    #                 coefficients = np.polyfit(time_stamps, flux, deg=degree)
    #                 # coefficients, cov = np.polyfit(time_stamps, flux, deg=degree, cov=True)
    #                 slopes[row, col] = coefficients[-2]
    #                 # slopes_var[row, col] = cov[-2, -2]
    #                 intercepts[row, col] = coefficients[-1]
    #                 # intercepts_var[row, col] = cov[-1, -1]
    #         logger.debug("Linear regression finished successfully")
    #
    #         # Evaluate statistics on the slopes
    #         clipped_mean, _, _ = sigma_clipped_stats(slopes)
    #         logger.info("Subtracting offsets and mean slope...")
    #         cube = np.subtract(cube, intercepts)
    #         cube = np.swapaxes(np.subtract(np.swapaxes(cube, 0, 2), time_stamps * clipped_mean), 2, 0)
    #
    #         # Create mask
    #         logger.info("Creating pixel mask from sigma-clipping the slopes...")
    #         clipped = sigma_clip(slopes, sigma=slope_mask_sigma)
    #         mask = clipped.mask.astype(int)
    #         mask_hdu = fits.ImageHDU(data=mask, name='MASK')
    #         logger.debug("Pixel mask created successfully")
    #         logger.debug("Appending pixel mask HDU..")
    #         hdu_list.append(mask_hdu)
    #
    #         # Overwriting data
    #         logger.info("Storing data to file...")
    #         hdu_list[self.extension].data = cube
    #         logger.debug(f"Updating HDU list in file {self.file!r}")
    #         hdu_list.flush()

    def estimate_frame_exposure_times(self, header, common_header_prefix, delta=1):
        """Estimate a mean exposure time per frame

        Args:
            header (fits.Header):
                FITS header to search for the time stamps.
            common_header_prefix (str):
                Common prefix among the header keywords storing time stamp information.
            delta (int, optional):
                Number of subsequent frames to differentiate from another

        Returns:
            mean_exposure_time (float):
                Mean of the exposure times per frame.
        """

        # Differentiate the time stamp values to obtain time deltas
        try:
            diff_times = np.diff(self.extract_time_stamps(header, common_header_prefix)[::delta])
        except IOError as e:
            sys.tracebacklimit = 0
            raise e

        # Report statistics
        logger.info(f"Exposure time is: {np.mean(diff_times):.3f} ({np.std(diff_times):.2e})")

        return np.mean(diff_times)

    @staticmethod
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

        # Check whether the time stamps have been found
        if len(time_stamps) == 0:
            raise IOError(f"Header is missing keywords matching the prefix {common_header_prefix!r}")

        # Remove the zero-th entry
        time_stamps.pop(0)

        return time_stamps
