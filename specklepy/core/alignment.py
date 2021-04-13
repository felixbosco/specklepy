import numpy as np
import os

from specklepy.core.sourceextraction import SourceExtractor
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.fits import get_data, get_frame_shape
from specklepy.logging import logger
from specklepy.plotting.utils import imshow
from specklepy.utils.array import peak_index


class ShiftEstimator(object):

    def __init__(self, reference_file=None, reference_image=None, in_dir=None):

        # Store reference image
        self.reference_file = reference_file
        self._reference_image = reference_image
        self.in_dir = in_dir

        # Initialize shifts
        self.shifts = None

    def __repr__(self):
        return f"ShiftEstimator ({self.reference_image_path})"

    @property
    def reference_image_path(self):
        if self.in_dir is None:
            return self.reference_file
        else:
            return os.path.join(self.in_dir, self.reference_file)

    @property
    def reference_image(self):
        if self._reference_image is None:
            self._reference_image = self.load_collapsed(self.reference_image_path)
        return self._reference_image

    @property
    def image_shape(self):
        return self._reference_image.shape

    @staticmethod
    def load_collapsed(file_name, path=None):
        """Load and collapse a data cube by integrating over time axis if not 2D image

        Arguments:
            path (str):
                Path to the file.

        Returns:
            data (np.ndarray):
                2D array, obtained from the input file.
        """

        logger.debug(f"Loading data from file {path!r}")
        data = get_data(file_name=file_name, path=path)

        if data.ndim == 3:
            data = np.sum(data, axis=0)

        return data

    def estimate_shifts(self, file_names, mode='correlation', in_dir=None, reference_file_index=None, debug=False):

        # Transform str-type file names into a list
        if isinstance(file_names, str):
            file_names = [file_names]

        # Update reference file and input directory if provided
        if reference_file_index is not None:
            self.reference_file = file_names[reference_file_index]
        if in_dir is not None:
            self.in_dir = in_dir

        # Reset stored shifts
        self.shifts = []

        # Be lazy and do not compute anything, if only one file is provided
        if len(file_names) == 1 and file_names[0] == os.path.basename(self.reference_file):
            logger.info("Only one data cube is provided, nothing to align.")
            self.shifts = [(0, 0)]
            return self.shifts

        # Estimate shifts in the requested modes
        if mode == 'peak':
            reference_peak = peak_index(self.reference_image)

            # Iterate through files
            for file in file_names:
                if file == self.reference_file:
                    shift = (0, 0)
                else:
                    image = self.load_collapsed(file, path=self.in_dir)
                    peak = peak_index(image)
                    shift = reference_peak[0] - peak[0], reference_peak[1] - peak[1]

                logger.info(f"Estimated shift {shift} for file {file!r}")
                self.shifts.append(shift)

        elif mode == 'correlation':
            # Get the Fourier transformed reference image for cross-correlation
            f_reference_image = np.fft.fft2(self.reference_image)

            # Iterate through files
            for file in file_names:
                if file == self.reference_file:
                    shift = (0, 0)
                else:
                    # Load and Fourier transform the image
                    image = self.load_collapsed(file, path=self.in_dir)
                    f_image = np.conjugate(np.fft.fft2(image))

                    # Compute the 2-dimensional correlation
                    correlation = np.fft.fftshift(np.fft.ifft2(np.multiply(f_reference_image, f_image)))
                    if debug:
                        imshow(np.abs(correlation), title='FFT shifted correlation')

                    # Derive the shift from the correlation
                    correlation_peak = peak_index(correlation)
                    shift = tuple(x - correlation.shape[i] // 2 for i, x in enumerate(correlation_peak))

                logger.info(f"Estimated shift {shift} for file {file!r}")
                self.shifts.append(shift)

        elif mode == 'sources':
            # Initialize SourceExtractor
            extractor = SourceExtractor(fwhm=10)
            extractor.initialize_image(source=self.reference_image_path)
            extractor.initialize_star_finder()
            extractor.find_sources()

            # Extract the position of the reference star in the reference image
            message = "Select the sources that shall be used for correlating the source lists!"
            reference_stars = extractor.select(message=message)
            if reference_stars is None or len(reference_stars) == 0:
                logger.warning("No sources have been selected! Please repeat the source selection!")
                reference_stars = extractor.select(message=message)
            logger.info(f"The sources reside at\n{reference_stars}\nin the reference image")

            # Iterate through files
            for file in file_names:
                if file == self.reference_file:
                    shift = (0, 0)
                else:
                    # Load the image into the source extractor
                    image = self.load_collapsed(file, path=self.in_dir)
                    extractor.initialize_image(source=image)
                    extractor.initialize_star_finder()
                    extractor.find_sources()

                    # Derive shift from cross matching the identified sources to the reference stars
                    shift = extractor.cross_match_table(reference_stars)
                    shift = tuple([int(round(x)) for x in shift])

                logger.info(f"Estimated shift {shift} for file {file!r}")
                self.shifts.append(shift)

        else:
            raise ValueError(f"Mode {mode!r} for estimating shifts is not understood! Choose from 'peak', "
                             f"'correlation' and 'sources'!")

        # Return
        return self.shifts


def derive_pad_vectors(shifts, cube_mode=False):
    """Computes padding vectors from the relative shifts between files.

    Args:
        shifts (list or np.ndarray):
            Shifts between files, relative to a reference image. See get_shifts function for details.
        cube_mode (bool, optional):
            If image is a cube, the estimated pad vectors will obtain pad_vector entries of (0, 0) for the zeroth axis.
            Default is False.

    Returns:
        pad_vectors (list):
            List of padding vectors for each shift in shifts.
        reference_image_pad_vector (list, optional):
            Pad vector of the reference image, which is needed for pad_array() in 'same' mode.
    """

    # Check input parameters
    if isinstance(shifts, (list, np.ndarray)):
        pass
    else:
        raise SpecklepyTypeError('get_pad_vectors()', argname='shifts', argtype=type(shifts), expected='list')
    if not isinstance(cube_mode, bool):
        raise SpecklepyTypeError('get_pad_vectors()', argname='cube_mode', argtype=type(cube_mode), expected='bool')

    # Initialize list
    pad_vectors = []

    # Get extreme points for 'full' padding
    y_max, x_max = np.max(np.array(shifts), axis=0)
    y_min, x_min = np.min(np.array(shifts), axis=0)

    # Iterate over Shifts
    for shift in shifts:

        if cube_mode:
            pad_vector = [(0, 0)]
        else:
            pad_vector = []

        pad_vector.append((shift[0] - y_min, y_max - shift[0]))
        pad_vector.append((shift[1] - x_min, x_max - shift[1]))

        pad_vectors.append(pad_vector)

    # In 'same' mode, pad_array needs also the pad vector of the reference image
    reference_image_pad_vector = [(np.abs(y_min), np.abs(y_max)), (np.abs(x_min), np.abs(x_max))]

    return pad_vectors, reference_image_pad_vector


def pad_array(array, pad_vector, mode='same', reference_image_pad_vector=None):
    """Pads an array according to the pad_vector and crops the image given the
    mode.

    Pads an array with zeros to match a desired field size. Intermediately, it always creates a 'full' image and only
    in 'same' mode it crops the edges such that the returned array covers only the field of the reference image.

    Args:
        array (np.ndarray):
            Input array that shall be padded to match the 'full' or 'same' fields.
        pad_vector (list):
            List of padding vectors, as obtained from get_pad_vectors().
        mode (str, optional):
            Define the size of the output image as 'same' to the reference image or expanding to include the 'full'
            covered field.
        reference_image_pad_vector (tuple or list, optional):
            Used in `same` mode to estimate the position of the reference image and crop beyond.

    Returns:
        padded (np.ndarray):
            Padded array, matching the field of the reference image in 'same'
            mode, or the complete field in 'full' mode.
    """

    # Check input parameters
    if isinstance(array, np.ndarray):
        if array.ndim not in [2, 3]:
            raise SpecklepyValueError('pad_array()', argname='array.ndim', argvalue=array.ndim, expected='2 or 3')
    else:
        raise SpecklepyTypeError('pad_array()', argname='array', argtype=type(array), expected='np.ndarray')

    #
    padded = np.pad(array, pad_vector, mode='constant')

    # Crop the image according to the desired mode
    if mode == 'same':
        # Take reference pad vector and adapt to correctly limit the image
        _r = reference_image_pad_vector
        # Pick only those pixels, covered by the reference image
        if array.ndim == 2:
            padded = padded[_r[0][0]: _adapt_max_coordinate(_r[0][1]), _r[1][0]: _adapt_max_coordinate(_r[1][1])]
        else:
            padded = padded[:, _r[0][0]: _adapt_max_coordinate(_r[0][1]), _r[1][0]: _adapt_max_coordinate(_r[1][1])]

    elif mode == 'full':
        # There is nothing to crop in 'full' mode
        pass

    elif mode == 'valid':
        raise NotImplementedError("specklepy.core.alignment.pad_array does not support the 'valid' mode yet!")

    return padded


def _adapt_max_coordinate(index):
    """Cast the upper interval border, such that indexing returns the correct
    entries.

    Args:
        index (int):
            Index.

    Returns:
        None or negative index (NoneType or int):
            ... depending on index.

    """
    if index == 0:
        return None
    else:
        return - index
