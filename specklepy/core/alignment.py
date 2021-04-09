import numpy as np
import os

from astropy.io import fits

from specklepy.core.sourceextraction import SourceExtractor
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
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
    def load_collapsed(path):
        """Load and collapse a data cube by integrating over time axis if not 2D image

        Arguments:
            path (str):
                Path to the file.

        Returns:
            data (np.ndarray):
                2D array, obtained from the input file.
        """

        logger.debug(f"Loading data from file {path!r}")
        data = fits.getdata(path).squeeze()

        if data.ndim == 3:
            data = np.sum(data, axis=0)

        return data

    @staticmethod
    def make_path(file_name, in_dir=None):
        try:
            return os.path.join(in_dir, file_name)
        except TypeError:
            return file_name

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
                    path = self.make_path(file, in_dir=self.in_dir)
                    image = self.load_collapsed(path=path)
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
                    path = self.make_path(file, in_dir=self.in_dir)
                    image = self.load_collapsed(path=path)
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
                    path = self.make_path(file, in_dir=self.in_dir)
                    image = self.load_collapsed(path=path)
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


def estimate_shifts(files, reference_file=None, mode='correlation', lazy_mode=True, return_image_shape=False,
                    in_dir=None, debug=False):
    """Computes the the relative shift of data cubes relative to a reference
    image.

    This function iterates over a list of files and uses the module function get_shift in 'correlation' mode to compute
    the relative shifts of files with respect to a reference file.

    Args:
        files (list or array_like):
            List of files to align.
        reference_file (str, int, optional):
            Path to a reference file or index of the file in files, relative to which the shifts are computed. Default
            is 0.
        mode (str, optional):
            Mode of the shift estimate. In 'correlation' mode, a 2D correlation is used to estimate the shift of the
            array. This is computationally much more expensive than the identical 'maximum' or 'peak' modes, which
            simply identify the coordinates of the emission peaks and return the difference. Though these modes may be
            fooled by reference sources of similar brightness. Passed to get_shift() function. Default is 'correlation'.
        lazy_mode (bool, optional):
            Set to False, to enforce the alignment of a single file with respect to the reference file. Default is True.
        return_image_shape (bool, optional):
            Set to True for for returning the shape of the anticipated output image. Default is False.
        in_dir (str, optional):
            Path to the files. `None` is substituted by an empty string.
        debug (bool, optional):
            If set to True, it shows the 2D correlation.

    Returns:
        shifts (list):
            List of shifts for each file relative to the reference file.
    """

    # Check input parameters
    if isinstance(files, str):
        files = [files]
    if not isinstance(files, (list, np.ndarray)):
        raise SpecklepyTypeError('get_shifts()', argname='files', argtype=type(files), expected='list')

    if reference_file is None:
        reference_file = files[0]
    elif isinstance(reference_file, int):
        reference_file = files[reference_file]
    elif not isinstance(reference_file, str):
        raise SpecklepyTypeError('get_shifts()', argname='reference_file', argtype=type(reference_file), expected='str')

    if not isinstance(mode, str):
        raise SpecklepyTypeError('get_shifts()', argname='mode', argtype=type(mode), expected='str')

    if not isinstance(lazy_mode, bool):
        raise SpecklepyTypeError('get_shifts()', argname='lazy_mode', argtype=type(lazy_mode), expected='bool')

    if not isinstance(return_image_shape, bool):
        raise SpecklepyTypeError('get_shifts()', argname='return_image_shape', argtype=type(return_image_shape),
                                 expected='bool')

    if in_dir is None:
        in_dir = ''

    # Skip computations if only one file is provided
    if lazy_mode and len(files) == 1:
        logger.info("Only one data cube is provided, nothing to align.")
        shifts = [(0, 0)]
        image_shape = fits.getdata(os.path.join(in_dir, files[0])).shape
        image_shape = (image_shape[-2], image_shape[-1])

    # Otherwise estimate shifts
    else:
        shifts = []

        # Identify reference file and Fourier transform the integrated image
        logger.info(f"Computing relative shifts between data cubes. Reference file is {reference_file!r}")
        reference_image = fits.getdata(os.path.join(in_dir, reference_file)).squeeze()
        image_shape = reference_image.shape
        logger.debug(f"Reference image has shape {image_shape}")

        # Collapse cube by integrating over time axis if reference image is a cube
        if reference_image.ndim == 3:
            reference_image = np.sum(reference_image, axis=0)

        # Fourier transform image for `correlation` mode
        if mode == 'correlation':
            reference_image = np.fft.fft2(reference_image)
            ref_image_ft = True
        else:
            ref_image_ft = False

        # Iterate through files and estimate shift via 2D correlation of the integrated cubes
        for index, file in enumerate(files):
            if file == reference_file:
                shift = (0, 0)
            else:
                image = fits.getdata(os.path.join(in_dir, file)).squeeze()
                if image.ndim == 3:
                    image = np.sum(image, axis=0)
                shift = estimate_shift(image, reference_image=reference_image, is_fourier_transformed=ref_image_ft,
                                       mode=mode, debug=debug)

            # Store estimate
            logger.info(f"Estimated a shift of {shift} for file {file!r}")
            shifts.append(shift)

    # Return
    if return_image_shape:
        return shifts, image_shape
    else:
        return shifts


def estimate_shift(image, reference_image, mode='correlation', is_fourier_transformed=False, debug=False):
    """Estimate the shift between an image and a reference image.

    Estimate the relative shift between an image and a reference image by means of a 2D correlation
    ('correlation' mode) or by comparison of the emission peaks ('peak' or 'maximum' modes).

    Args:
        image (np.ndarray):
            2D array of the image to be shifted.
        reference_image (np.ndarray):
            2D array of the reference image of the shift.
        is_fourier_transformed (bool, optional):
            Indicate whether the reference image is already Fourier transformed. This is implemented to save
            computation by computing that transform only once.
        mode (str, optional):
            Mode of the shift estimate. In 'correlation' mode, a 2D correlation is used to estimate the shift of the
            array. This is computationally much more expensive than the identical 'maximum' or 'peak' modes, which
            simply identify the coordinates of the emission peaks and return the difference. Though these modes may be
            fooled by reference sources of similar brightness. Default is 'correlation'.
        debug (bool, optional):
            Set to True to inspect intermediate results. Default is False.

    Returns:
        shift (tuple):
            Tuple of shift indices for each axis.
    """

    # Check input parameters
    if not isinstance(image, np.ndarray) or image.ndim != 2:
        raise TypeError(f"Image input must be 2-dim numpy.ndarray, but was provided as {type(image)}")
    if not isinstance(reference_image, np.ndarray) or image.ndim != 2:
        raise TypeError(f"Image input must be 2-dim numpy.ndarray, but was provided as {type(reference_image)}")
    if not isinstance(is_fourier_transformed, bool):
        raise SpecklepyTypeError('estimate_shift()', argname='is_fourier_transformed',
                                 argtype=type(is_fourier_transformed), expected='bool')
    if not isinstance(mode, str):
        raise SpecklepyTypeError('estimate_shift()', argname='mode', argtype=type(mode), expected='str')

    # Simple comparison of the peaks in the images
    if mode == 'maximum' or mode == 'peak':
        peak_image = peak_index(image)
        peak_ref_image = peak_index(reference_image)
        return peak_ref_image[0] - peak_image[0], peak_ref_image[1] - peak_image[1]

    # Using correlation of the two images
    elif mode == 'correlation':
        # Get the Fourier transformed reference image for cross-correlation
        if not is_fourier_transformed:
            f_reference_image = np.fft.fft2(reference_image)
        else:
            f_reference_image = reference_image

        # Fourier transform the image
        f_image = np.conjugate(np.fft.fft2(image))

        # Compute the 2-dimensional correlation
        correlation = np.fft.ifft2(np.multiply(f_reference_image, f_image))
        correlation = np.fft.fftshift(correlation)
        if debug:
            imshow(np.abs(correlation), title='FFT shifted correlation')

        # Derive the shift from the correlation
        shift = np.unravel_index(np.argmax(correlation), correlation.shape)
        shift = tuple(x - int(correlation.shape[i] / 2) for i, x in enumerate(shift))
        return shift

    else:
        raise NotImplementedError(f"Estimating shift in {mode!r} mode is not implemented!")


def derive_pad_vectors(shifts, cube_mode=False, return_reference_image_pad_vector=False):
    """Computes padding vectors from the relative shifts between files.

    Args:
        shifts (list or np.ndarray):
            Shifts between files, relative to a reference image. See get_shifts function for details.
        cube_mode (bool, optional):
            If image is a cube, the estimated pad vectors will obtain pad_vector entries of (0, 0) for the zeroth axis.
            Default is False.
        return_reference_image_pad_vector (bool, optional):
            In 'same' mode, pad_array needs also the pad vector of the reference image. This is returned if this arg is
            set True. Default is False.

    Returns:
        pad_vectors (list):
            List of padding vectors for each shift in shifts.
        reference_image_pad_vector (list, optional):
            Pad vector of the reference image, which is needed for pad_array()
            in 'same' mode, thus is only returned in 'same' mode.
    """

    # Check input parameters
    if isinstance(shifts, (list, np.ndarray)):
        pass
    else:
        raise SpecklepyTypeError('get_pad_vectors()', argname='shifts', argtype=type(shifts), expected='list')
    if not isinstance(cube_mode, bool):
        raise SpecklepyTypeError('get_pad_vectors()', argname='cube_mode', argtype=type(cube_mode), expected='bool')
    if not isinstance(return_reference_image_pad_vector, bool):
        raise SpecklepyTypeError('get_pad_vectors()', argname='return_reference_image_pad_vector',
                                 argtype=type(return_reference_image_pad_vector), expected='bool')

    # Initialize list
    pad_vectors = []

    # Get extreme points for 'full' padding
    xmax, ymax = np.max(np.array(shifts), axis=0)
    xmin, ymin = np.min(np.array(shifts), axis=0)

    # Iterate over Shifts
    for shift in shifts:

        if cube_mode:
            pad_vector = [(0, 0)]
        else:
            pad_vector = []

        pad_vector.append((shift[0] - xmin, xmax - shift[0]))
        pad_vector.append((shift[1] - ymin, ymax - shift[1]))

        pad_vectors.append(pad_vector)

    # In 'same' mode, pad_array needs also the pad vector of the reference image
    if return_reference_image_pad_vector:
        reference_image_pad_vector = [(np.abs(xmin), np.abs(xmax)), (np.abs(ymin), np.abs(ymax))]
        return pad_vectors, reference_image_pad_vector
    else:
        return pad_vectors


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
