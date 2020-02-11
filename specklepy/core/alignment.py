import numpy as np
from astropy.io import fits

from specklepy.logging import logging
from specklepy.utils.plot import imshow



def get_shifts(files, reference_file=0, mode='correlation', lazy_mode=True, return_image_shape=False, debug=False):
    """Computes the the relative shift of data cubes relative to a reference
    image.

    This function iterates over a list of files and uses the module function
    get_shift in 'correlation' mode to compute the relative shifts of files with
    respect to a reference file.

    Args:
        files (list):
            List of files to align.
        reference_file (str, int, optional):
            Path to a reference file or index of the file in files, relative to
            which the shifts are computed. Default is 0.
        mode (str, optional):
            Mode of the shift estimate. In 'correlation' mode, a 2D correlation
            is used to estimate the shift of the array. This is computationally
            much more expensive than the identical 'maximum' or 'peak' modes,
            which simply identify the coordinates of the emission peaks and
            return the difference. Though these modes may be fooled by reference
            sources of similar brightness. Passed to get_shift() function.
            Default is 'correlation'.
        lazy_mode (bool, optional):
            Set to False, to enforce the alignment of a single file with respect
            to the reference file. Default is True.
        return_image_shape (bool, optional):
            Set to True for for returning the shape of the anticipated output
            image. Default is False.
        debug (bool, optional):
            If set to True, it shows the 2D correlation.

    Returns:
        shifts (list):
            List of shifts for each file relative to the reference file.
    """

    # Check input parameters
    # TODO: Add a check, whether files is indeed a string before putting it into a list
    if not isinstance(files, list):
        files = [files]

    if isinstance(reference_file, int):
        reference_file = files[reference_file]
    elif reference_file is None:
        reference_file = files[0]
    elif not isinstance(reference_file, str):
        raise TypeError("The function get_shifts received reference_file argument of type {}, but needs be int or str, i.e. a file name.".format(type(reference_file)))

    # Skip computations if only one file is provided
    if lazy_mode and len(files) == 1:
        logging.info("Only one data cube is provided, nothing to align.")
        shifts = [(0, 0)]
        image_shape = fits.getdata(files[0]).shape
        image_shape = (image_shape[-2], image_shape[-1])

    # Otherwise estimate shifts
    else:
        shifts = []

        # Identify reference file and Fourier transform the integrated image
        logging.info(f"Computing relative shifts between data cubes. Reference file is {reference_file}")
        reference_image = fits.getdata(reference_file)
        if reference_image.ndim == 3:
            # Integrating over time axis if reference image is a cube
            reference_image = np.sum(reference_image, axis=0)
        Freference_image = np.fft.fft2(reference_image)
        image_shape = reference_image.shape
        del reference_image

        # Iterate over files and estimate shift via 2D correlation of the integrated cubes
        for index, file in enumerate(files):
            if file == reference_file:
                shift = (0, 0)
            else:
                image = fits.getdata(file)
                if image.ndim == 3:
                    image = np.sum(image, axis=0)
                shift = get_shift(image, reference_image=Freference_image, is_Fourier_transformed=True, mode=mode, debug=debug)
            shifts.append(shift)
            logging.info(f"Identified a shift of {shift} for file {file}")
        logging.info(f"Identified the following shifts:\n\t{shifts}")

    if return_image_shape:
        return shifts, image_shape
    else:
        return shifts



def get_shift(image, reference_image=None, is_Fourier_transformed=False, mode='correlation', debug=False):
    """Estimate the shift between an image and a reference image.

    Estimate the relative shift between an image and a reference image by means
    of a 2D correlation ('correlation' mode) or by comparison of the emission
    peaks ('peak' or 'maximum' modes).

    Args:
        image (np.ndarray):
            2D array of the image to be shifted.
        reference_image (np.ndarray):
            2D array of the reference image of the shift.
        is_Fourier_transformed (bool):
            Indicate whether the reference image is already Fourier transformed.
            This is implemented to save computation by computing that transform
            only once.
        mode (str, optional):
            Mode of the shift estimate. In 'correlation' mode, a 2D correlation
            is used to estimate the shift of the array. This is computationally
            much more expensive than the identical 'maximum' or 'peak' modes,
            which simply identify the coordinates of the emission peaks and
            return the difference. Though these modes may be fooled by reference
            sources of similar brightness. Default is 'correlation'.
        debug (bool, optional):
            Set to True to inspect intermediate results. Default is False.

    Returns:
        shift (tuple):
            Tuple of shift indizes for each axis.
    """

    # Check input parameters
    if not isinstance(image, np.ndarray) or image.ndim is not 2:
        raise TypeError(f"Image input must be 2D numpy.ndarray, but was provided as {type(image)}")
    if not isinstance(reference_image, np.ndarray) or image.ndim is not 2:
        raise TypeError(f"Image input must be 2D numpy.ndarray, but was provided as {type(reference_image)}")

    # Simple comparison of the peaks in the images
    if mode == 'maximum' or mode == 'peak':
        peak_image = np.unravel_index(np.argmax(image, axis=None), image.shape)
        peak_ref_image = np.unravel_index(np.argmax(reference_image, axis=None), reference_image.shape)
        return (peak_ref_image[0] - peak_image[0] , peak_ref_image[1] - peak_image[1])

    # Using correlation of the two images
    elif mode == 'correlation':
        if not is_Fourier_transformed:
            Freference_image = np.fft.fft2(reference_image)
        else:
            Freference_image = reference_image
        # if reference_image is not None and Freference_image is None:
        #     Freference_image = np.fft.fft2(reference_image)
        # elif Freference_image is not None and reference_image is None:
        #     pass
        # else:
        #     raise ValueError("Exactly one of reference_image or Freference_image needs be provided to get_shift!")
        Fimage = np.conjugate(np.fft.fft2(image))
        correlation = np.fft.ifft2(np.multiply(Freference_image, Fimage))
        correlation = np.fft.fftshift(correlation)
        if debug:
            imshow(np.abs(correlation), title='FFT shifted correlation of file {}'.format(index))
        shift = np.unravel_index(np.argmax(correlation), correlation.shape)
        shift = tuple(x - int(correlation.shape[i] / 2) for i, x in enumerate(shift))
        return shift

    else:
        raise ValueError("get_shift received unknown mode {}".format(mode))



def get_pad_vectors(shifts, cube_mode=False, return_reference_image_pad_vector=False):
    """Computes padding vectors from the relative shifts between files.

    Args:
        shifts (list):
            Shifts between files, relative to a reference image. See get_shifts
            function for details.
        cube_mode (bool, optional):
            If image is a cube, the estimated pad vectors will obtain pad_vector
            entries of (0, 0) for the zeroth axis. Default is False.
        return_reference_image_pad_vector (bool, optional):
            In 'same' mode, pad_array needs also the pad vector of the reference
            image. This is returned if this arg is set True. Default is False.

    Returns:
        pad_vectors (list):
            List of padding vectors for each shift in shifts.
        reference_image_pad_vector (list, optional):
            Pad vector of the reference image, which is needed for pad_array()
            in 'same' mode, thus is only returned in 'same' mode.
    """

    # Check input types


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

        pad_vector.append( (shift[0] - xmin, xmax - shift[0]) )
        pad_vector.append( (shift[1] - ymin, ymax - shift[1]) )

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

    Pads an array with zeros to match a desired field size. Intermediately, it
    always creates a 'full' image and only in 'same' mode it crops the edges
    such that the returned array covers only the field of the reference image.

    Args:
        array (np.ndarray):
            Input array that shall be padded to match the 'full' or 'same'
            fields.
        pad_vector (list):
            List of padding vectors, as obtained from get_pad_vectors().
        mode (str, optional):
            Define the size of the output image as 'same' to the reference image
            or expanding to include the 'full' covered field.

    Returns:
        padded (np.ndarray):
            Padded array, matching the field of the reference image in 'same'
            mode, or the complete field in 'full' mode.
    """

    if not isinstance(array, np.ndarray):
        raise TypeError("specklepy.core.alignment.pad_array received array argument \
                            of type {}, but must be np.ndarray.".format(type(array)))
    if array.ndim not in [2, 3]:
        raise ValueError("specklepy.core.alignment.pad_array received array argument \
                            of dimension {}, but must be 2 or 3.".format(array.ndim))

    padded = np.pad(array, pad_vector, mode='constant')

    # Crop the image according to the desired mode
    if mode is 'same':
        # Take reference pad vector and adapt to correctly limit the image
        _r = reference_image_pad_vector
        # Pick only those pixels, covered by the reference image
        if array.ndim == 2:
            padded = padded[_r[0][0] : _adapt_max_coordinate(_r[0][1]) , _r[1][0] : _adapt_max_coordinate(_r[1][1])]
        else:
            padded = padded[: , _r[0][0] : _adapt_max_coordinate(_r[0][1]) , _r[1][0] : _adapt_max_coordinate(_r[1][1])]

    elif mode is 'full':
        # There is nothing to crop in 'full' mode
        pass

    elif mode is 'valid':
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
