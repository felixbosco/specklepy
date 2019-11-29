import numpy as np
from astropy.io import fits

from holopy.logging import logging
from holopy.utils.plot import imshow



def get_shifts(files, reference_file=None, reference_file_index=0, lazy_mode=True, return_image_shape=False, debug=False):
    """Computes the the relative shift of data cubes relative to a reference
    image.

    Long description...

    Args:
        files (list): List of files to align.
        reference_file (str, optional): Path to a reference file, relative to
            which the shifts are computed. If not provided, the reference file
            index is used, see below. Default is None.
        reference_file_index (str, optional): Index of the file in the file
            list, relative to which the shifts are cpmputed. Default is 0.
        lazy_mode (bool, optional): If set to True and  Default is True.
        debug (bool, optional): If set to True, it shows the 2D correlation.

    Returns:
        shifts (list): List of shifts for each file relative to the reference
            file.
    """

    if not isinstance(files, list):
        files = [files]

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
        if reference_file is None:
            reference_file = files[reference_file_index]
        else:
            reference_file_index = None

        logging.info("Computing relative shifts between data cubes. Reference file is {}".format(reference_file))
        reference_image = fits.getdata(reference_file)
        if reference_image.ndim == 3:
            # Integrating over time axis if reference image is a cube
            reference_image = np.sum(reference_image, axis=0)
        Freference_image = np.fft.fft2(reference_image)
        image_shape = reference_image.shape
        del reference_image

        # Iterate over files and estimate shift via 2D correlation of the integrated cubes
        for index, file in enumerate(files):
            if index == reference_file_index:
                shift = (0, 0)
            else:
                image = fits.getdata(file)
                if image.ndim == 3:
                    image = np.sum(image, axis=0)
                shift = estimate_shift(image, Freference_image=Freference_image, mode='correlation', debug=debug)
            shifts.append(shift)
        logging.info("Identified the following shifts:\n\t{}".format(shifts))

    if return_image_shape:
        return shifts, image_shape
    else:
        return shifts



def estimate_shift(image, reference_image=None, Freference_image=None, mode='maximum', debug=False):
    """Estimate the shift between an image and a refernce image.

    Long description ...

    Args:
        image (np.ndarray):
        reference_image (np.ndarray):
        Freference_image (np.ndarray):
        mode (str, optional):
        debug (bool):

    Returns:
        shift (tuple): Tuple of shift indizes for each axis.
    """

    # Simple comparison of the peaks in the images
    if mode == 'maximum' or mode == 'peak':
        peak_image = np.unravel_index(np.argmax(image, axis=None), image.shape)
        peak_ref_image = np.unravel_index(np.argmax(ref_image, axis=None), ref_image.shape)
        return (peak_image[0] - peak_ref_image[0], peak_image[1] - peak_ref_image[1])

    # Using correlation of the two images
    elif mode == 'correlation':
        if reference_image is not None and Freference_image is None:
            Freference_image = np.fft.fft2(reference_image)
        elif Freference_image is not None and reference_image is None:
            pass
        else:
            raise ValueError("Exactly one of reference_image or Freference_image needs be provided to estimate_shift!")
        Fimage = np.conjugate(np.fft.fft2(image))
        correlation = np.fft.ifft2(np.multiply(Freference_image, Fimage))
        correlation = np.fft.fftshift(correlation)
        if debug:
            imshow(np.abs(correlation), title='FFT shifted correlation of file {}'.format(index))
        shift = np.unravel_index(np.argmax(correlation), correlation.shape)
        shift = tuple(x - int(correlation.shape[i] / 2) for i, x in enumerate(shift))
        shift = tuple(-1 * i for i in shift)
        return shift

    else:
        raise ValueError("estimate_shift received unknown mode {}".format(mode))



def get_pad_vectors(shifts, image_shape, reference_image_shape, mode='same'):
    """Computes padding vectors from the relative shifts between files.

    Long description...

    Args:
        shifts (list): Shifts between files, relative to a reference image. See
            holopy.alignment.get_shifts for details.
        mode (str, optional): Define the size of the output image as 'same'
            to the reference image or expanding to include the 'full'
            covered field.
        reference_file (str, optional):
        reference_file_index (str, optional):
        debug (bool, optional):

    Returns:
        pad_vectors (list):
    """

    # Turn shifts into pad vectors
    if mode == 'same' or mode == 'full':
    #     pad_vectors = [len(files) * (0, 0)]
    # elif mode == 'full':
        xmax, ymax = np.max(-1 * np.array(shifts), axis=0)
        xmin, ymin = np.min(-1 * np.array(shifts), axis=0)
        shift_limits = {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax}
        # print('>>>>>>>>>', shift_limits)
        pad_vectors = []
        for shift in shifts:
            padding_x = (shift[0] - shift_limits['xmin'], shift_limits['xmax'] - shift[0])
            padding_y = (shift[1] - shift_limits['ymin'], shift_limits['ymax'] - shift[1])
            pad_vectors.append(((0, 0), padding_x, padding_y))
        for pad_vector in pad_vectors:
            print('>>>>>>>>>>>', pad_vector)
    else:
        raise ValueError("Mode '{}' not defined for {}.align_cubes()".format(mode))



def pad_files(files, shifts, mode='same', sum=False):
    """

    """
    pass
