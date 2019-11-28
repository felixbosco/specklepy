import numpy as np
from astropy.io import fits

from holopy.logging import logging


def compute_shifts(files, reference_file=None, reference_file_index=0, debug=False):
    """Align the data cubes relative to a reference image.

    Long description...

    Args:
        files (list):
        reference_file (str, optional):
        reference_file_index (str, optional):
        debug (bool, optional):

    Returns:
        outfile_shape (tuple):
        pad_vectors (sequence):
    """

    files = np.atleast_1d(files)
    shifts = []

    # Skip computations if only one data cube is provided
    if len(files) == 1:
        logging.info("Only one data cube is provided, nothing to align.")
        shifts = [(0, 0)]
        pad_vectors = [(0, 0)]

    # Otherwise estimate shifts and compute pad vectors
    else:
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
        del reference_image

        # Iterate over inFiles and estimate shift via 2D correlation of the integrated cubes
        for index, file in enumerate(files):
            if index == reference_file_index:
                shift = (0, 0)
            else:
                image = np.sum(fits.getdata(file), axis=0)
                Fimage = np.conjugate(np.fft.fft2(image))
                correlation = np.fft.ifft2(np.multiply(Freference_image, Fimage))
                correlation = np.fft.fftshift(correlation)
                if debug:
                    imshow(np.abs(correlation), title='FFT shifted correlation of file {}'.format(index))
                shift = np.unravel_index(np.argmax(correlation), correlation.shape)
                shift = tuple(x - int(correlation.shape[i] / 2) for i, x in enumerate(shift))
                shift = tuple(-1 * i for i in shift)
            shifts.append(shift)
        logging.info("Identified the following shifts:\n\t{}".format(shifts))
    return shifts


def compute_pad_vectors(shifts, mode='same'):
    """Align the data cubes relative to a reference image.

    Long description...

    Args:
        mode (str, optional): Define the size of the output image as 'same'
            to the reference image or expanding to include the 'full'
            covered field.
        reference_file (str, optional):
        reference_file_index (str, optional):
        debug (bool, optional):

    Returns:
        outfile_shape (tuple):
        pad_vectors (sequence):
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
