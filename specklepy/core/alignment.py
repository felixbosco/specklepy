import numpy as np
import os

from specklepy.core.sourceextraction import SourceExtractor
from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.io.fits import get_data, get_frame_shape
from specklepy.logging import logger
from specklepy.plotting.utils import imshow
from specklepy.utils.array import peak_index


class FrameAlignment(object):

    def __init__(self, reference_file=None, reference_image=None, file_path=None):

        # Store reference image
        self.reference_file = reference_file
        self._reference_image = reference_image
        self.file_path = file_path

        # Initialize future attributes
        self.shifts = None
        self.pad_vectors = None
        self.reference_image_pad_vector = None

    def __repr__(self):
        return f"ShiftEstimator ({self.reference_image_path})"

    @property
    def reference_image_path(self):
        if self.file_path is None:
            return self.reference_file
        else:
            return os.path.join(self.file_path, self.reference_file)

    @property
    def reference_image(self):
        if self._reference_image is None:
            self._reference_image = self.load_collapsed(self.reference_file, path=self.file_path)
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

    def estimate_shifts(self, file_names, mode='correlation', in_dir=None, reference_file_index=None,
                        source_extractor_kwargs=None, debug=False):
        """

        Args:
            file_names:
            mode:
            in_dir:
            reference_file_index:
            source_extractor_kwargs (dict, optional):
                Used only if `mode=='sources'`. Parsed to the `SourceExtractor` instance.
            debug (bool, optional):
                Switch to `debug` mode.

        Returns:
            shifts (list):
                List of tuples of shifts between the images in the files of `file_names`.
        """

        # Transform str-type file names into a list
        if isinstance(file_names, str):
            file_names = [file_names]

        # Update reference file and input directory if provided
        if reference_file_index is not None:
            self.reference_file = file_names[reference_file_index]
        if in_dir is not None:
            self.file_path = in_dir

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
                    image = self.load_collapsed(file, path=self.file_path)
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
                    image = self.load_collapsed(file, path=self.file_path)
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
            if source_extractor_kwargs is None:
                source_extractor_kwargs = {'fwhm': 10}
            if 'fwhm' not in source_extractor_kwargs:
                source_extractor_kwargs['fwhm'] = 10

            # Initialize SourceExtractor
            extractor = SourceExtractor(**source_extractor_kwargs)
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
                    image = self.load_collapsed(file, path=self.file_path)
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

    def derive_pad_vectors(self, shifts=None):
        """Computes padding vectors from the relative shifts between files.

        Args:
            shifts (list or np.ndarray):
                Shifts between files, relative to a reference image. See get_shifts function for details.

        Returns:
            pad_vectors (list):
                List of padding vectors for each shift in shifts.
            reference_image_pad_vector (list, optional):
                Pad vector of the reference image, which is needed for pad_array() in 'same' mode.
        """

        # Update shifts if provided
        if shifts is not None:
            self.shifts = shifts

        # Initialize list
        self.pad_vectors = []

        # Get extreme points for 'full' padding
        y_max, x_max = np.max(np.array(self.shifts), axis=0)
        y_min, x_min = np.min(np.array(self.shifts), axis=0)

        # Iterate over Shifts
        for shift in self.shifts:

            pad_vector = [(shift[0] - y_min, y_max - shift[0]), (shift[1] - x_min, x_max - shift[1])]
            self.pad_vectors.append(pad_vector)

        # In 'same' mode, pad_array needs also the pad vector of the reference image
        self.reference_image_pad_vector = [(np.abs(y_min), np.abs(y_max)), (np.abs(x_min), np.abs(x_max))]

        return self.pad_vectors, self.reference_image_pad_vector

    def pad_array(self, array, pad_vector_index, mode='same', reference_image_pad_vector=None):
        """Pads an array according to the pad_vector and crops the image given the mode.

        Pads an array with zeros to match a desired field size. Intermediately, it always creates a 'full' image and
        only in 'same' mode it crops the edges such that the returned array covers only the field of the reference
        image.

        Args:
            array (np.ndarray):
                Input array that shall be padded to match the 'full' or 'same' fields.
            pad_vector_index (int):
                Index of the padding vector to apply from the stored list of pad vectors.
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

        # Update reference pad vector if provided
        if reference_image_pad_vector is not None:
            self.reference_image_pad_vector = reference_image_pad_vector

        # Pad the array to full size
        padded = np.pad(array, self.pad_vectors[pad_vector_index], mode='constant')

        # Crop the image according to the desired mode
        if mode == 'same':
            # Take reference pad vector and adapt to correctly limit the image
            cv = self.get_crop_vector()
            # Pick only those pixels, covered by the reference image
            if array.ndim == 2:
                padded = padded[cv[0][0]: cv[0][1], cv[1][0]:cv[1][1]]
            else:
                padded = padded[:, cv[0][0]: cv[0][1], cv[1][0]: cv[1][1]]

        elif mode == 'full':
            # There is nothing to crop in 'full' mode
            pass

        elif mode == 'valid':
            raise NotImplementedError("Alignment.pad_array does not support the 'valid' mode yet!")

        return padded

    def get_crop_vector(self):
        crop_vector = []
        for i in range(2):
            crop_vector_tuple = [self.reference_image_pad_vector[i][0]]
            max_index = self.reference_image_pad_vector[i][1]
            if max_index == 0:
                crop_vector_tuple.append(None)
            else:
                crop_vector_tuple.append(- max_index)
            crop_vector.append(crop_vector_tuple)
        return crop_vector
