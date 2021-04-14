import numpy as np
import os

from specklepy.core import FrameAlignment
from specklepy.io import FileStream
from specklepy.io.fits import get_data, get_header
from specklepy.logging import logger
from specklepy.reduction.filter import bad_pixel_mask, fill_hot_pixels, mask_hot_pixels
from specklepy.utils.array import frame_number, frame_shape, peak_index
from specklepy.utils.bootstrap import random_draw_vectors


class SpeckleCube(object):

    variance_extension = 'VAR'
    mask_extension = 'MASK'
    fill_value = 0

    def __init__(self, file_name, variance_extension=None, mask_extension=None, in_dir=None, out_dir=None):

        # Store input arguments
        self.in_dir, self.file_name = os.path.split(file_name)
        if variance_extension is not None:
            self.variance_extension = variance_extension
        if mask_extension is not None:
            self.mask_extension = mask_extension
        if in_dir:
            self.in_dir = in_dir
        self.out_dir = out_dir if out_dir else ''

        # Initialize attributes
        self.fits_header = None
        self.image = None
        self.image_var = None
        self.method = None
        self.box = None
        self._cube = None

    @property
    def cube(self):
        if self._cube is None:
            self._cube = get_data(self.file_name, path=self.in_dir)
        return self._cube

    @property
    def variance(self):
        if self.variance_extension is None:
            return None
        return get_data(self.file_name, path=self.in_dir, extension=self.variance_extension,
                        ignore_missing_extension=True)

    @property
    def mask(self):
        if self.mask_extension is None:
            return None
        return get_data(self.file_name, path=self.in_dir, extension=self.mask_extension, dtype=bool,
                        ignore_missing_extension=True)

    @property
    def single_frame_mode(self):
        return frame_number(self.cube) == 1

    @property
    def frame_number(self):
        return frame_number(self.cube)

    @property
    def frame_shape(self):
        return frame_shape(self.cube)

    def fill_masked(self, mask, fill_value=None):

        # Update fill value
        if fill_value is not None:
            self.fill_value = fill_value

        # Fill masked pixels
        self.image[mask] = self.fill_value
        if self.image_var is not None:
            self.image_var[mask] = 0  # Variance is independent on the fill value

    def collapse(self, mask=False, fill_value=None):

        # Compute collapsed image from the cube
        if self.cube.ndim == 2:
            self.image = self.cube
        elif self.cube.ndim == 3:
            self.image = np.sum(self.cube, axis=0)
        else:
            raise NotImplementedError(f"Frame alignment is not defined for FITS cubes with {self.cube.ndim!r} axes!")

        # Extract variance image
        if self.variance_extension:
            if self.variance.ndim == 2:
                self.image_var = self.variance
            else:
                logger.warning('Image variance data is not an image!')

        # Fill masked pixels
        self.fill_masked(mask=mask, fill_value=fill_value)

        # Set method attribute
        self.method = 'collapse'

        return self.image, self.image_var

    def estimate_peak_indexes(self, mask=None):
        # Compute shifts
        peak_indexes = np.zeros((self.frame_number, 2), dtype=int)
        for f, frame in enumerate(self.cube):
            if mask is not None:
                frame[mask] = self.fill_value
            if self.box is not None:
                frame = self.box(frame)
            peak_indexes[f] = np.array(peak_index(frame), dtype=int)
        return peak_indexes

    def estimate_frame_shifts(self, mask=None):
        peak_indexes = self.estimate_peak_indexes(mask=mask).transpose()
        mean_index = np.mean(np.array(peak_indexes), axis=1)
        shifts = np.array([int(mean_index[0]) - peak_indexes[0], int(mean_index[1]) - peak_indexes[1]])
        return shifts.transpose()

    def align_frames(self, mask=None, bootstrap_number=None):
        """Use the SSA algorithm to align the frames of the data cube.

        Args:
            mask (np.ndarray, optional):
                Bad pixel mask, parsed down to the method that estimates the shifts between frames.
            bootstrap_number (int, optional):
                Number of bootstrap images to create for sampling the uncertainty of the aligned image.

        Returns:
            aligned (np.ndarray):
                Aligned image.
            variance (np.ndarray):
                Variance map of the aligned image. This is either the aligned variance of the frames or the bootstrap
                resampling variance of the image, if `bootstrap_number` is provided.
        """

        # Estimate shifts between frames
        shifts = self.estimate_frame_shifts(mask=mask)

        # Initialize output images
        aligned = np.zeros(self.frame_shape)
        aligned_var = np.zeros(self.frame_shape) if self.variance is not None else None
        bootstrap_var = np.zeros(self.frame_shape) if bootstrap_number is not None else None
        bootstrap_images = [np.zeros(self.frame_shape)] * bootstrap_number if bootstrap_number is not None else None
        bootstrap_vector = random_draw_vectors(number_draws=bootstrap_number, number_frames=self.frame_number)

        # Shift frames and add to `aligned`
        alignment = FrameAlignment()
        alignment.derive_pad_vectors(shifts=shifts)
        for f, frame in enumerate(self.cube):
            padded = alignment.pad_array(frame, pad_vector_index=f, mode='same')
            aligned += padded

            # Co-add variance
            if aligned_var is not None:
                aligned_var += alignment.pad_array(self.variance, pad_vector_index=f, mode='same')

            # Co-add bootstrap frames
            for b in range(bootstrap_number):
                bootstrap_coefficient = bootstrap_vector[b, f]
                if bootstrap_coefficient != 0:
                    bootstrap_images[b] += bootstrap_coefficient * padded

        # Evaluate variance
        if bootstrap_var is None:
            # Account for variance averaging
            variance = aligned_var / self.frame_number
        else:
            # Evaluate variance of bootstrap images
            variance = np.var(bootstrap_images, axis=0)

        return aligned, variance

    def ssa(self, box=None, mask_bad_pixels=False, mask=False, fill_value=None, bootstrap_number=None):

        # Update box attribute
        if box:
            self.box = box
        if fill_value is not None:
            self.fill_value = fill_value

        # Compute SSA'ed image from the cube
        if self.cube.ndim == 2:
            self.image = self.cube
            self.image_var = self.variance
            self.fill_masked(mask=mask, fill_value=fill_value)

        elif self.cube.ndim == 3:
            # Build bad pixel mask and combine masks
            bpm = bad_pixel_mask(cube=self.cube, var=self.variance) if mask_bad_pixels else False
            mask = np.logical_or(self.mask, bpm) | mask

            # Coadd frames
            self.image, self.image_var = self.align_frames(mask=mask, bootstrap_number=bootstrap_number)
        else:
            raise NotImplementedError(f"Frame alignment is not defined for FITS cubes with {self.cube.ndim!r} axes!")
        self.method = 'ssa'

        return self.image, self.image_var

    def mask_hot_pixels(self, fill_value=0):
        if fill_value is None:
            self.image = mask_hot_pixels(image=self.image)
        else:
            self.image = fill_hot_pixels(self.image, fill_value=fill_value)

    def default_save_path(self):
        if self.method == 'ssa':
            prefix = 'ssa_'
        elif self.method == 'collapse':
            prefix = 'int_'
        else:
            raise ValueError
        return os.path.join(self.out_dir, prefix + self.file_name)

    def make_header(self):
        return get_header(path=self.in_dir, file_name=self.file_name)

    def write(self, path=None):

        # Stop if image is not computed yet
        if self.image is None:
            raise ValueError("SpeckleCube instance contains no image, so nothing to write!")

        # Update out_dir path if new provided
        if path:
            self.out_dir = path

        # Setup file stream
        file_stream = FileStream(self.default_save_path())
        file_stream.initialize(data=self.image, header=self.make_header(), overwrite=True)
        if self.image_var is not None:
            file_stream.new_extension(data=self.image_var, name=self.variance_extension)

        # Return target path
        return self.default_save_path()
