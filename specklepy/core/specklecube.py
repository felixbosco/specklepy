import numpy as np
import os
import sys

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.core import alignment
from specklepy.io.filestream import FileStream
from specklepy.io.fits import get_data
from specklepy.logging import logger
from specklepy.reduction.filter import bad_pixel_mask, fill_hot_pixels, mask_hot_pixels
from specklepy.utils.array import frame_number, frame_shape, peak_index


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
    def in_file(self):
        return os.path.join(self.in_dir, self.file_name)

    @property
    def cube(self):
        if self._cube is None:
            self._cube = get_data(self.in_file)
        return self._cube

    @property
    def variance(self):
        if self.variance_extension is None:
            return None
        return get_data(self.in_file, extension=self.variance_extension, ignore_missing_extension=True)

    @property
    def mask(self):
        if self.mask_extension is None:
            return None
        return get_data(self.in_file, extension=self.mask_extension, dtype=bool, ignore_missing_extension=True)

    @property
    def single_frame_mode(self):
        return frame_number(self.cube) == 1

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

    def estimate_peak_indexes(self, cube, box=None, mask=None):
        # Compute shifts
        peak_indexes = np.zeros((frame_number(cube), 2), dtype=int)
        for f, frame in enumerate(cube):
            if mask is not None:
                frame[mask] = self.fill_value
            if box is not None:
                frame = box(frame)
            peak_indexes[f] = np.array(peak_index(frame), dtype=int)
        return peak_indexes

    def estimate_frame_shifts(self, cube, box=None, mask=None):
        peak_indexes = self.estimate_peak_indexes(cube=cube, box=box, mask=mask).transpose()
        mean_index = np.mean(np.array(peak_indexes), axis=1)
        shifts = np.array([int(mean_index[0]) - peak_indexes[0], int(mean_index[1]) - peak_indexes[1]])
        return shifts.transpose()

    def align_frames(self, mask=None):

        # Estimate shifts between frames
        shifts = self.estimate_frame_shifts(self.cube, box=self.box, mask=mask)

        # Initialize output images
        aligned = np.zeros(self.frame_shape)
        aligned_var = np.zeros(self.frame_shape) if self.variance is not None else None

        # Shift frames and add to `aligned`
        pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(shifts, cube_mode=False,
                                                                   return_reference_image_pad_vector=True)
        for f, frame in enumerate(self.cube):
            aligned += alignment.pad_array(frame, pad_vectors[f], mode='same',
                                           reference_image_pad_vector=ref_pad_vector)
            if aligned_var is not None:
                aligned_var += alignment.pad_array(self.variance, pad_vectors[f], mode='same',
                                                   reference_image_pad_vector=ref_pad_vector)
        return aligned, aligned_var

    def ssa(self, box=None, mask_bad_pixels=False, mask=False, fill_value=None):

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
            self.image, self.image_var = self.align_frames(mask=mask)
            # self.image, self.image_var = coadd_frames(cube=self.cube, var_cube=self.variance, box=self.box, mask=mask,
            #                                           fill_value=self.fill_value)
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
        hdr = fits.getheader(os.path.join(self.in_dir, self.file_name))
        return hdr

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

        # Set up FITS HDU list
        # self.fits_header = self.make_header()
        # hdu = fits.PrimaryHDU(data=self.image, header=self.fits_header)
        # hdu_list = fits.HDUList(hdus=[hdu])
        #
        # # Append variance extension if available
        # if self.image_var is not None:
        #     var_hdu = fits.ImageHDU(data=self.image_var, name=self.variance_extension)
        #     hdu_list.append(var_hdu)
        #
        # # Store data to file
        # try:
        #     hdu_list.writeto(self.default_save_path(), overwrite=True)
        # except fits.verify.VerifyError as e:
        #     # Something is weird in the FITS header, so reset
        #     logger.error(e)
        #     logger.error("Amend this and try again!")
        #     hdu_list.writeto(self.default_save_path(), overwrite=True)

        # Return target path
        return self.default_save_path()


def coadd_frames(cube, var_cube=None, box=None, mask=None, fill_value=0):
    """Compute the simple shift-and-add (SSA) reconstruction of a data cube.

    This function uses the SSA algorithm to coadd frames of a cube. If provided, this function coadds the variances
    within a var cube considering the exact same shifts.

    Args:
        cube (np.ndarray, ndim=3):
            Data cube which is integrated along the zero-th axis.
        var_cube (np.ndarray, ndim=3, optional):
            Data cube of variances which is integrated along the zero-th axis with the same shifts as the cube.
        box (Box object, optional):
            Constraining the search for the intensity peak to the specified box. Searching the full frames if not
            provided.
        mask (np.ndarray, optional):
            Bad pixel mask in the shape of a frame of the cube.
        fill_value (float-like, optional):
            Masked pixels will obtain this value prior to co-addition.

    Returns:
        coadded (np.ndarray, ndim=2):
            SSA-integrated frames of the input cube.
        var_coadded (np.ndarray, ndim=2):
            SSA-integrated variances of the input cube or the variance map itself if provided as a 2D cube.
    """

    # Assert that cube is a 3-dimensional np.ndarray
    try:
        if cube.ndim != 3:
            raise SpecklepyValueError('coadd_frames()', argname='cube.ndim', argvalue=cube.ndim, expected='2 or 3')
    except AttributeError:
        raise SpecklepyTypeError('coadd_frames()', argname='cube', argtype=type(cube), expected='np.ndarray')

    # Check on the shape of var_cube
    if var_cube is not None:
        if not isinstance(var_cube, np.ndarray):
            raise SpecklepyTypeError('coadd_frames()', argname='var_cube', argtype=type(var_cube),
                                     expected='np.ndarray')
        if var_cube.ndim == cube.ndim and var_cube.shape != cube.shape:
            raise SpecklepyValueError('coadd_frames()', argname='var_cube.shape', argvalue=str(var_cube.shape),
                                      expected=str(cube.shape))
        elif var_cube.ndim == cube.ndim-1 and frame_shape(var_cube) != frame_shape(cube):
            raise SpecklepyValueError('coadd_frames()', argname='var_cube.shape', argvalue=str(var_cube.shape),
                                      expected=str(cube.shape))

    # Compute shifts
    peak_indexes = np.zeros((frame_number(cube), 2), dtype=int)
    for f, frame in enumerate(cube):
        if mask is not None:
            frame[mask] = fill_value
        if box is not None:
            frame = box(frame)
        try:
            peak_indexes[f] = np.array(peak_index(frame), dtype=int)
        except ValueError as e:
            sys.tracebacklimit = 0
            logger.warning(f"The box {box} might be set inappropriate!")
            raise e

    # Compute shifts from indexes
    peak_indexes = peak_indexes.transpose()
    mean_index = np.mean(np.array(peak_indexes), axis=1)
    shifts = np.array([int(mean_index[0]) - peak_indexes[0], int(mean_index[1]) - peak_indexes[1]])
    shifts = shifts.transpose()

    # Shift frames and add to coadded
    coadded = np.zeros(frame_shape(cube))
    pad_vectors, ref_pad_vector = alignment.derive_pad_vectors(shifts, cube_mode=False,
                                                               return_reference_image_pad_vector=True)
    for f, frame in enumerate(cube):
        coadded += alignment.pad_array(frame, pad_vectors[f], mode='same', reference_image_pad_vector=ref_pad_vector)

    # Coadd variance cube (if not an image itself)
    if var_cube is not None:
        if var_cube.ndim == 3:
            var_coadded = np.zeros(coadded.shape)
            for f, frame in enumerate(var_cube):
                var_coadded += alignment.pad_array(frame, pad_vectors[f], mode='same',
                                                   reference_image_pad_vector=ref_pad_vector)
        elif var_cube.ndim == 2:
            var_coadded = var_cube
        else:
            raise RuntimeError(f"var_cube has unexpected shape: {var_cube.shape}")
    else:
        var_coadded = None

    return coadded, var_coadded
