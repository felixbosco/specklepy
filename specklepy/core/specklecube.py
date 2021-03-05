from IPython import embed
import numpy as np
import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.core import alignment
from specklepy.logging import logger


class SpeckleCube(object):

    def __init__(self, file_name, variance_extension=None, mask_extension='MASK', in_dir=None, out_dir=None):

        # Store input arguments
        self.in_dir, self.file_name = os.path.split(file_name)
        self.var_ext = variance_extension
        self.mask_ext = mask_extension
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
            self._cube = fits.getdata(self.in_file).squeeze()
        return self._cube

    @property
    def variance(self):
        if self.var_ext is None:
            return None
        try:
            return fits.getdata(self.in_file, self.var_ext).squeeze()
        except:
            return None

    @property
    def mask(self):
        if self.mask_ext is None:
            return None
        try:
            return fits.getdata(self.in_file, self.mask_ext)
        except KeyError:
            return None

    def collapse(self):

        # Compute collapsed image from the cube
        if self.cube.ndim == 2:
            self.image = self.cube
        elif self.cube.ndim == 3:
            self.image = np.sum(self.cube, axis=0)
        else:
            raise NotImplementedError(f"Frame alignment is not defined for FITS cubes with {self.cube.ndim!r} axes!")

        # Extract variance image
        if self.var_ext:
            if self.variance.ndim == 2:
                self.image_var = self.variance
            else:
                logger.warning('Image variance data is not an image!')

        # Set method attribute
        self.method = 'collapse'

    def ssa(self, box=None):

        # Update box attribute
        if box:
            self.box = box

        # Compute SSA'ed image from the cube
        if self.cube.ndim == 2:
            self.image = self.cube
            self.image_var = self.variance
        elif self.cube.ndim == 3:
            self.image, self.image_var = coadd_frames(cube=self.cube, var_cube=self.variance, box=self.box,
                                                      mask=self.mask)
        else:
            raise NotImplementedError(f"Frame alignment is not defined for FITS cubes with {self.cube.ndim!r} axes!")
        self.method = 'ssa'

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
        # TODO: implement using an OutFile instance here.

        # Stop if image is not computed yet
        if self.image is None:
            raise ValueError

        # Update out_dir path if new provided
        if path:
            self.out_dir = path

        # Set up FITS HDU list
        self.fits_header = self.make_header()
        hdu = fits.PrimaryHDU(data=self.image, header=self.fits_header)
        hdu_list = fits.HDUList(hdus=[hdu])

        # Append variance extension if available
        if self.image_var is not None:
            var_hdu = fits.ImageHDU(data=self.image_var, name=self.var_ext)
            hdu_list.append(var_hdu)

        # Store data to file
        try:
            hdu_list.writeto(self.default_save_path(), overwrite=True)
        except fits.verify.VerifyError as e:
            # Something is weird in the FITS header, so reset
            logger.error(e)
            logger.error("Amend this and try again!")
            embed()
            hdu_list.writeto(self.default_save_path(), overwrite=True)

        # Return target path
        return self.default_save_path()


def coadd_frames(cube, var_cube=None, box=None, mask=None):
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
        elif var_cube.ndim == cube.ndim-1:
            if var_cube.shape[0] != cube.shape[1] or var_cube.shape[1] != cube.shape[2]:
                raise SpecklepyValueError('coadd_frames()', argname='var_cube.shape', argvalue=str(var_cube.shape),
                                          expected=str(cube.shape))

    # Compute shifts
    peak_indexes = np.zeros((cube.shape[0], 2), dtype=int)
    for f, frame in enumerate(cube):
        if box is not None:
            frame = box(frame)
        if mask is not None:
            frame = np.ma.masked_array(frame, mask=mask).filled(0)
        peak_indexes[f] = np.array(alignment.peak_index(frame), dtype=int)

    # Compute shifts from indexes
    peak_indexes = peak_indexes.transpose()
    mean_index = np.mean(np.array(peak_indexes), axis=1)
    shifts = np.array([int(mean_index[0]) - peak_indexes[0], int(mean_index[1]) - peak_indexes[1]])
    shifts = shifts.transpose()

    # Shift frames and add to coadded
    coadded = np.zeros(cube[0].shape)
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
