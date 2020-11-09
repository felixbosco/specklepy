from IPython import embed
import numpy as np
import os

from astropy.io import fits

from specklepy.exceptions import SpecklepyTypeError, SpecklepyValueError
from specklepy.core import alignment
from specklepy.logging import logger


class SpeckleCube(object):

    def __init__(self, file_name, variance_extension=None, in_dir=None, out_dir=None):

        # Store input arguments
        self.in_dir, self.file_name = os.path.split(file_name)
        self.var_ext = variance_extension
        if in_dir:
            self.in_dir = in_dir
        self.out_dir = out_dir if out_dir else ''

        # Initialize attributes
        self.fits_header = None
        self.image = None
        self.image_var = None
        self.method = None
        self.box = None

    @property
    def in_file(self):
        return os.path.join(self.in_dir, self.file_name)

    def collapse(self):

        # Read data from file
        cube = fits.getdata(self.in_file)

        # Compute collapsed image from the cube
        self.image = np.sum(cube, axis=0)

        # Extract variance image
        if self.var_ext:
            image_var = fits.getdata(self.in_file, self.var_ext)
            if image_var.ndim == 2:
                self.image_var = image_var
            else:
                logger.warning('Image variance data is not an image!')

        # Set method attribute
        self.method = 'collapse'

    def ssa(self, box=None):

        if box:
            self.box = box

        # Read data from file
        cube = fits.getdata(self.in_file)
        if self.var_ext:
            var_cube = fits.getdata(self.in_file, self.var_ext)
        else:
            var_cube = None

        # Compute SSA'ed image from the cube
        self.image, self.image_var = coadd_frames(cube=cube, var_cube=var_cube, box=self.box)
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

    def store(self, path=None):
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
        if self.image_var:
            var_hdu = fits.ImageHDU(data=self.image_var, name=self.var_ext)
            hdu_list.append(var_hdu)

        # Store data to file
        hdu_list.writeto(self.default_save_path(), overwrite=True)

        # Return target path
        return self.default_save_path()


def coadd_frames(cube, var_cube=None, box=None):
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

    Returns:
        coadded (np.ndarray, ndim=2):
            SSA-integrated frames of the input cube.
        var_coadded (np.ndarray, ndim=2):
            SSA-integrated variances of the input cube or the variance map itself if provided as a 2D cube.
    """

    if not isinstance(cube, np.ndarray):
        raise SpecklepyTypeError('coadd_frames()', argname='cube', argtype=type(cube), expected='np.ndarray')
    if cube.ndim not in [2, 3]:
        raise SpecklepyValueError('coadd_frames()', argname='cube.ndim', argvalue=cube.ndim, expected='2 or 3')

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
    peak_indizes = np.zeros((cube.shape[0], 2), dtype=int)
    for index, frame in enumerate(cube):
        if box is not None:
            frame = box(frame)
        peak_indizes[index] = np.array(np.unravel_index(np.argmax(frame, axis=None), frame.shape), dtype=int)

    # Compute shifts from indizes
    peak_indizes = peak_indizes.transpose()
    xmean, ymean = np.mean(np.array(peak_indizes), axis=1)
    xmean = int(xmean)
    ymean = int(ymean)
    shifts = np.array([xmean - peak_indizes[0], ymean - peak_indizes[1]])
    shifts = shifts.transpose()

    # Shift frames and add to coadded
    coadded = np.zeros(cube[0].shape)
    pad_vectors, ref_pad_vector = alignment.get_pad_vectors(shifts, cube_mode=False,
                                                            return_reference_image_pad_vector=True)
    for index, frame in enumerate(cube):
        coadded += alignment.pad_array(frame, pad_vectors[index], mode='same',
                                       reference_image_pad_vector=ref_pad_vector)

    # Coadd variance cube (if not an image itself)
    if var_cube is not None:
        if var_cube.ndim == 3:
            var_coadded = np.zeros(coadded.shape)
            for index, frame in enumerate(var_cube):
                var_coadded += alignment.pad_array(frame, pad_vectors[index], mode='same',
                                                   reference_image_pad_vector=ref_pad_vector)
        elif var_cube.ndim == 2:
            var_coadded = var_cube
        else:
            raise RuntimeError(f"var_cube has unexpected shape: {var_cube.shape}")
    else:
        var_coadded = None

    return coadded, var_coadded
